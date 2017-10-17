/*
Copyright 2017 Andrew R. Booker

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

static uint64_t *pseudoprime_table;
static int32_t *pseudoprime_index;
static uint32_t *collision_table;
static int32_t *collision_index;
static int32_t *factor_index;
static uint16_t *factor_table;

int initfactor64(const char *db) {
	int fd;
	uint64_t length;
	void *map;

	if ((fd=open(db,O_RDONLY)) < 0) return -1;
	length = 31894014u*sizeof(*pseudoprime_table)+64u*sizeof(*pseudoprime_index)
		+203280220u*sizeof(*collision_table)+852034u*sizeof(*collision_index)
		+510510u*sizeof(*factor_index)+775350505u*sizeof(*factor_table);
	map = mmap((void *)0,length,PROT_READ,MAP_SHARED,fd,0);
	if (map == MAP_FAILED) return -2;
	madvise(map,length,MADV_WILLNEED);

	pseudoprime_table = (uint64_t *)map;
	pseudoprime_index = (int32_t *)&pseudoprime_table[31894014];
	collision_table = (uint32_t *)&pseudoprime_index[64];
	collision_index = (int32_t *)&collision_table[203280220];
	factor_index = (int32_t *)&collision_index[852034];
	factor_table = (uint16_t *)&factor_index[510510];

	return 0;
}

#if defined(__x86_64__) || defined(_M_X64)
// for x86 we optimize the low level arithmetic in assembly language

// |a-b|
static inline uint64_t absdiff(uint64_t a,uint64_t b) {
	uint64_t t = b-a;
	asm("sub %2, %0\n\t"    // a -= b;
			"cmovc %1, %0\n\t"  // if (carry) a = t;
			:"+&r" (a)
			:"r" (t), "r" (b)
			:"cc"
	);
	return a;
}

// addmod from Kruppa 2010 page 67
static inline uint64_t addmod(uint64_t a,uint64_t b,uint64_t n) {
	uint64_t t = a-n;
	a += b;
	asm("add %2, %1\n\t"    // t += b
			"cmovc %1, %0\n\t"  // if (carry) a = t
			:"+r" (a), "+&r" (t)
			:"r" (b)
			:"cc"
	);
	return a;
}

// Ben Buhrow's x86 mulredc functions
static inline uint64_t mulredc64(uint64_t a,uint64_t b,uint64_t n,uint64_t npi) {
	asm("mulq %1 \n\t"
			"movq %%rax, %%r10 \n\t"
			"movq %%rdx, %%r11 \n\t"
			"movq $0, %%r12 \n\t"
			"mulq %2 \n\t"
			"mulq %3 \n\t"
			"addq %%r10, %%rax \n\t"
			"adcq %%r11, %%rdx \n\t"
			"cmovae %3, %%r12 \n\t"
			"xorq %%rax, %%rax \n\t"
			"subq %3, %%rdx \n\t"
			"cmovc %%r12, %%rax \n\t"
			"addq %%rdx, %%rax \n\t"
			: "+&a"(a)
			: "r"(b), "r"(npi), "r"(n)
			: "rdx", "r10", "r11", "r12", "cc"
	);
  return a;
}

// assumes n < 2^63
static inline uint64_t mulredc63(uint64_t a,uint64_t b,uint64_t n,uint64_t npi) {
	asm("mulq %2 \n\t"
			"movq %%rax, %%r10 \n\t"
			"movq %%rdx, %%r11 \n\t"
			"mulq %3 \n\t"
			"mulq %4 \n\t"
			"addq %%r10, %%rax \n\t"
			"adcq %%r11, %%rdx \n\t"
			"xorq %%rax, %%rax \n\t"
			"subq %4, %%rdx \n\t"
			"cmovc %4, %%rax \n\t"
			"addq %%rdx, %%rax \n\t"
			: "=a"(a)
			: "0"(a), "r"(b), "r"(npi), "r"(n)
			: "rdx", "r10", "r11", "cc"
	);
  return a;
}
#else
// plain C versions to support non-x86 processors
static inline uint64_t absdiff(uint64_t a,uint64_t b) {
	return (a < b) ? b-a : a-b;
}

static inline uint64_t addmod(uint64_t a,uint64_t b,uint64_t n) {
	a += b;
	return (a < b || a >= n) ? a-n : a;
}

static inline uint64_t mulredc63(uint64_t x,uint64_t y,uint64_t n,uint64_t npi) {
	union { __uint128_t q; uint64_t l[2]; } u;

	u.q = (__uint128_t)x*y;
	u.q += (u.l[0]*npi)*(__uint128_t)n;
	u.l[1] -= n;
	return u.l[1]+(n&((int64_t)u.l[1]>>63));
}

static inline uint64_t mulredc64(uint64_t x,uint64_t y,uint64_t n,uint64_t npi) {
	union { __uint128_t q; uint64_t l[2]; } u;
	__uint128_t t;

	u.q = (__uint128_t)x*y;
	t = (u.l[0]*npi)*(__uint128_t)n;
	u.q += t;
	return (u.q < t || u.l[1] >= n) ? u.l[1]-n : u.l[1];
}
#endif

// compute -n^{-1} mod 2^64 by Newton's method
static inline uint64_t montgomery_inverse(uint64_t n) {
	uint32_t k = (3*n)^12; // correct mod 2^4
	k *= 2+k*(uint32_t)n; k *= 2+k*(uint32_t)n; k *= 2+k*(uint32_t)n;
	return (uint64_t)k*(2+(uint64_t)k*n);
}

// Montgomery representation of 1
#define montgomery_one(n) (1+(~(uint64_t)0)%(n))

// returns true if n passes a strong pseudoprime test to base 2
// assumes n is odd and > 1
static inline int isprobableprime(uint64_t n) {
	uint64_t one,minusone,x,nbar;
	int64_t q;
	int k;

	q = n>>1; k = __builtin_ctzl(q);
	q <<= __builtin_clzl(q);
	nbar = montgomery_inverse(n);

	if ((int64_t)n < 0) {
		// when n > 2^63, montgomery_one(n) = 2^64-n
		one = -n; minusone = n+n;
		x = addmod(one,one,n);
		for (q<<=1;q;q<<=1) {
			x = mulredc64(x,x,n,nbar);
			if (q < 0) x = addmod(x,x,n);
		}
		if (x == one || x == minusone) return 1;
		while (--k >= 0) {
			x = mulredc64(x,x,n,nbar);
			if (x == minusone) return 1;
		}
	} else {
		one = montgomery_one(n); minusone = n-one;
		x = addmod(one,one,n);
		for (q<<=1;q;q<<=1) {
			x = mulredc63(x,x,n,nbar);
			if (q < 0) x = addmod(x,x,n);
		}
		if (x == one || x == minusone) return 1;
		while (--k >= 0) {
			x = mulredc63(x,x,n,nbar);
			if (x == minusone) return 1;
		}
	}
	return 0;
}

// skips trial division
// n must be odd and > 1
static int oddisprime64(uint64_t n) {
	int i,j,t;

	if (!isprobableprime(n)) return 0;
	t = 62-__builtin_clzl(n); // n has t+2 bits
	i = pseudoprime_index[t]; j = pseudoprime_index[t+1];
	while (j-i > 1) {
		t = (i+j)>>1;
		if (pseudoprime_table[t] > n) j = t; else i = t;
	}
	return (pseudoprime_table[i] != n);
}

#define M ((uint64_t)1<<32)
int isprime64(uint64_t n) {
	int32_t r;

	if (n < M) {
		if (n <= 1) return 0;
		r = factor_index[(uint32_t)n%510510];
		if (r < 0) return ((int32_t)n == -r);
		return !factor_table[(uint32_t)n/510510*92160+r];
	}
	if (factor_index[n%510510] < 0) return 0;
	return oddisprime64(n);
}

// given that f is coprime to 510510, less than M and divides n,
// return a list of the prime factors of f and remove them from n
static int smallfactors(uint64_t *p,int *e,uint32_t f,uint64_t *n) {
	uint64_t x;
	uint32_t y;
	int k;

	k = 0;
	if (f > 1) {
		*n /= f;
		do {
			if (!(y=factor_table[f/510510*92160+factor_index[f%510510]])) {
				p[k] = y = f, e[k] = 1;
				f = 1;
			} else {
				p[k] = y, e[k] = 0;
				do f /= y, e[k]++; while (f % y == 0);
			}
			while (*n % y == 0)
				*n /= y, e[k]++;
			k++;
		} while (f > 1);
	}
	return k;
}

// assumes x is odd
static uint64_t oddgcd(uint64_t x,uint64_t y) {
	uint64_t t;

	if (y) {
		y >>= __builtin_ctzl(y);
		while (t=absdiff(x,y))
			if (x < y)
				y = t>>__builtin_ctzl(t);
			else
				x = t>>__builtin_ctzl(t);
	}
	return x;
}

int factor64(uint64_t *p,int *e,uint64_t n) {
// iterations can be any number > 14232 such that
//   iterations % maxstride == 0 && collision_index[iterations-1] != 0
// maxstride must be a power of 2
#define iterations 473344
#define maxstride 256
	uint64_t m,nbar,x,y,y0,f,one;
	uint32_t *ptr,s;
	int i,j,k,mask;

	k = 0;
	if (!(n & 1)) {
		if (!n) return 0;
		p[k] = 2, e[k] = __builtin_ctzl(n);
		n >>= e[k];
		k++;
	}

	f = oddgcd(n,(uint64_t)3*5*7*11*13*17*19*23*29*31*37*41*43*47*53);
	if (f > 1) {
		n /= f;
#define smallprime(s)\
		if (f % s == 0) {\
			f /= s; p[k] = s; e[k] = 1;\
			while (n % s == 0) n /= s, e[k]++;\
			k++; if (f == 1) goto done;\
		}

		// unroll the loop so the compiler can optimize the mod operations
		smallprime(3);  smallprime(5);  smallprime(7);
		smallprime(11); smallprime(13); smallprime(17);
		smallprime(19); smallprime(23); smallprime(29);
		smallprime(31); smallprime(37); smallprime(41);
		smallprime(43); smallprime(47); smallprime(53);
	}
done:

	if (n < M) {
		k += smallfactors(p+k,e+k,n,&n);
		return k;
	}
	if (oddisprime64(n)) {
		p[k] = n, e[k] = 1;
		k++;
		return k;
	}

	// Pollard rho
	nbar = montgomery_inverse(n);
	one = montgomery_one(n);
	m = n, y = f = one;
	for (i=1;i<iterations;i<<=1) {
		mask = (i < maxstride) ? i-1 : maxstride-1;
		x = y0 = y; j = 0;
		do {
			if ((int64_t)n < 0) {
				// one <= 2^64-n, so this does not overflow
				y = mulredc64(y,y+one,n,nbar);
				f = mulredc64(f,absdiff(y,x),n,nbar);
			} else {
				y = mulredc63(y,y+one,n,nbar);
				f = mulredc63(f,absdiff(y,x),n,nbar);
			}
			if (++j & mask) continue;
			if ((f=oddgcd(m,f))==1) {
				y0 = y;
				continue;
			}
			if (f >= M) {
				// backtrack to find exact cycle time
				y = y0, j -= (mask+1);
				if ((int64_t)n < 0)
					do {
						y = mulredc64(y,y+one,n,nbar);
						f = oddgcd(m,absdiff(y,x));
						j++;
					} while (f == 1);
				else
					do {
						y = mulredc63(y,y+one,n,nbar);
						f = oddgcd(m,absdiff(y,x));
						j++;
					} while (f == 1);

				ptr = collision_table+collision_index[i+j-2];
				while (f >= M)
					if (oddisprime64(f)) {
						p[k] = f, e[k] = 1;
						m /= f; k++; f = 1;
					} else {
						while (f % *ptr) ptr++;
						p[k] = *ptr++; e[k] = 0;
						do f /= p[k], m /= p[k], e[k]++; while (f % p[k] == 0);
						while (m % p[k] == 0) m /= p[k], e[k]++;
						k++;
					}
			}
			k += smallfactors(p+k,e+k,f,&m);
			if (m < M) {
				k += smallfactors(p+k,e+k,m,&m);
				return k;
			}
			if (oddisprime64(m)) {
				p[k] = m, e[k] = 1;
				k++;
				return k;
			}

			if (!(j&mask)) y0 = y;
			f = one;
		} while (j < i && i+j < iterations);
	}

	// the only remaining primes are those with long cycles
	// since iterations > 14232, every prime below 2^(64/3) has been
	//   considered, so m is either a prime square or semiprime
	ptr = collision_table+collision_index[iterations-1];
	while (m % *ptr) ptr++;
	p[k] = *ptr; m /= p[k];
	if (m == p[k])
		e[k] = 2;
	else {
		e[k++] = 1;
		p[k] = m, e[k] = 1;
	}
	k++;
	return k;
}
