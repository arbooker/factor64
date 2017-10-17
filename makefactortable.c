#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define M ((uint64_t)1<<32)
int main(void) {
	uint64_t n;
	uint32_t k;
	uint16_t *factor;
	int32_t r;

	factor = calloc(M>>1,sizeof(*factor));
	for (k=(uint32_t)sqrtl((long double)M);k>=3;k--)
		if (k%2 && k%3 && k%5 && k%7 && k%11 && k%13 && k%17)
			for (n=k*3;n<M;n+=k*2)
				factor[n>>1] = (uint16_t)k;

	for (n=0,k=0;n<510510;n++) {
		if (n % 2 == 0) r = -2;
		else if (n % 3 == 0) r = -3;
		else if (n % 5 == 0) r = -5;
		else if (n % 7 == 0) r = -7;
		else if (n % 11 == 0) r = -11;
		else if (n % 13 == 0) r = -13;
		else if (n % 17 == 0) r = -17;
		else r = k++;
		fwrite(&r,sizeof(r),1,stdout);
	}
			
	for (n=1,k=0;n<M;n+=2)
		if (n % 3 && n % 5 && n % 7 && n % 11 && n % 13 && n % 17) {
			fwrite(&factor[n>>1],sizeof(factor[0]),1,stdout);
			k++;
		}

	return 0;
}
