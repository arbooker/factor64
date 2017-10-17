This is the documentation for factor64, a library for primality testing
and factoring of integers less than 2^64, by Andrew R. Booker. It is free
software, distributed under the terms of the GPLv3. For more information,
see the file LICENSE.txt.

CREDITS: Thanks to Ben Buhrow and Dana Jacobsen for contributing code
and ideas. The primality test relies on the list of 64-bit base-2 strong
pseudoprimes computed by Jan Feitsma et al.; see
http://www.janfeitsma.nl/math/psp2/index

The library is supported by a data file (called factor.bin by default)
that is downloaded in compressed form (about 846MB) and reconstituted at
compile time; see the Makefile for details. A sample program for calling
the library, emulating the basic usage of the Unix factor utility,
is included in the file sample.c.

The source code for the library consists of a single C file, factor64.c.
It exports three functions:

int initfactor64(const char *filename);
Initializes the data table pointed to by filename. The table uses about
2.5GB of memory, but the file is accessed via memory mapping, so separate
processes/threads will share a single copy. The return value is zero
on success and a negative number if an error occurs.

int isprime64(uint64_t n);
Returns 1 if n is prime, and 0 otherwise.

int factor64(uint64_t p[],int e[],uint64_t n);
Places the prime factors of n in p and their exponents in e, and returns
the number of factors. The maximum return value is 15, and the output
arrays should be sized accordingly. The prime factors are not necessarily
returned in order; see sample.c for example code to sort the output
if needed.
