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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int initfactor64(const char *);
int factor64(uint64_t *,int *,uint64_t);
int isprime64(uint64_t);

static void sort_factors(int k,uint64_t *p,int *e) {
	int i,j,te;
	uint64_t tp;

	for (i=1;i<k;i++)
		for (j=i;j>0 && p[j-1]>p[j];j--) {
			tp=p[j-1]; p[j-1]=p[j]; p[j]=tp;
			te=e[j-1]; e[j-1]=e[j]; e[j]=te;
		}
}

int main(int argc,char *argv[]) {
	uint64_t n,p[15];
	int i,j,k,e[15];
	char line[1024];

	if (initfactor64("factor.bin") < 0) {
		fprintf(stderr, "Cannot read factor data\n");
		return -1;
	}

	for (i=1;;) {
		if (argc <= 1) { // no args, read from stdin
			if (!fgets(line,sizeof(line),stdin)) break;
			n = (uint64_t)strtoull(line,NULL,10);
		} else {
			if (i >= argc) break;
			n = (uint64_t)strtoull(argv[i++],NULL,10);
		}
		k = factor64(p,e,n);
		sort_factors(k,p,e);
		printf("%llu:",(unsigned long long)n);
		for (j=0;j<k;j++)
			while (--e[j] >= 0)
				printf(" %llu",(unsigned long long)p[j]);
		printf("\n");
	}

	return 0;
}
