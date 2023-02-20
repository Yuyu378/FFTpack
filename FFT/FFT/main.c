// 
// main.c
// 
//		Examples & Tests
// 

#include <conio.h>		// for _getch() function
#include <stdio.h>
#include <stdlib.h>		// for malloc(), calloc(), free() function
#include <complex.h>

#include "fft.h"
#include "numeric.h"
#include "excomplex.h"


//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Main
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

int main(void) {

	unsigned N = 24;

	_Dcomplex* x = 0;
	if ((x = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	for (unsigned i = 0; i < N; i++) {
		x[i]._Val[0] = i;
		x[i]._Val[1] = 0.;
		cprint(x[i]);
	}
	printf("\n");

	_Dcomplex* y = 0, * z = 0;

	y = fft(x, N, backward);
	
	for (unsigned i = 0; i < N; i++) {
		cprint(*(y + i));
	}
	printf("\n");

	z = ifft(y, N, backward);
	for (unsigned i = 0; i < N; i++) {
		cprint(*(z + i));
	}
	printf("\n");

	free(x);
	if (N != 1) free(y);
	if (N != 1) free(z);

	printf("\nPress any key to continue . . . ");
	if (_getch()) printf("\n");
	return 0;
}
