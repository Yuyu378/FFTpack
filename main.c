//	
// Jan 28, 2023
// Feb 06, 2023
//
// main.c
//

#include <conio.h>		// for _getch() function
#include <stdio.h>
#include <stdlib.h>		// for malloc(), calloc(), free() function

#include "fft.h"		// for dft(), idft(), fft(), ifft() function
#include "mymath.h"		// for cprint(), complex data computing, fft computing


// Main ----------------------------------------------------------------------------------

int main(void) {

	// An example of use
	
	unsigned N = 32;

	_Dcomplex* x = 0;
	if ((x = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	for (unsigned i = 0; i < N; i++) {
		x[i]._Val[0] = i;
		x[i]._Val[1] = 0.;
		cprint(x[i]);			// print complex data
	}
	printf("\n");

	_Dcomplex* y = 0, * z = 0;

	// N must be the same length as the x array
	y = fft(x, N, backward);
	
	for (unsigned i = 0; i < N; i++) {
		cprint(*(y + i));
	}
	printf("\n");

	// The use of y = ifft(y, N, backward) should be avoided, 
	// which may lead to memory leaks.
	z = ifft(y, N, backward);
	for (unsigned i = 0; i < N; i++) {
		cprint(*(z + i));
	}
	printf("\n");

	// If the array length is 1, free() is redundant.
	if (N != 1) free(x);
	if (N != 1) free(y);
	if (N != 1) free(z);

	printf("\nPress any key to continue . . . ");
	if (_getch()) printf("\n");
	return 0;
}
