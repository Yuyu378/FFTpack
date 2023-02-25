// 
// main.c
// 
//		Examples
// 

#include <conio.h>		// for _getch() function
#include <stdio.h>
#include <stdlib.h>		// for malloc(), calloc(), free() function
#include <complex.h>
#include <Windows.h>

#include "fft.h"
#include "test.h"
#include "numeric.h"
#include "excomplex.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Main
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

int main(void) {

	//Compare(31, true);
	//return 0;

	LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
	LARGE_INTEGER Frequency;

	QueryPerformanceFrequency(&Frequency);

	unsigned N = 32;
	bool isprint = true;

	_Dcomplex* x = 0;
	if ((x = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	printf("Raw data :\n");
	QueryPerformanceCounter(&StartingTime);

	for (unsigned i = 0; i < N; i++) {
		x[i]._Val[0] = i;
		x[i]._Val[1] = 0.;
		if (isprint) cprint(x[i]);
	}
	QueryPerformanceCounter(&EndingTime);
	ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
	ElapsedMicroseconds.QuadPart *= 1000000;
	ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
	printf("run_time: %lld \u03bcs\n", ElapsedMicroseconds.QuadPart);
	printf("\n");

	{
		_Dcomplex* y = 0, * z = 0;

		printf("Transformed data :\n");
		QueryPerformanceCounter(&StartingTime);

		y = fft(x, N, backward);

		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
		if (isprint) {
			for (unsigned i = 0; i < N; i++) {
				cprint(*(y + i));
			}
		}
		printf("run_time: %lld \u03bcs\n", ElapsedMicroseconds.QuadPart);
		printf("\n");

		printf("Inverse Transformed data :\n");
		QueryPerformanceCounter(&StartingTime);

		z = ifft(y, N, backward);

		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
		if (isprint) {
			for (unsigned i = 0; i < N; i++) {
				cprint(*(z + i));
			}
		}
		printf("run_time: %lld \u03bcs\n", ElapsedMicroseconds.QuadPart);
		printf("\n");

		if (N != 1) free(y);
		if (N != 1) free(z);
	}

	free(x);

	printf("\nPress any key to continue . . . ");
	if (_getch()) printf("\n");
	return 0;
}
