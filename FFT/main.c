// 
// main.c
// 
//		Examples
// 

#include <conio.h>		// for _getch() function
#include <stdio.h>
#include <stdlib.h>		// for malloc(), calloc(), free() function
#include <Windows.h>

#include "fft.h"
#include "test.h"
#include "numeric.h"
#include "complexlib.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Main
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#define UNICODE

int main(void) {

	//Compare(31, true);
	//return 0;

	LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
	LARGE_INTEGER Frequency;

	QueryPerformanceFrequency(&Frequency);

	unsigned N = 32;
	bool isprint = true;

	dcomplex* x = 0;
	if ((x = (dcomplex*)calloc(N, sizeof(dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	printf("Raw data :\n");
	QueryPerformanceCounter(&StartingTime);

	for (unsigned i = 0; i < N; i++) {
		x[i].re = i;
		x[i].im = 0.;
		if (isprint) printf("%lf, %lf\n", x[i].re, x[i].im);
	}
	QueryPerformanceCounter(&EndingTime);
	ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
	ElapsedMicroseconds.QuadPart *= 1000000;
	ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
	printf("run_time: %lld us\n", ElapsedMicroseconds.QuadPart);
	printf("\n");

	{
		dcomplex* y = 0, * z = 0;

		printf("Transformed data :\n");
		QueryPerformanceCounter(&StartingTime);

		y = fft(x, N, backward);

		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
		if (isprint) {
			for (unsigned i = 0; i < N; i++) {
				printf("%lf, %lf\n", y[i].re, y[i].im);
			}
		}
		printf("run_time: %lld us\n", ElapsedMicroseconds.QuadPart);
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
				printf("%lf, %lf\n", z[i].re, z[i].im);
			}
		}
		printf("run_time: %lld us\n", ElapsedMicroseconds.QuadPart);
		printf("\n");

		if (N != 1) free(y);
		if (N != 1) free(z);
	}

	free(x);

	printf("\nPress any key to continue . . . ");
	if (_getch()) printf("\n");
	return 0;
}
