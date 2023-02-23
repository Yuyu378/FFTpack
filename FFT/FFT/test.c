// 
// test.c
// 
//		Testing and comparing old functions and new functions
// 

#include "test.h"

static _Dcomplex* butterfly_8point_fft(_Dcomplex* data, bool is_forward) {

	int k = (is_forward) ? 1 : -1;

	_Dcomplex wn[4] = {
		_Cbuild(1., 0.),										// exp(-i0PI/8)
		_Cbuild(cos(2. * PI / 8.), -sin(2. * PI / 8.) * k),		// exp(-i2PI/8)
		_Cbuild(cos(4. * PI / 8.), -sin(4. * PI / 8.) * k),		// exp(-i4PI/8)
		_Cbuild(cos(6. * PI / 8.), -sin(6. * PI / 8.) * k)		// exp(-i6PI/8)
	};

	_Dcomplex tmp1[8] = { 0 };
	_Dcomplex tmp2[8] = { 0 };

	// Decimation In Time - Fast Fourier Transform

	tmp1[0] = cadd(data[0], data[4]);
	tmp1[1] = csub(data[0], data[4]);
	tmp1[2] = cadd(data[2], data[6]);
	tmp1[3] = cmul(csub(data[2], data[6]), wn[2]);
	tmp1[4] = cadd(data[1], data[5]);
	tmp1[5] = csub(data[1], data[5]);
	tmp1[6] = cadd(data[3], data[7]);
	tmp1[7] = cmul(csub(data[3], data[7]), wn[2]);

	tmp2[0] = cadd(tmp1[0], tmp1[2]);
	tmp2[1] = cadd(tmp1[1], tmp1[3]);
	tmp2[2] = csub(tmp1[0], tmp1[2]);
	tmp2[3] = csub(tmp1[1], tmp1[3]);
	tmp2[4] = cadd(tmp1[4], tmp1[6]);
	tmp2[5] = cmul(cadd(tmp1[5], tmp1[7]), wn[1]);
	tmp2[6] = cmul(csub(tmp1[4], tmp1[6]), wn[2]);
	tmp2[7] = cmul(csub(tmp1[5], tmp1[7]), wn[3]);

	_Dcomplex* result = 0;
	if ((result = (_Dcomplex*)calloc(8, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	*(result + 0) = cadd(tmp2[0], tmp2[4]);
	*(result + 1) = cadd(tmp2[1], tmp2[5]);
	*(result + 2) = cadd(tmp2[2], tmp2[6]);
	*(result + 3) = cadd(tmp2[3], tmp2[7]);
	*(result + 4) = csub(tmp2[0], tmp2[4]);
	*(result + 5) = csub(tmp2[1], tmp2[5]);
	*(result + 6) = csub(tmp2[2], tmp2[6]);
	*(result + 7) = csub(tmp2[3], tmp2[7]);

	return result;
}

_Dcomplex* CooleyTukey_Radix2_Recursion(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input) {
	
	if (N == 1) return data;

	_Dcomplex* Peven = 0;
	if ((Peven = (_Dcomplex*)calloc((N / 2), sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	_Dcomplex* Podd = 0;
	if ((Podd = (_Dcomplex*)calloc((N / 2), sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	for (unsigned i = 0; i < N / 2; i++) {
		*(Peven + i) = *(data + i * 2);
		*(Podd + i) = *(data + i * 2 + 1);
	}
	if (!keep_input) free(data);

	Peven = CooleyTukey_Radix2_Recursion(Peven, N / 2, is_forward, false);
	Podd = CooleyTukey_Radix2_Recursion(Podd, N / 2, is_forward, false);

	_Dcomplex* result = 0;
	if ((result = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	_Dcomplex w = { 0 };
	int k = (is_forward) ? 1 : -1;
	for (unsigned i = 0; i < N / 2; i++) {
		w = cexp(_Cbuild(0., -2. * PI * i * k / N));
		_Dcomplex wPodd = cmul(w, *(Podd + i));
		*(result + i) = cadd(*(Peven + i), wPodd);
		*(result + i + N / 2) = csub(*(Peven + i), wPodd);
	}

	free(Peven);
	free(Podd);

	return result;
}

_Dcomplex* CooleyTukey_Radix2_Recursion_8pointAccelerate(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input) {

	if (N == 8) return butterfly_8point_fft(data, is_forward);
	if (N == 1) return data;

	_Dcomplex* Peven = 0;
	if ((Peven = (_Dcomplex*)calloc((N / 2), sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	_Dcomplex* Podd = 0;
	if ((Podd = (_Dcomplex*)calloc((N / 2), sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	for (unsigned i = 0; i < N / 2; i++) {
		*(Peven + i) = *(data + i * 2);
		*(Podd + i) = *(data + i * 2 + 1);
	}
	if (!keep_input) free(data);

	Peven = CooleyTukey_Radix2_Recursion(Peven, N / 2, is_forward, false);
	Podd = CooleyTukey_Radix2_Recursion(Podd, N / 2, is_forward, false);

	_Dcomplex* result = 0;
	if ((result = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	_Dcomplex w = { 0 };
	int k = (is_forward) ? 1 : -1;
	for (unsigned i = 0; i < N / 2; i++) {
		w = cexp(_Cbuild(0., -2. * PI * i * k / N));
		_Dcomplex wPodd = cmul(w, *(Podd + i));
		*(result + i) = cadd(*(Peven + i), wPodd);
		*(result + i + N / 2) = csub(*(Peven + i), wPodd);
	}

	free(Peven);
	free(Podd);

	return result;
}

_Dcomplex* CooleyTukey_Radix2_nonRecursion(_Dcomplex* data, unsigned N, bool is_forward) {
	int b = (is_forward) ? 1 : -1;
	_Dcomplex* result = bitReverseSorting(data, N);
	for (unsigned i = 2; i <= N; i <<= 1) {
		_Dcomplex step = cexp(_Cbuild(0., -2. * PI * b / i));
		for (unsigned j = 0; j < N; j += i) {
			_Dcomplex w = { 1, 0 };
			for (unsigned k = j; k < j + i / 2; k++) {
				_Dcomplex Peven = result[k];
				_Dcomplex wPodd = cmul(w, result[k + i / 2]);
				result[k] = cadd(Peven, wPodd);
				result[k + i / 2] = csub(Peven, wPodd);
				w = cmul(w, step);
			}
		}
	}
	return result;
}

void Compare(unsigned N) {

	if (!isPowerOf2(N)) {
		printf("N must be the power of 2.\n");
		return;
	}

	LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
	LARGE_INTEGER Frequency;

	QueryPerformanceFrequency(&Frequency);

	_Dcomplex* x = 0;
	if ((x = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return;
	}

	printf("Raw data :\n");
	QueryPerformanceCounter(&StartingTime);

	for (unsigned i = 0; i < N; i++) {
		x[i]._Val[0] = i;
		x[i]._Val[1] = 0.;
		cprint(x[i]);
	}
	QueryPerformanceCounter(&EndingTime);
	ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
	ElapsedMicroseconds.QuadPart *= 1000000;
	ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
	printf("run_time: %lld \u03bcs\n", ElapsedMicroseconds.QuadPart);
	printf("\n");

	{

		_Dcomplex* y = 0;

		QueryPerformanceCounter(&StartingTime);

		y = CooleyTukey_Radix2_Recursion_8pointAccelerate(x, N, true, true);

		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;

		printf("Cooley-Tukey radix-2 algorithm with recursion with Butterfly Diagram 8-point FFT :\n");
		for (unsigned i = 0; i < N; i++) {
			cprint(*(y + i));
		}
		printf("run_time: %lld \u03bcs\n", ElapsedMicroseconds.QuadPart);
		printf("\n");

		if (N != 1) free(y);
	}

	{
		_Dcomplex* y = 0;

		QueryPerformanceCounter(&StartingTime);

		y = CooleyTukey_Radix2_Recursion(x, N, true, true);

		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;

		printf("Cooley-Tukey radix-2 algorithm with recursion :\n");
		for (unsigned i = 0; i < N; i++) {
			cprint(*(y + i));
		}
		printf("run_time: %lld \u03bcs\n", ElapsedMicroseconds.QuadPart);
		printf("\n");

		if (N != 1) free(y);
	}

	{
		_Dcomplex* y = 0;

		QueryPerformanceCounter(&StartingTime);

		y = CooleyTukey_Radix2_nonRecursion(x, N, true);

		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;

		printf("Cooley-Tukey radix-2 algorithm without recursion :\n");
		for (unsigned i = 0; i < N; i++) {
			cprint(*(y + i));
		}
		printf("run_time: %lld \u03bcs\n", ElapsedMicroseconds.QuadPart);
		printf("\n");

		if (N != 1) free(y);
	}
	
	free(x);
	return;
}
