// 
// test.c
// 
//		Testing and comparing old functions and new functions
// 

#include "test.h"

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

static _Dcomplex* Butterfly_8point_fft(_Dcomplex* data, bool is_forward, bool keep_input) {

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
	if (!keep_input) free(data);

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

_Dcomplex* CooleyTukey_Radix2_Recursion_8pointAccelerate(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input) {

	if (N == 8) return Butterfly_8point_fft(data, is_forward, keep_input);
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

	Peven = CooleyTukey_Radix2_Recursion_8pointAccelerate(Peven, N / 2, is_forward, false);
	Podd = CooleyTukey_Radix2_Recursion_8pointAccelerate(Podd, N / 2, is_forward, false);

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

_Dcomplex* CooleyTukey_Radix2_nonRecursion(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input) {
	int b = (is_forward) ? 1 : -1;
	_Dcomplex* result = bitReverseSorting(data, N);
	if (!keep_input) free(data);
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

_Dcomplex* Radar(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input, _Dcomplex* (*fft)()) {

	if (N == 1) return data;

	unsigned M = nextPowerOf2(2 * (N - 1) - 1);
	unsigned g = findPrimitiveRoot(N);

	// a[q] = x_g^q = x[ pow(g, q) (mod N) ]
	_Dcomplex* a = 0;
	if ((a = (_Dcomplex*)calloc(M, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	// b[q] = x_g^-q = x[ pow(g, -q) (mod N) ]
	_Dcomplex* b = 0;
	if ((b = (_Dcomplex*)calloc(M, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	unsigned* g_mod_minus_q = 0;
	if ((g_mod_minus_q = (unsigned*)calloc((size_t)N - 1, sizeof(unsigned))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	int k = (is_forward) ? 1 : -1;
	*g_mod_minus_q = expModuloNInverse(g, N);
	*a = *(data + expModuloN(g, 0, N));
	*b = cexp(_Cbuild(0., -2 * PI * k / N));
	for (unsigned q = 1; q < N - 1; q++) {
		*(a + q + M - N + 1) = *(data + expModuloN(g, q, N));
		*(b + q) = cexp(_Cbuild(0., -2 * PI * *(g_mod_minus_q + q - 1) * k / N));
		// pow(g, -2) (mod N)
		//      = { [pow(g, -1) (mod N)] * [pow(g, -1) (mod N)] } % N
		//		= [ pow(g, -1) * pow(g, -1) ] (mod N)
		*(g_mod_minus_q + q) = (*(g_mod_minus_q + q - 1) * *g_mod_minus_q) % N;
	}
	for (unsigned i = N - 1; i < M; i++) {
		*(b + i) = *(b + i - N + 1);
	}

	_Dcomplex* fft_a = fft(a, M, true, true);
	free(a);
	_Dcomplex* fft_b = fft(b, M, true, true);
	free(b);

	_Dcomplex* fft_ab = 0;
	if ((fft_ab = (_Dcomplex*)calloc(M, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned i = 0; i < M; i++) {
		*(fft_ab + i) = cmul(*(fft_a + i), *(fft_b + i));
	}
	free(fft_a);
	free(fft_b);

	_Dcomplex* conv_ab = fft(fft_ab, M, false, true);
	free(fft_ab);

	_Dcomplex* result = 0;
	if ((result = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned n = 0; n < N; n++) {
		*result = cadd(*result, *(data + n));
	}
	*(result + 1) = cadd(*data, _Cmulcr(*conv_ab, 1. / M));
	for (unsigned q = 1; q < N - 1; q++) {
		*(result + *(g_mod_minus_q + q - 1)) = cadd(*data, _Cmulcr(*(conv_ab + q), 1. / M));
	}
	free(g_mod_minus_q);
	free(conv_ab);
	return result;
}

_Dcomplex* CooleyTukey_HybridRadix(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input, _Dcomplex* (*fft)()) {

	if (N == 1) return data;
	if (N == 2) return fft(data, N, is_forward, false);
	if (isprime(N)) return Radar(data, N, is_forward, false, fft);
	if (isPowerOf2(N)) return fft(data, N, is_forward, false);

	unsigned* factors = factor(N);
	unsigned N1 = *(factors + 1), N2 = N / N1;
	free(factors);

	_Dcomplex** X = 0;
	if ((X = (_Dcomplex**)calloc(N1, sizeof(_Dcomplex*))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned i = 0; i < N1; i++) {
		if ((*(X + i) = (_Dcomplex*)calloc(N2, sizeof(_Dcomplex))) == NULL) {
			printf("Allocation failed!\n");
			return 0;
		}
	}

	for (unsigned n1 = 0; n1 < N1; n1++) {
		for (unsigned n2 = 0; n2 < N2; n2++) {
			X[n1][n2] = data[N1 * n2 + n1];
		}
	}
	if (!keep_input) free(data);

	int k = (is_forward) ? 1 : -1;
	for (unsigned n1 = 0; n1 < N1; n1++) {
		_Dcomplex* inner = 0;
		if ((inner = (_Dcomplex*)calloc(N2, sizeof(_Dcomplex))) == NULL) {
			printf("Allocation failed!\n");
			return 0;
		}
		for (unsigned i = 0; i < N2; i++) {
			*(inner + i) = X[n1][i];
		}
		_Dcomplex* fft_inner = CooleyTukey_HybridRadix(inner, N2, is_forward, false, fft);
		for (unsigned n2 = 0; n2 < N2; n2++) {
			X[n1][n2] = cmul(fft_inner[n2], cexp(_Cbuild(0., -2 * PI * n1 * n2 * k / N)));
		}
		free(fft_inner);
	}

	for (unsigned n2 = 0; n2 < N2; n2++) {
		_Dcomplex* outer = 0;
		if ((outer = (_Dcomplex*)calloc(N1, sizeof(_Dcomplex))) == NULL) {
			printf("Allocation failed!\n");
			return 0;
		}
		for (unsigned i = 0; i < N1; i++) {
			*(outer + i) = X[i][n2];
		}
		_Dcomplex* fft_outer = CooleyTukey_HybridRadix(outer, N1, is_forward, false, fft);
		for (unsigned n1 = 0; n1 < N1; n1++) {
			X[n1][n2] = fft_outer[n1];
		}
		free(fft_outer);
	}

	_Dcomplex* result = 0;
	if ((result = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned n1 = 0; n1 < N1; n1++) {
		for (unsigned n2 = 0; n2 < N2; n2++) {
			*(result + N2 * n1 + n2) = X[n1][n2];
		}
	}

	for (unsigned i = 0; i < N1; i++) {
		free(*(X + i));
	}
	free(X);

	return result;
}

_Dcomplex* Bluestein(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input, _Dcomplex* (*fft)()) {

	unsigned M = nextPowerOf2(2 * N - 1);

	_Dcomplex* a = 0;
	if ((a = (_Dcomplex*)calloc(M, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	_Dcomplex* b = 0;
	if ((b = (_Dcomplex*)calloc(M, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	int k = (is_forward) ? 1 : -1;
	double common = k * PI / N;
	for (unsigned n = 0; n < N; n++) {
		double twiddle = common * n * n;
		*(a + n) = cmul(*(data + n), cexp(_Cbuild(0, -twiddle)));
		*(b + n) = cexp(_Cbuild(0, twiddle));
	}
	for (unsigned n = M - N + 1; n < M; n++) {
		unsigned p = M - n;
		*(b + n) = cexp(_Cbuild(0, common * p * p));
	}

	_Dcomplex* fft_a = fft(a, M, true, true);
	free(a);
	_Dcomplex* fft_b = fft(b, M, true, true);
	free(b);

	_Dcomplex* fft_ab = 0;
	if ((fft_ab = (_Dcomplex*)calloc(M, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned i = 0; i < M; i++) {
		*(fft_ab + i) = cmul(*(fft_a + i), *(fft_b + i));
	}
	free(fft_a);
	free(fft_b);

	_Dcomplex* conv_ab = fft(fft_ab, M, false, true);
	free(fft_ab);

	_Dcomplex* result = 0;
	if ((result = (_Dcomplex*)realloc(conv_ab, N * sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		free(conv_ab);
		return 0;
	}
	conv_ab = NULL;

	for (unsigned n = 0; n < N; n++) {
		*(result + n) = _Cmulcr(cmul(*(result + n), cexp(_Cbuild(0, -common * n * n))), 1. / M);
	}

	return result;
}

void test(_Dcomplex* x, unsigned N, bool is_display, _Dcomplex* (*basic_fft)(), _Dcomplex* (*target_fft)()) {

	LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
	LARGE_INTEGER Frequency;

	QueryPerformanceFrequency(&Frequency);

	{
		_Dcomplex* y = 0;

		QueryPerformanceCounter(&StartingTime);

		y = target_fft(x, N, true, true, basic_fft);

		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;

		if (is_display) {
			for (unsigned i = 0; i < N; i++) {
				cprint(*(y + i));
			}
		}
		printf("run_time: %lld \u03bcs\n", ElapsedMicroseconds.QuadPart);
		printf("\n");

		if (N != 1) free(y);
	}

	return;
}

void Compare(unsigned N, bool is_display) {

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
		if (is_display) cprint(x[i]);
	}
	QueryPerformanceCounter(&EndingTime);
	ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
	ElapsedMicroseconds.QuadPart *= 1000000;
	ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
	printf("run_time: %lld \u03bcs\n", ElapsedMicroseconds.QuadPart);
	printf("\n");

	if (isPowerOf2(N)) {
		printf("N is power of 2\n\n");
		
		printf("Cooley-Tukey radix-2 algorithm with recursion :\n");
		test(x, N, is_display, NULL, CooleyTukey_Radix2_Recursion);

		printf("Cooley-Tukey radix-2 algorithm with Butterfly Diagram 8-point accelerate recursion :\n");
		test(x, N, is_display, NULL, CooleyTukey_Radix2_Recursion_8pointAccelerate);

		printf("Cooley-Tukey radix-2 algorithm without recursion :\n");
		test(x, N, is_display, NULL, CooleyTukey_Radix2_nonRecursion);

	}
	else if (isprime(N)) {
		printf("N is a prime number\n\n");

		printf("Radar's algorithm with recursion Cooley-Tukey radix-2 algorithm :\n");
		test(x, N, is_display, CooleyTukey_Radix2_Recursion, Radar);

		printf("Radar's algorithm with 8-point butterfly diagram accelerated recursion Cooley-Tukey radix-2 algorithm :\n");
		test(x, N, is_display, CooleyTukey_Radix2_Recursion_8pointAccelerate, Radar);

		printf("Radar's algorithm with non-recursion Cooley-Tukey radix-2 algorithm :\n");
		test(x, N, is_display, CooleyTukey_Radix2_nonRecursion, Radar);
	}
	else {
		printf("N is a composite number\n\n");

		printf("Cooley-Tukey Hybrid-Radix algorithm with Radar's Algorithm having recursion Cooley-Tukey radix-2 algorithm :\n");
		test(x, N, is_display, CooleyTukey_Radix2_Recursion, CooleyTukey_HybridRadix);

		printf("Cooley-Tukey Hybrid-Radix algorithm with Radar's Algorithm having 8-point butterfly diagram accelerated Cooley-Tukey radix-2 algorithm :\n");
		test(x, N, is_display, CooleyTukey_Radix2_Recursion_8pointAccelerate, CooleyTukey_HybridRadix);

		printf("Cooley-Tukey Hybrid-Radix algorithm with Radar's Algorithm having non-recursion Cooley-Tukey radix-2 algorithm :\n");
		test(x, N, is_display, CooleyTukey_Radix2_nonRecursion, CooleyTukey_HybridRadix);
	}

	printf("Bluestein's algorithm with recursion Cooley-Tukey radix-2 algorithm :\n");
	test(x, N, is_display, CooleyTukey_Radix2_Recursion, Bluestein);

	printf("Bluestein's algorithm with 8-point butterfly diagram accelerated recursion Cooley-Tukey radix-2 algorithm :\n");
	test(x, N, is_display, CooleyTukey_Radix2_Recursion_8pointAccelerate, Bluestein);

	printf("Bluestein's algorithm with non-recursion Cooley-Tukey radix-2 algorithm :\n");
	test(x, N, is_display, CooleyTukey_Radix2_nonRecursion, Bluestein);

	free(x);
	return;
}
