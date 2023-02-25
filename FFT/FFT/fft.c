// 
// fft.c
// 
//		Fourier Transform algorithms implement
// 

#include "fft.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Functions
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

static double _get_forward_norm(unsigned N, norm_mode mode) {
	switch (mode) {
	case backward:
		return 1.;
	case ortho:
		return sqrt(N);
	case forward:
		return 1. * N;
	default:
		return 1.;
	}
}

static double _get_backward_norm(unsigned N, norm_mode mode) {
	switch (mode) {
	case backward:
		return 1. * N;
	case ortho:
		return sqrt(N);
	case forward:
		return 1.;
	default:
		return 1. * N;
	}
}

static _Dcomplex* _execute_dft(_Dcomplex* data, unsigned int N, bool is_forward) {

	_Dcomplex* r = 0;
	if ((r = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	int b = (is_forward) ? 1 : -1;
	for (unsigned k = 0; k < N; k++) {
		for (unsigned n = 0; n < N; n++) {
			*(r + k) = cadd(*(r + k),
				cmul(*(data + n), cexp(_Cbuild(0., -2. * PI * k * n * b / N)))
			);
		}
	}

	return r;
}

static _Dcomplex* _raw_dft(_Dcomplex* data, unsigned int N, bool is_forward, double fct) {
	_Dcomplex* result = _execute_dft(data, N, is_forward);
	for (unsigned i = 0; i < N; i++) {
		*(result + i) = cmul(*(result + i), _Cbuild(fct, 0.));
	}
	return result;
}

_Dcomplex* dft(_Dcomplex* data, unsigned int N, norm_mode mode) {
	double fct = 1. / _get_forward_norm(N, mode);
	return _raw_dft(data, N, true, fct);
}

_Dcomplex* idft(_Dcomplex* data, unsigned int N, norm_mode mode) {
	double fct = 1. / _get_backward_norm(N, mode);
	return _raw_dft(data, N, false, fct);
}

static _Dcomplex* _execute_fft(_Dcomplex* data, unsigned int N, bool is_forward) {
	if (N == 1) return data;
	if (N == 2) return _raw_radix2fft(data, N, is_forward, true);
	if (isprime(N)) return _raw_radarfft(data, N, is_forward, true);
	if (isPowerOf2(N)) return _raw_radix2fft(data, N, is_forward, true);
	else return _raw_bluesteinfft(data, N, is_forward, true);
}

static _Dcomplex* _raw_fft(_Dcomplex* data, unsigned int N, bool is_forward, double fct) {
	_Dcomplex* result = _execute_fft(data, N, is_forward);
	for (unsigned i = 0; i < N; i++) {
		*(result + i) = _Cmulcr(*(result + i), fct);
	}
	return result;
}

_Dcomplex* fft(_Dcomplex* data, unsigned int N, norm_mode mode) {
	double fct = 1. / _get_forward_norm(N, mode);
	return _raw_fft(data, N, true, fct);
}

_Dcomplex* ifft(_Dcomplex* data, unsigned int N, norm_mode mode) {
	double fct = 1. / _get_backward_norm(N, mode);
	return _raw_fft(data, N, false, fct);
}

static _Dcomplex* _raw_radix2fft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input) {
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

static _Dcomplex* _raw_radarfft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input) {

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

	_Dcomplex* fft_a = _raw_radix2fft(a, M, true, false);
	_Dcomplex* fft_b = _raw_radix2fft(b, M, true, false);

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

	_Dcomplex* conv_ab = _raw_radix2fft(fft_ab, M, false, false);

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
	if (!keep_input) free(data);
	free(g_mod_minus_q);
	free(conv_ab);
	return result;
}

static _Dcomplex* _raw_bluesteinfft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input) {

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
	if (!keep_input) free(data);
	for (unsigned n = M - N + 1; n < M; n++) {
		unsigned p = M - n;
		*(b + n) = cexp(_Cbuild(0, common * p * p));
	}

	_Dcomplex* fft_a = _raw_radix2fft(a, M, true, false);
	_Dcomplex* fft_b = _raw_radix2fft(b, M, true, false);

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

	_Dcomplex* conv_ab = _raw_radix2fft(fft_ab, M, false, false);

	_Dcomplex* result = 0;
	if ((result = (_Dcomplex*)realloc(conv_ab, N * sizeof(_Dcomplex))) == NULL) {
		printf("Reallocation failed!\n");
		free(conv_ab);
		return 0;
	}
	conv_ab = NULL;

	for (unsigned n = 0; n < N; n++) {
		*(result + n) = _Cmulcr(cmul(*(result + n), cexp(_Cbuild(0, -common * n * n))), 1. / M);
	}

	return result;
}
