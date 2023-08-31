// 
// fft.c
// 
//		Fourier Transform algorithms implement
// 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

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

static dcomplex* _execute_dft(dcomplex* data, unsigned int N, bool is_forward) {

	dcomplex* r = 0;
	if ((r = (dcomplex*)calloc(N, sizeof(dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	int b = (is_forward) ? 1 : -1;
	for (unsigned k = 0; k < N; k++) {
		for (unsigned n = 0; n < N; n++) {
			*(r + k) = cadd(*(r + k),
				cmul(*(data + n), cexp(cbuild(0., -2. * PI * k * n * b / N)))
			);
		}
	}

	return r;
}

static dcomplex* _raw_dft(dcomplex* data, unsigned int N, bool is_forward, double fct) {
	dcomplex* result = _execute_dft(data, N, is_forward);
	for (unsigned i = 0; i < N; i++) {
		*(result + i) = cmul(*(result + i), cbuild(fct, 0.));
	}
	return result;
}

dcomplex* dft(dcomplex* data, unsigned int N, norm_mode mode) {
	double fct = 1. / _get_forward_norm(N, mode);
	return _raw_dft(data, N, true, fct);
}

dcomplex* idft(dcomplex* data, unsigned int N, norm_mode mode) {
	double fct = 1. / _get_backward_norm(N, mode);
	return _raw_dft(data, N, false, fct);
}

static dcomplex* _execute_fft(dcomplex* data, unsigned int N, bool is_forward) {
	if (N == 1) return data;
	if (N == 2) return _raw_radix2fft(data, N, is_forward, true);
	if (isprime(N)) return _raw_radarfft(data, N, is_forward, true);
	if (isPowerOf2(N)) return _raw_radix2fft(data, N, is_forward, true);
	else return _raw_bluesteinfft(data, N, is_forward, true);
}

static dcomplex* _raw_fft(dcomplex* data, unsigned int N, bool is_forward, double fct) {
	dcomplex* result = _execute_fft(data, N, is_forward);
	for (unsigned i = 0; i < N; i++) {
		*(result + i) = _Cmulcr(*(result + i), fct);
	}
	return result;
}

dcomplex* fft(dcomplex* data, unsigned int N, norm_mode mode) {
	double fct = 1. / _get_forward_norm(N, mode);
	return _raw_fft(data, N, true, fct);
}

dcomplex* ifft(dcomplex* data, unsigned int N, norm_mode mode) {
	double fct = 1. / _get_backward_norm(N, mode);
	dcomplex* ans = _raw_fft(data, N, false, fct);
	for (unsigned int i = 0; i < N; i++) {
		if (ans[i].im < __DBL_EPSILON__) ans[i].im = 0.0;
	}
	return ans;
}

static dcomplex* _raw_radix2fft(dcomplex* data, unsigned N, bool is_forward, bool keep_input) {
	int b = (is_forward) ? 1 : -1;
	dcomplex* result = bitReverseSorting(data, N);
	if (!keep_input) free(data);
	for (unsigned i = 2; i <= N; i <<= 1) {
		dcomplex step = cexp(cbuild(0., -2. * PI * b / i));
		for (unsigned j = 0; j < N; j += i) {
			dcomplex w = { 1, 0 };
			for (unsigned k = j; k < j + i / 2; k++) {
				dcomplex Peven = result[k];
				dcomplex wPodd = cmul(w, result[k + i / 2]);
				result[k] = cadd(Peven, wPodd);
				result[k + i / 2] = csub(Peven, wPodd);
				w = cmul(w, step);
			}
		}
	}
	return result;
}

static dcomplex* _raw_radarfft(dcomplex* data, unsigned N, bool is_forward, bool keep_input) {

	if (N == 1) return data;

	unsigned M = nextPowerOf2(2 * (N - 1) - 1);
	unsigned g = findPrimitiveRoot(N);

	// a[q] = x_g^q = x[ pow(g, q) (mod N) ]
	dcomplex* a = 0;
	if ((a = (dcomplex*)calloc(M, sizeof(dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	// b[q] = x_g^-q = x[ pow(g, -q) (mod N) ]
	dcomplex* b = 0;
	if ((b = (dcomplex*)calloc(M, sizeof(dcomplex))) == NULL) {
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
	*b = cexp(cbuild(0., -2 * PI * k / N));
	for (unsigned q = 1; q < N - 1; q++) {
		*(a + q + M - N + 1) = *(data + expModuloN(g, q, N));
		*(b + q) = cexp(cbuild(0., -2 * PI * *(g_mod_minus_q + q - 1) * k / N));
		// pow(g, -2) (mod N)
		//      = { [pow(g, -1) (mod N)] * [pow(g, -1) (mod N)] } % N
		//		= [ pow(g, -1) * pow(g, -1) ] (mod N)
		*(g_mod_minus_q + q) = (*(g_mod_minus_q + q - 1) * *g_mod_minus_q) % N;
	}
	for (unsigned i = N - 1; i < M; i++) {
		*(b + i) = *(b + i - N + 1);
	}

	dcomplex* fft_a = _raw_radix2fft(a, M, true, false);
	dcomplex* fft_b = _raw_radix2fft(b, M, true, false);

	dcomplex* fft_ab = 0;
	if ((fft_ab = (dcomplex*)calloc(M, sizeof(dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned i = 0; i < M; i++) {
		*(fft_ab + i) = cmul(*(fft_a + i), *(fft_b + i));
	}
	free(fft_a);
	free(fft_b);

	dcomplex* conv_ab = _raw_radix2fft(fft_ab, M, false, false);

	dcomplex* result = 0;
	if ((result = (dcomplex*)calloc(N, sizeof(dcomplex))) == NULL) {
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

static dcomplex* _raw_bluesteinfft(dcomplex* data, unsigned N, bool is_forward, bool keep_input) {

	unsigned M = nextPowerOf2(2 * N - 1);

	dcomplex* a = 0;
	if ((a = (dcomplex*)calloc(M, sizeof(dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	dcomplex* b = 0;
	if ((b = (dcomplex*)calloc(M, sizeof(dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	int k = (is_forward) ? 1 : -1;
	double common = k * PI / N;
	for (unsigned n = 0; n < N; n++) {
		double twiddle = common * n * n;
		*(a + n) = cmul(*(data + n), cexp(cbuild(0, -twiddle)));
		*(b + n) = cexp(cbuild(0, twiddle));
	}
	if (!keep_input) free(data);
	for (unsigned n = M - N + 1; n < M; n++) {
		unsigned p = M - n;
		*(b + n) = cexp(cbuild(0, common * p * p));
	}

	dcomplex* fft_a = _raw_radix2fft(a, M, true, false);
	dcomplex* fft_b = _raw_radix2fft(b, M, true, false);

	dcomplex* fft_ab = 0;
	if ((fft_ab = (dcomplex*)calloc(M, sizeof(dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned i = 0; i < M; i++) {
		*(fft_ab + i) = cmul(*(fft_a + i), *(fft_b + i));
	}
	free(fft_a);
	free(fft_b);

	dcomplex* conv_ab = _raw_radix2fft(fft_ab, M, false, false);

	dcomplex* result = 0;
	if ((result = (dcomplex*)realloc(conv_ab, N * sizeof(dcomplex))) == NULL) {
		printf("Reallocation failed!\n");
		free(conv_ab);
		return 0;
	}
	conv_ab = NULL;

	for (unsigned n = 0; n < N; n++) {
		*(result + n) = _Cmulcr(cmul(*(result + n), cexp(cbuild(0, -common * n * n))), 1. / M);
	}

	return result;
}
