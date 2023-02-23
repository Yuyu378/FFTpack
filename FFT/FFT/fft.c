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
	for (int k = 0; k < (signed)N; k++) {
		for (int n = 0; n < (signed)N; n++) {
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
	else return _raw_hybrid_radixfft(data, N, is_forward, true);
}

static _Dcomplex* _raw_fft(_Dcomplex* data, unsigned int N, bool is_forward, double fct) {
	_Dcomplex* result = _execute_fft(data, N, is_forward);
	for (unsigned i = 0; i < N; i++) {
		*(result + i) = cmul(*(result + i), _Cbuild(fct, 0.));
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

static _Dcomplex* radix2fft(_Dcomplex* data, unsigned N) {
	return _raw_radix2fft(data, N, true, true);
}

static _Dcomplex* radix2ifft(_Dcomplex* data, unsigned N) {
	return _raw_radix2fft(data, N, false, true);
}

static _Dcomplex* _raw_radarfft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input) {

	if (N == 1) return data;

	unsigned M = nextPowerOf2(2 * (N - 1) - 1);
	unsigned g = findPrimitiveRoot(N);

	// a[q] = x_g^q = x[ pow(g, q) (mod N) ]
	_Dcomplex* a = 0;
	if ((a = (_Dcomplex*)calloc((size_t)N - 1, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	// b[q] = x_g^-q = x[ pow(g, -q) (mod N) ]
	_Dcomplex* b = 0;
	if ((b = (_Dcomplex*)calloc((size_t)N - 1, sizeof(_Dcomplex))) == NULL) {
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
		*(a + q) = *(data + expModuloN(g, q, N));
		*(b + q) = cexp(_Cbuild(0., -2 * PI * *(g_mod_minus_q + q - 1) * k / N));
		// pow(g, -2) (mod N)
		//      = { [pow(g, -1) (mod N)] * [pow(g, -1) (mod N)] } % N
		//		= [ pow(g, -1) * pow(g, -1) ] (mod N)
		*(g_mod_minus_q + q) = (*(g_mod_minus_q + q - 1) * *g_mod_minus_q) % N;
	}

	// Padding zero between the 0th element and the 1st element 
	// until the length of a is extended to M.
	_Dcomplex* new_a = 0;
	if ((new_a = (_Dcomplex*)realloc(a, M * sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		free(a);
		return 0;
	}
	a = NULL;
	// Repeat the array until the length increases to M
	_Dcomplex* new_b = 0;
	if ((new_b = (_Dcomplex*)realloc(b, M * sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		free(b);
		return 0;
	}
	b = NULL;

	for (unsigned i = 1; i < N - 1; i++) {
		*(new_a + i + M - N + 1) = _Cbuild(0., 0.);
		cswap(new_a + i, new_a + i + M - N + 1);
	}
	for (unsigned i = N - 1; i < M; i++) {
		if(i < (M - N + 1) + 1) *(new_a + i) = _Cbuild(0., 0.);
		*(new_b + i) = *(new_b + i - N + 1);
	}

	_Dcomplex* fft_a = radix2fft(new_a, M);
	free(new_a);
	_Dcomplex* fft_b = radix2fft(new_b, M);
	free(new_b);

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

	_Dcomplex* conv_ab = radix2ifft(fft_ab, M);
	free(fft_ab);

	_Dcomplex* result = 0;
	if ((result = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned n = 0; n < N; n++) {
		*result = cadd(*result, *(data + n));
	}
	*(result + 1) = cadd(*data, cdiv(*conv_ab, _Cbuild((double)M, 0.)));
	for (unsigned q = 1; q < N - 1; q++) {
		*(result + *(g_mod_minus_q + q - 1)) = cadd(*data, cdiv(*(conv_ab + q), _Cbuild((double)M, 0.)));
	}
	if (!keep_input) free(data);
	free(g_mod_minus_q);
	free(conv_ab);
	return result;
}

static _Dcomplex* radarfft(_Dcomplex* data, unsigned N) {
	return _raw_radarfft(data, N, true, true);
}

static _Dcomplex* radarifft(_Dcomplex* data, unsigned N) {
	return _raw_radarfft(data, N, false, true);
}

static _Dcomplex* _raw_hybrid_radixfft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input) {

	if (N == 1) return data;
	if (N == 2) return _raw_radix2fft(data, N, is_forward, false);
	if (isprime(N)) return _raw_radarfft(data, N, is_forward, false);
	if (isPowerOf2(N)) return _raw_radix2fft(data, N, is_forward, false);

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
		_Dcomplex* fft_inner = _raw_hybrid_radixfft(inner, N2, is_forward, false);
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
		_Dcomplex* fft_outer = _raw_hybrid_radixfft(outer, N1, is_forward, false);
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

static _Dcomplex* hybrid_radixfft(_Dcomplex* data, unsigned N) {
	return _raw_hybrid_radixfft(data, N, true, true);
}

static _Dcomplex* hybrid_radixifft(_Dcomplex* data, unsigned N) {
	return _raw_hybrid_radixfft(data, N, false, true);
}
