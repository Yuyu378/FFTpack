// 
// fft.c
// 
//		Fourier Transform algorithms implementation
// 

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "fft.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#if defined(M_PI) && !defined(PI)
#define PI M_PI
#endif

#ifndef PI
#define PI (3.1415926535897932384626433832795028841971693993751058209749445923)
#endif

static double _get_forward_norm(unsigned N, norm_mode mode);
static double _get_backward_norm(unsigned N, norm_mode mode);
static dcomplex* _execute_dft(dcomplex* data, unsigned int N, _Bool is_forward);
static dcomplex* _raw_dft(dcomplex* data, unsigned int N, _Bool is_forward, double fct);
static _Bool _isprime(unsigned num);
static _Bool _is_power_of_2(unsigned num);
static dcomplex* _execute_fft(dcomplex* data, unsigned int N, _Bool is_forward);
static dcomplex* _raw_fft(dcomplex* data, unsigned int N, _Bool is_forward, double fct);
static void _bit_reverse_sorting(dcomplex* arr, unsigned len, dcomplex** pres);
static dcomplex* _cooley_tukey(dcomplex* data, unsigned N, _Bool is_forward);
static unsigned _next_power_of_2(unsigned num);
static unsigned* _factor(unsigned num);
static int _ex_euclidean(unsigned a, unsigned b, int* x, int* y);
static int* _exgcd(unsigned a, unsigned b);
static unsigned _expmodn(unsigned g, unsigned k, unsigned n);
static unsigned _iexpmodn(unsigned a, unsigned n);
static _Bool _is_primitive_root(unsigned g, unsigned p);
static unsigned _find_primitive_root(unsigned n);
static dcomplex* _rader(dcomplex* data, unsigned N, _Bool is_forward);
static dcomplex* _bluestein(dcomplex* data, unsigned N, _Bool is_forward);

/**
 * @brief Computes the normalization factor for the forward transform.
 *
 * @param N The size of the input data array.
 * @param mode The normalization mode.
 * @return The normalization factor.
 */
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

/**
 * @brief Computes the normalization factor for the backward (inverse) transform.
 *
 * @param N The size of the input data array.
 * @param mode The normalization mode.
 * @return The normalization factor.
 */
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

/**
 * @brief Executes the Discrete Fourier Transform (DFT) on the input data.
 *        This is a core function that realizes the theoretical algorithm.
 * @cite 
 * Wikipedia contributors. (2024, April 30). Discrete Fourier transform. 
 * In Wikipedia, The Free Encyclopedia. Retrieved 06:18, May 30, 2024, from
 * https://en.wikipedia.org/w/index.php?title=Discrete_Fourier_transform&oldid=1221523509
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.
 * @param is_forward Boolean indicating the direction of the transform.
 *                   True for forward transform, False for inverse transform.
 * @return Pointer to the transformed complex data array. 
 *         `NULL` if memory allocation fails.
 * @note This function allocates memory for the result, which should be freed by the caller.
 */
static dcomplex* _execute_dft(dcomplex* data, unsigned int N, _Bool is_forward) {
    dcomplex* res = NULL;
    res = (dcomplex*)calloc(N, sizeof(dcomplex));
    if (res == NULL) {
        printf("Allocation failed.\n");
        return NULL;
    }
    int b = (is_forward) ? 1 : -1;
    for (unsigned k = 0; k < N; k++) {
        for (unsigned n = 0; n < N; n++) {
            res[k] += data[n] * cexp(-2.0 * I * PI * k * n * b / N);
        }
    }
    return res;
}

/**
 * @brief Executes the raw Discrete Fourier Transform (DFT) and applies a scaling factor.
 *        This is a transition function that adjusts the output of the core function.
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.
 * @param is_forward Boolean indicating the direction of the transform.
 *                   True for forward transform, False for inverse transform.
 * @param fct The scaling factor to be applied to the result.
 * @return Pointer to the transformed complex data array. 
 *         `NULL` if memory allocation fails.
 * @note This function allocates memory for the result, which should be freed by the caller.
 */
static dcomplex* _raw_dft(dcomplex* data, unsigned int N, _Bool is_forward, double fct){
    dcomplex* res = _execute_dft(data, N, is_forward);
    for (unsigned i = 0; i < N; i++) {
        res[i] *= fct;
    }
    return res;
}

dcomplex* dft(dcomplex* data, unsigned int N, norm_mode mode) {
    double fct = 1. / _get_forward_norm(N, mode);
    return _raw_dft(data, N, 1, fct);
}

dcomplex* idft(dcomplex* data, unsigned int N, norm_mode mode) {
    double fct = 1. / _get_backward_norm(N, mode);
    return _raw_dft(data, N, 0, fct);
}

/**
 * @brief Determines if a given number is prime.
 *
 * @param num The number to check.
 * @return True if the number is prime, False otherwise.
 */
static _Bool _isprime(unsigned num) {
    if (num <= 1) return 0;
    if (num == 2) return 1;
    if (num % 2 == 0) return 0;
    for (unsigned i = 3; i * i <= num; i += 2) {
        if (num % i == 0) return 0;
    }
    return 1;
}

/**
 * @brief Determines if a given number is a power of 2.
 *
 * @param num The number to check.
 * @return True if the number is a power of 2, False otherwise.
 */
static _Bool _is_power_of_2(unsigned num) {
    if (num && !(num & (num - 1))) return 1;
    return 0;
}

/**
 * @brief Executes the Fast Fourier Transform (FFT) on the input data.
 *        This is a transition function that calls the corresponding core function for 
 *        different data lengths.
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.
 * @param is_forward Boolean indicating the direction of the transform.
 *                   True for forward transform, False for inverse transform.
 * @return Pointer to the transformed complex data array. 
 * @note The returned pointer should be freed by the caller.
 */
static dcomplex* _execute_fft(dcomplex* data, unsigned int N, _Bool is_forward) {
    if (N == 1) {
        dcomplex* ans = NULL;
        ans = (dcomplex*)calloc(N, sizeof(dcomplex));
        if (ans == NULL) {
            printf("Allocation failed.\n");
            return NULL;
        }
        ans[0] = data[0];
        return ans;
    }
    if (_is_power_of_2(N)) return _cooley_tukey(data, N, is_forward);
    if (_isprime(N)) return _rader(data, N, is_forward);
    return _bluestein(data, N, is_forward);
}

/**
 * @brief Executes the raw Fast Fourier Transform (FFT) and applies a scaling factor.
 *        This is a transition function that adjusts the output of the core function.
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.
 * @param is_forward Boolean indicating the direction of the transform.
 *                   True for forward transform, False for inverse transform.
 * @param fct The scaling factor to be applied to the result.
 * @return Pointer to the transformed complex data array. `NULL` if memory allocation fails.
 *
 * @note This function allocates memory for the result, which should be freed by the caller.
 */
static dcomplex* _raw_fft(dcomplex* data, unsigned int N, _Bool is_forward, double fct) {
    dcomplex* res = _execute_fft(data, N, is_forward);
    for (unsigned i = 0; i < N; i++) {
        res[i] *= fct;
    }
    return res;
}

dcomplex* fft(dcomplex* data, unsigned int N, norm_mode mode) {
    double fct = 1. / _get_forward_norm(N, mode);
    return _raw_fft(data, N, 1, fct);
}

dcomplex* ifft(dcomplex* data, unsigned int N, norm_mode mode) {
    double fct = 1. / _get_backward_norm(N, mode);
    dcomplex* ans = _raw_fft(data, N, 0, fct);
    for (unsigned int i = 0; i < N; i++) {
        if (fabs(cimag(ans[i])) < __DBL_EPSILON__) {
            ans[i] = creal(ans[i]);
        }
    }
    return ans;
}

/**
 * @brief Performs bit-reverse sorting on the input data array.
 *
 * @param arr Pointer to the input complex data array.
 * @param len Size of the input data array.❗Must be power of 2.
 * @param res Pointer to the output result complex data array.❗Length must be same as `arr`.
 * @return None
 */
static void _bit_reverse_sorting(dcomplex* arr, unsigned len, dcomplex** pres) {
    dcomplex* res = *pres;
    memcpy_s(res, len * sizeof(dcomplex), arr, len * sizeof(dcomplex));
    unsigned j = len >> 1;
    for (unsigned i = 1; i < len - 1; i++) {
        if (i < j) {
            dcomplex tmp = res[i];
            res[i] = res[j];
            res[j] = tmp;
        }
        unsigned k = len >> 1;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    return;
}

/**
 * @brief Implements the Cooley-Tukey FFT algorithm.
 *        This is a core function that realizes the theoretical algorithm.
 * @cite
 * Wikipedia contributors. (2024, February 21). Cooley–Tukey FFT algorithm. 
 * In Wikipedia, The Free Encyclopedia. Retrieved 06:25, May 30, 2024, from 
 * https://en.wikipedia.org/w/index.php?title=Cooley%E2%80%93Tukey_FFT_algorithm&oldid=1209251205
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.❗Must be power of 2.
 * @param is_forward Boolean indicating the direction of the transform.
 *                   True for forward transform, False for inverse transform.
 * @param pres Pointer of pointer to the output result complex data array.
 * @return None
 */
static void _cooley_tukey_core(dcomplex* data, unsigned N, _Bool is_forward, dcomplex** pres) {
    _bit_reverse_sorting(data, N, pres);
    dcomplex* res = *pres;
    int b = (is_forward) ? 1 : -1;
    for (unsigned i = 2; i <= N; i <<= 1) {
        dcomplex step = cexp(-2.0 * I * PI * b / i);
        for (unsigned j = 0; j < N; j += i) {
            dcomplex w = 1;
            for (unsigned k = j; k < j + i / 2; k++) {
                dcomplex Peven = res[k];
                dcomplex wPodd = w * res[k + i / 2];
                res[k] = Peven + wPodd;
                res[k + i / 2] = Peven - wPodd;
                w *= step;
            }
        }
    }
    return;
}

/**
 * @brief Cooley-Tukey FFT transform.
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.❗Must be power of 2.
 * @param is_forward Boolean indicating the direction of the transform.
 *                   True for forward transform, False for inverse transform.
 * @return Pointer to the transformed complex data array. NULL if memory allocation fails.
 *
 * @note This function allocates memory for the result if `res` variable input `NULL`, 
 *       which should be freed by the caller.
 */
static dcomplex* _cooley_tukey(dcomplex* data, unsigned N, _Bool is_forward) {
    dcomplex* res = (dcomplex*)calloc(N, sizeof(dcomplex));
    if (res == NULL) {
        printf("Allocation failed.\n");
        return NULL;
    }
    _cooley_tukey_core(data, N, is_forward, &res);
    return res;
}

/**
 * @brief Finds the next power of 2 greater than or equal to a given number.
 *
 * @param num The input number.
 * @return The next power of 2.
 */
static unsigned _next_power_of_2(unsigned num) {
    unsigned val = 1;
    if (_is_power_of_2(num)) return num;
    while (val < num) val <<= 1;
    return val;
}

/**
 * @brief Factors a given number into its prime factors.
 *
 * @param num The number to be factored.
 * @return An array of prime factors, terminated by a zero. NULL if memory allocation fails.
 */
static unsigned* _factor(unsigned num) {
    // Handle special cases 0 and 1
    if (num == 0 || num == 1) {
        unsigned* factors = (unsigned*)calloc(2, sizeof(unsigned));
        if (factors == NULL) {
            printf("Allocation failed.\n");
            return NULL;
        }
        factors[0] = 2; // Two elements: count and the factor itself
        factors[1] = num;
        return factors;
    }

    // Estimate the initial size of the array to minimize reallocations
    unsigned estimated_size = (unsigned)(log2(num) + 1) + 1; // +1 for storing the count
    unsigned* factors = (unsigned*)calloc(estimated_size, sizeof(unsigned));
    if (factors == NULL) {
        printf("Allocation failed.\n");
        return NULL;
    }
    factors[0] = 1; // Start with count 1 to include the count storage itself

    unsigned count = 1; // Initial count (1 because factors[0] is reserved for count)
    unsigned limit = (unsigned)sqrt(num) + 1;

    for (unsigned i = 2; i <= limit; i++) {
        while (num % i == 0) {
            if (count >= estimated_size) {
                estimated_size *= 2;
                factors = realloc(factors, estimated_size * sizeof(unsigned));
                if (factors == NULL) {
                    printf("Allocation failed.\n");
                    return NULL;
                }
            }
            factors[count++] = i;
            num /= i;
        }
    }

    if (num > 1) {
        if (count >= estimated_size) {
            estimated_size += 1;
            factors = realloc(factors, estimated_size * sizeof(unsigned));
            if (!factors) {
                printf("Allocation failed.\n");
                return NULL;
            }
        }
        factors[count++] = num;
    }

    // Reallocate memory to fit the actual number of factors found
    factors = realloc(factors, count * sizeof(unsigned));
    if (factors == NULL) {
        printf("Allocation failed.\n");
        return NULL;
    }
    
    factors[0] = count; // Store the number of factors found
    return factors;
}

/**
 * @brief Performs the extended Euclidean algorithm to find the greatest common divisor (GCD) 
 *        and the coefficients of Bézout's identity.
 *
 * @param a First integer.
 * @param b Second integer.
 * @param x Pointer to store the coefficient for `a`.
 * @param y Pointer to store the coefficient for `b`.
 * @return The GCD of `a` and `b`.
 */
static int _ex_euclidean(unsigned a, unsigned b, int* x, int* y) {
    if (b == 0) {
        *x = 1, * y = 0;
        return a;
    }
    int d = _ex_euclidean(b, a % b, y, x);
    *y -= a / b * *x;
    return d;
}

/**
 * @brief Finds the coefficients of Bézout's identity for two integers using the extended 
 *        Euclidean algorithm.
 *
 * @param a First integer.
 * @param b Second integer.
 * @return An array containing the GCD and the coefficients. NULL if memory allocation fails.
 */
static int* _exgcd(unsigned a, unsigned b) {
    int x = 0, y = 0, g = 0;
    g = _ex_euclidean(a, b, &x, &y);

    int* r = (int*)calloc(4, sizeof(int));
    if (r == NULL) {
        printf("Allocation failed.\n");
        return 0;
    }
    r[0] = 4, r[1] = x, r[2] = y, r[3] = g;
    return r;
}

/**
 * @brief Computes `g^k mod n` using modular exponentiation.
 *
 * @param g The base.
 * @param k The exponent.
 * @param n The modulus.
 * @return The result of `g^k mod n`.
 */
static unsigned _expmodn(unsigned g, unsigned k, unsigned n) {
    unsigned a = 1;
    while (k) {
        if (k & 1) a = (a * g) % n;
        k >>= 1;
        g = (g * g) % n;
    }
    return a;
}

/**
 * @brief Finds the modular multiplicative inverse of `a` modulo `n`.
 *
 * @param a The number to find the inverse of.
 * @param n The modulus.
 * @return The modular multiplicative inverse of `a` modulo `n`.
 */
static unsigned _iexpmodn(unsigned a, unsigned n) {
    int* arr = _exgcd(a, n);
    int x = 0;
    if (arr[3] == 1) {
        if (arr[1] > 0) {
            x = arr[1];
            free(arr);
            return x;
        }
        else {
            x = arr[1] + n;
            free(arr);
            return x;
        }
    }
    return 0;
}

/**
 * @brief Checks if a given number is a primitive root modulo `p`.
 *
 * @param g The number to check.
 * @param p The modulus.
 * @return True if `g` is a primitive root modulo `p`, False otherwise.
 */
static _Bool _is_primitive_root(unsigned g, unsigned p) {
    unsigned* r = _factor(p - 1);
    if (r == NULL) return 0;
    for (unsigned i = 1; i < r[0]; i++) {
        if (_expmodn(g, (p - 1) / r[i], p) == 1) {
            free(r);
            return 0;
        }
    }
    free(r);
    return 1;
}

/**
 * @brief Finds a primitive root modulo `n`.
 * 
 * This function finds a primitive root of a prime number `n`. A primitive root `g` of `n` is 
 * a number such that for every integer `k` from 1 to `n-1`, there exists an integer `m` such 
 * that `g^m ≡ k(mod n)`. In other words, the powers of `g` generate all the integers 
 * from 1 to `n-1` modulo `n`.
 *
 * @param n The prime number for which to find the primitive root.
 * @return The primitive root of `n`, or 0 if `n` is not prime or if a primitive root is not found.
 *         Note: The function prints an error message if `n` is not prime or if a primitive root is not found.
 *
 * @note This function assumes that the `_isprime` and `_is_primitive_root` helper functions 
 *       are defined elsewhere in the code.
 */
static unsigned _find_primitive_root(unsigned n) {
    if (!_isprime(n)) {
        printf("n should be prime!\n");
        return 0;
    }
    for (unsigned rootCandidate = 2; rootCandidate < n; rootCandidate++) {
        if (_is_primitive_root(rootCandidate, n)) return rootCandidate;
    }
    printf("Primitive root not found");
    return 0;
}

/**
 * @brief Implements the Rader's FFT algorithm.
 *        This is a core function that realizes the theoretical algorithm.
 * @cite
 * Wikipedia contributors. (2022, July 12). Rader's FFT algorithm. 
 * In Wikipedia, The Free Encyclopedia. Retrieved 06:35, May 30, 2024, from 
 * https://en.wikipedia.org/w/index.php?title=Rader%27s_FFT_algorithm&oldid=1097739627
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.❗Must be a prime number.
 * @param is_forward Boolean indicating the direction of the transform.
 *                   True for forward transform, False for inverse transform.
 * @return Pointer to the transformed complex data array. NULL if memory allocation fails.
 *
 * @note This function allocates memory for the result, which should be freed by the caller.
 */
static dcomplex* _rader(dcomplex* data, unsigned N, _Bool is_forward) {

    unsigned M = _next_power_of_2(2 * (N - 1) - 1);
    unsigned g = _find_primitive_root(N);
    size_t total_size = 3 * M * sizeof(dcomplex) + (N - 1) * sizeof(unsigned);

    void* mem_pool = calloc(1, total_size);
    if (mem_pool == NULL) {
        printf("Allocation failed.\n");
        return NULL;
    }

    // a[q] = x_g^q = x[ pow(g, q) (mod N) ]
    dcomplex* a = (dcomplex*)mem_pool;
    // b[q] = x_g^-q = x[ pow(g, -q) (mod N) ]
    dcomplex* b = (dcomplex*)mem_pool + M;
    dcomplex* fft_a = (dcomplex*)mem_pool + 2 * M;
    unsigned* g_mod_minus_q = (unsigned*)((dcomplex*)mem_pool + 3 * M);

    int k = (is_forward) ? 1 : -1;
    *g_mod_minus_q = _iexpmodn(g, N);
    *a = data[_expmodn(g, 0, N)];
    *b = cexp(-2.0 * I * PI * k / N);
    for (unsigned q = 1; q < N - 1; q++) {
        a[q + M - N + 1] = data[_expmodn(g, q, N)];
        b[q] = cexp(-2.0 * I * PI * g_mod_minus_q[q - 1] * k / N);
        // pow(g, -2) (mod N)
        //      = { [pow(g, -1) (mod N)] * [pow(g, -1) (mod N)] } % N
        //		= [ pow(g, -1) * pow(g, -1) ] (mod N)
        g_mod_minus_q[q] = (g_mod_minus_q[q - 1] * *g_mod_minus_q) % N;
    }
    for (unsigned i = N - 1; i < M; i++) {
        b[i] = b[i - N + 1];
    }

    _cooley_tukey_core(a, M, 1, &fft_a);
    dcomplex* fft_b = (dcomplex*)mem_pool;
    if (fft_a == NULL) {
        printf(
            "Runtime error -\n"
            "Cooley Tukey transform failed when execute _rader()"
            "\n"
            "    _cooley_tukey_core(a, M, 1, &fft_a);\n"
            "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n"
            "\n"
        );
        free(mem_pool);
        return NULL;
    }
    
    _cooley_tukey_core(b, M, 1, &fft_b);
    dcomplex* fft_ab = (dcomplex*)mem_pool + M;
    if (fft_b == NULL) {
        printf(
            "Runtime error -\n"
            "Cooley Tukey transform failed when execute _rader()"
            "\n"
            "    _cooley_tukey_core(b, M, 1, &fft_b);\n"
            "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n"
            "\n"
        );
        free(mem_pool);
        return NULL;
    }

    for (unsigned i = 0; i < M; i++) {
        fft_ab[i] = fft_a[i] * fft_b[i];
    }
    dcomplex* conv_ab = (dcomplex*)mem_pool;

    _cooley_tukey_core(fft_ab, M, 0, &conv_ab);
    if (conv_ab == NULL) {
        printf(
            "Runtime error - \n"
            "Cooley Tukey transform failed when execute _rader()"
            "\n"
            "    _cooley_tukey_core(fft_ab, M, 0, &conv_ab);\n"
            "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n"
            "\n"
        );
        free(mem_pool);
        return NULL;
    }

    dcomplex* res = NULL;
    res = (dcomplex*)calloc(N, sizeof(dcomplex));
    if (res == NULL) {
        printf("Allocation failed.\n");
        free(mem_pool);
        return NULL;
    }
    for (unsigned n = 0; n < N; n++) {
        *res = *res + data[n];
    }
    res[1] = *data + *conv_ab / M;
    for (unsigned q = 1; q < N - 1; q++) {
        res[g_mod_minus_q[q - 1]] = *data + conv_ab[q] / M;
    }
    free(mem_pool);
    return res;
}

/**
 * @brief Implements the Bluestein's FFT algorithm.
 *        This is a core function that realizes the theoretical algorithm.
 * @cite
 * Wikipedia contributors. (2023, August 18). Chirp Z-transform. 
 * In Wikipedia, The Free Encyclopedia. Retrieved 06:37, May 30, 2024, from 
 * https://en.wikipedia.org/w/index.php?title=Chirp_Z-transform&oldid=1171007425
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array. 
 *          (No size limitation, but faster by using Cooley-Tukey's algorithm when size is power of 2 
 *          and by Rader's when size is prime)
 * @param is_forward Boolean indicating the direction of the transform.
 *                   True for forward transform, False for inverse transform.
 * @param keep_input Boolean indicating whether to keep the input data array.
 * @return Pointer to the transformed complex data array. NULL if memory allocation fails.
 *
 * @note This function allocates memory for the result, which should be freed by the caller.
 */
static dcomplex* _bluestein(dcomplex* data, unsigned N, _Bool is_forward) {
    unsigned M = _next_power_of_2(2 * N - 1);
    size_t total_size = 3 * M * sizeof(dcomplex); // Space for a, b, and etc.

    dcomplex* mem_pool = (dcomplex*)calloc(1, total_size);
    if (mem_pool == NULL) {
        printf("Allocation failed.\n");
        return NULL;
    }

    dcomplex* a = mem_pool;
    dcomplex* b = mem_pool + M;
    dcomplex* fft_a = mem_pool + 2 * M;

    int k = (is_forward) ? 1 : -1;
    double common = k * PI / N;
    for (unsigned n = 0; n < N; n++) {
        double twiddle = common * n * n;
        a[n] = data[n] * cexp(-twiddle * I);
        b[n] = cexp(twiddle * I);
    }
    for (unsigned n = M - N + 1; n < M; n++) {
        unsigned p = M - n;
        b[n] = cexp(common * I * p * p);
    }

    _cooley_tukey_core(a, M, 1, &fft_a);
    dcomplex* fft_b = mem_pool;
    if (fft_a == NULL) {
        printf(
            "Runtime error -\n"
            "Cooley Tukey transform failed when executing _bluestein()"
            "\n"
            "    _cooley_tukey_core(a, M, 1, &fft_a);\n"
            "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n"
            "\n"
        );
        free(mem_pool);
        return NULL;
    }

    _cooley_tukey_core(b, M, 1, &fft_b);
    dcomplex* fft_ab = mem_pool + M;
    if (fft_b == NULL) {
        printf(
            "Runtime error -\n"
            "Cooley Tukey transform failed when executing _bluestein()"
            "\n"
            "    _cooley_tukey_core(b, M, 1, &fft_b);\n"
            "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n"
            "\n"
        );
        free(mem_pool);
        return NULL;
    }

    for (unsigned i = 0; i < M; i++) {
        fft_ab[i] = fft_a[i] * fft_b[i];
    }
    dcomplex* conv_ab = mem_pool;

    _cooley_tukey_core(fft_ab, M, 0, &conv_ab);
    if (conv_ab == NULL) {
        printf(
            "Runtime error - \n"
            "Cooley Tukey transform failed when executing _bluestein()"
            "\n"
            "    _cooley_tukey_core(fft_ab, M, 0, &conv_ab);\n"
            "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n"
            "\n"
        );
        free(mem_pool);
        return NULL;
    }

    dcomplex* res = (dcomplex*)realloc(conv_ab, N * sizeof(dcomplex));
    if (res == NULL) {
        printf("Reallocation failed!\n");
        free(mem_pool);
        return NULL;
    }

    for (unsigned n = 0; n < N; n++) {
        res[n] *= cexp(-common * I * n * n) / M;
    }

    return res;
}
