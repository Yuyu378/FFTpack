// 
// fft.h
// 
//      Fourier Transform algorithms header
// 
//      v0.1.0
//      ----------------------------
//      * Implement basic Discrete Fourier Transform algorithm
//      * Integration as dft algorithm
// 
//      * Implement 8-points butterfly diagram Decimation-in-Time algorithm
//      * Implement radix-2 Cooley-Tukey Fast Fourier Transform algorithm
//      * Implement Rader's Fast Fourier Transform algorithm
//      * Implement Mixed-radix (Hybrid-radix) Fast Fourier Transform algorithm
//      * Integration as fft algorithm
// 
//      v0.1.1
//      ----------------------------
//      * Simple optimized Rader's algorithm
//      * Fix memory leaks caused by recursion
// 
//      v0.2.0
//      ----------------------------
//      * Rewrite _cooley_tukey by bit reversal method
//      * Deprecate _raw_radix8fft (see test.h & test.c)
//      * Deprecate radix8fft (see test.h & test.c)
//      * Deprecate radix8fft (see test.h & test.c)
// 
//      v0.3.0
//      ----------------------------
//      * Implement _bluestein
//      * Deprecate (no use)
//          ~ radix2fft	/ radix2ifft
//          ~ raderfft / raderifft
//          ~ hybrid_radixfft / hybrid_radixifft
//      * Deprecate (replaced)
//          ~ _raw_hybrid_radixfft (replaced by _bluestein)
// 
//      v0.4.0
//      ----------------------------
//      * Implement memory pooling mechanism
//        (which can avoid memory management difficulties and fragmentation issues)
//      * Combine numeric's header and source into fft
//      * Rename functions to lower snake case
//      * Rewrite factor function using preallocated techniques
//      * Privatize the functions, making only dft/idft and fft/ifft public
//      * Rename all private functions to the form _xxx 
//        e.g. void foo(void) -> static void _foo(void)
// 

#pragma once

/* Includes ------------------------------------------------------------------------ */
#if defined(ARDUINO)
// and other C compilers with no <complex.h>
#include "complexlib.h"
#else
#include <complex.h>
#endif

/* Type Definitions ---------------------------------------------------------------- */

#ifdef __GNUC__
typedef double complex dcomplex;
#endif
#ifdef _MSC_VER
typedef _Dcomplex dcomplex;
#endif

typedef enum {
    backward,
    ortho,
    forward
} norm_mode;

/* Function Declarations ----------------------------------------------------------- */

/**
 * @brief Compute the Discrete Fourier Transform (DFT) on the input data.
 * 
 * @cite 
 * Wikipedia contributors. (2024, April 30). Discrete Fourier transform. 
 * In Wikipedia, The Free Encyclopedia. Retrieved 06:18, May 30, 2024, from
 * https://en.wikipedia.org/w/index.php?title=Discrete_Fourier_transform&oldid=1221523509
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.
 * @param norm_mode Mode to normalized.
 *                  `backward` - 1;
 *                  `ortho`    - 1/sqrt(N);
 *                  `forward`  - 1/N
 * @return Pointer to the transformed complex data array. `NULL` if memory allocation fails.
 * @note This function allocates memory for the result, which should be freed by the caller.
 */
dcomplex* dft(dcomplex* data, unsigned int N, norm_mode mode);

/**
 * @brief Compute the Inverse Discrete Fourier Transform (IDFT) on the input data.
 * 
 * @cite 
 * Wikipedia contributors. (2024, April 30). Discrete Fourier transform. 
 * In Wikipedia, The Free Encyclopedia. Retrieved 06:18, May 30, 2024, from
 * https://en.wikipedia.org/w/index.php?title=Discrete_Fourier_transform&oldid=1221523509
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.
 * @param norm_mode Mode to normalized.
 *                  `backward` - 1/N;
 *                  `ortho`    - 1/sqrt(N);
 *                  `forward`  - 1
 * @return Pointer to the transformed complex data array. `NULL` if memory allocation fails.
 * @note This function allocates memory for the result, which should be freed by the caller.
 */
dcomplex* idft(dcomplex* data, unsigned int N, norm_mode mode);

/**
 * @brief Compute the Fast Fourier Transform (FFT) on the input data.
 * 
 * @cite 
 * Wikipedia contributors. (2024, May 1). Fast Fourier transform. 
 * In Wikipedia, The Free Encyclopedia. Retrieved 20:28, May 30, 2024, from 
 * https://en.wikipedia.org/w/index.php?title=Fast_Fourier_transform&oldid=1221693816
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.
 * @param norm_mode Mode to normalized.
 *                  `backward` - 1;
 *                  `ortho`    - 1/sqrt(N);
 *                  `forward`  - 1/N
 * @return Pointer to the transformed complex data array. `NULL` if memory allocation fails.
 * @note This function allocates memory for the result, which should be freed by the caller.
 */
dcomplex* fft(dcomplex* data, unsigned int N, norm_mode mode);

/**
 * @brief Compute the Inverse Fast Fourier Transform (IFFT) on the input data.
 * 
 * @cite 
 * Wikipedia contributors. (2024, May 1). Fast Fourier transform. 
 * In Wikipedia, The Free Encyclopedia. Retrieved 20:28, May 30, 2024, from 
 * https://en.wikipedia.org/w/index.php?title=Fast_Fourier_transform&oldid=1221693816
 *
 * @param data Pointer to the input complex data array.
 * @param N Size of the input data array.
 * @param norm_mode Mode to normalized.
 *  				`backward` - 1/N;
 *                  `ortho`    - 1/sqrt(N);
 *                  `forward`  - 1
 * @return Pointer to the transformed complex data array. `NULL` if memory allocation fails.
 * @note This function allocates memory for the result, which should be freed by the caller.
 */
dcomplex* ifft(dcomplex* data, unsigned int N, norm_mode mode);
