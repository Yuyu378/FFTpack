// 
// fft.h
// 
//		Fourier Transform algorithms header
// 
//		v0.1.0
//		----------------------------
//		* Implement basic Discrete Fourier Transform algorithm
//		* Integration as dft algorithm
// 
//		* Implement 8-points butterfly diagram Decimation-in-Time algorithm
//		* Implement radix-2 Cooley-Tukey Fast Fourier Transform algorithm
//		* Implement Radar's Fast Fourier Transform algorithm
//		* Implement Mixed-radix (Hybrid-radix) Fast Fourier Transform algorithm
//		* Integration as fft algorithm
// 
//		v0.1.1
//		----------------------------
//		* Simple optimized Radar's algorithm
//		* Fix memory leaks caused by recursion
// 
//		v0.2.0
//		----------------------------
//		* Rewrite _raw_radix2fft by bit reversal method
//		* Deprecate _raw_radix8fft (see test.h & test.c)
//		* Deprecate radix8fft (see test.h & test.c)
//		* Deprecate radix8fft (see test.h & test.c)
// 
//		v0.3.0
//		----------------------------
//		* Implement _raw_bluesteinfft
//		* Deprecate (no use)
//			~ radix2fft	/ radix2ifft
//			~ radarfft / radarifft
//			~ hybrid_radixfft / hybrid_radixifft
//		* Deprecate (replaced)
//			~ _raw_hybrid_radixfft (replaced by _raw_bluesteinfft)
// 

#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdbool.h>

#include "numeric.h"
#include "excomplex.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Types
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

typedef enum {
	backward,
	ortho,
	forward
} norm_mode;

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Functions
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

static double _get_forward_norm(unsigned N, norm_mode mode);

static double _get_backward_norm(unsigned N, norm_mode mode);

static _Dcomplex* _execute_dft(_Dcomplex* data, unsigned int N, bool is_forward);

static _Dcomplex* _raw_dft(_Dcomplex* data, unsigned int N, bool is_forward, double fct);

/// <summary>
/// Calculate Decrete Fourier Transform
/// <para>.</para>
/// <para> N must be same as data size </para>
/// <para> mode have { backward / ortho / forward } </para>
/// </summary>
/// <param name="data  "> _Dcomplex* </param>
/// <param name="N       "> unsigned int </param>
/// <param name="mode"> norm_mode </param>
/// <returns> 
/// <para> mode = backward, 1 * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// <para> mode = ortho, (1 / sqrt(N)) * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// <para> mode = forward, (1 / N) * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// </returns>
_Dcomplex* dft(_Dcomplex* data, unsigned int N, norm_mode mode);

/// <summary>
/// Calculate Inverse Decrete Fourier Transform
/// <para>.</para>
/// <para> N must be same as data size </para>
/// <para> mode have { backward / ortho / forward } </para>
/// </summary>
/// <param name="data  "> _Dcomplex* </param>
/// <param name="N       "> unsigned int </param>
/// <param name="mode"> norm_mode </param>
/// <returns> 
/// <para> mode = backward, (1 / N) * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// <para> mode = ortho, (1 / sqrt(N)) * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// <para> mode = forward, 1 * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// </returns>
_Dcomplex* idft(_Dcomplex* data, unsigned int N, norm_mode mode);

static _Dcomplex* _execute_fft(_Dcomplex* data, unsigned int N, bool is_forward);

static _Dcomplex* _raw_fft(_Dcomplex* data, unsigned int N, bool is_forward, double fct);

/// <summary>
/// Calculate Fast Fourier Transform
/// <para>.</para>
/// <para> N must be same as data size </para>
/// <para> mode have { backward / ortho / forward } </para>
/// </summary>
/// <param name="data  "> _Dcomplex* </param>
/// <param name="N       "> unsigned int </param>
/// <param name="mode"> norm_mode </param>
/// <returns> 
/// <para> mode = backward, 1 * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// <para> mode = ortho, (1 / sqrt(N)) * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// <para> mode = forward, (1 / N) * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// </returns>
_Dcomplex* fft(_Dcomplex* data, unsigned int N, norm_mode mode);

/// <summary>
/// Calculate Inverse Fast Fourier Transform
/// <para>.</para>
/// <para> N must be same as data size </para>
/// <para> mode have { backward / ortho / forward } </para>
/// </summary>
/// <param name="data  "> _Dcomplex* </param>
/// <param name="N       "> unsigned int </param>
/// <param name="mode"> norm_mode </param>
/// <returns> 
/// <para> mode = backward, (1 / N) * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// <para> mode = ortho, (1 / sqrt(N)) * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// <para> mode = forward, 1 * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// </returns>
_Dcomplex* ifft(_Dcomplex* data, unsigned int N, norm_mode mode);

static _Dcomplex* _raw_radix2fft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);

static _Dcomplex* _raw_radarfft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);

static _Dcomplex* _raw_bluesteinfft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);
