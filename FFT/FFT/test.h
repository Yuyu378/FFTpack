// 
// test.h
// 
//		Test header
// 
//		v0.2.0
//		----------------------------
//		* Implement function Compare() to compare execute time between 
//			~ 8-point butterfly diagram
//			~ Cooley-Tukey radix-2 (recursion)
//			~ Cooley-Tukey radix-2 (recursion) with 8-point butterfly diagram accelerate
//			~ Cooley-Tukey radix-2 with bit reversal
// 
//		v0.3.0
//		----------------------------
//		* Rewrite butterfly_8point_fft to Butterfly_8point_fft function
//			~ fix memory leak when used in recursion functions
//		* Rewrite CooleyTukey_Radix2_nonRecursion function
//			~ avoid memory leak when used in recursion functions
//		* Add Radar function to test Radar's Algorithm
//		* Add CooleyTukey_HybridRadix function to test Variation Cooley-Tukey Algorithm
//		* Add Bluestein function to test Bluestein's Algorithm
//		* Rewrite Compare function
//			~ Modularize the execution of the FFT transformation into a function named test
//			~ Add the comparison functions of Radar's/Variation Cooley-Tukey/Bluestein's algorithm
// 

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdbool.h>
#include <Windows.h>

#include "numeric.h"
#include "excomplex.h"

_Dcomplex* CooleyTukey_Radix2_Recursion(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);

_Dcomplex* Butterfly_8point_fft(_Dcomplex* data, bool is_forward, bool keep_input);

_Dcomplex* CooleyTukey_Radix2_Recursion_8pointAccelerate(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);

_Dcomplex* CooleyTukey_Radix2_nonRecursion(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);

_Dcomplex* Radar(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input, _Dcomplex* (*fft)());

_Dcomplex* CooleyTukey_HybridRadix(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input, _Dcomplex* (*fft)());

_Dcomplex* Bluestein(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input, _Dcomplex* (*fft)());

void test(_Dcomplex* x, unsigned N, bool is_display, _Dcomplex* (*basic_fft)(), _Dcomplex* (*target_fft)());

void Compare(unsigned N, bool is_display);
