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

#include <stdbool.h>

#include "complexlib.h"

dcomplex* CooleyTukey_Radix2_Recursion(dcomplex* data, unsigned N, bool is_forward, bool keep_input);

dcomplex* Butterfly_8point_fft(dcomplex* data, bool is_forward, bool keep_input);

dcomplex* CooleyTukey_Radix2_Recursion_8pointAccelerate(dcomplex* data, unsigned N, bool is_forward, bool keep_input);

dcomplex* CooleyTukey_Radix2_nonRecursion(dcomplex* data, unsigned N, bool is_forward, bool keep_input);

dcomplex* Radar(dcomplex* data, unsigned N, bool is_forward, bool keep_input, dcomplex* (*fft)());

dcomplex* CooleyTukey_HybridRadix(dcomplex* data, unsigned N, bool is_forward, bool keep_input, dcomplex* (*fft)());

dcomplex* Bluestein(dcomplex* data, unsigned N, bool is_forward, bool keep_input, dcomplex* (*fft)());

void test(dcomplex* x, unsigned N, bool is_display, dcomplex* (*basic_fft)(), dcomplex* (*target_fft)());

void Compare(unsigned N, bool is_display);
