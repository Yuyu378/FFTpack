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

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdbool.h>
#include <Windows.h>

#include "numeric.h"
#include "excomplex.h"

_Dcomplex* butterfly_8point_fft(_Dcomplex* data, bool is_forward);

_Dcomplex* CooleyTukey_Radix2_Recursion(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);

_Dcomplex* CooleyTukey_Radix2_Recursion_8pointAccelerate(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);

_Dcomplex* CooleyTukey_Radix2_nonRecursion(_Dcomplex* data, unsigned N, bool is_forward);

void Compare(unsigned N);


