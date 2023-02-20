// 
// excomplex.h
// 
//		Extend complex header
// 
//		v0.1.1
//		----------------------------
//		* Split the complex algorithms from mymath
//      * Add polar coordinate system types and algorithms
//      * Simple optimized self-defined complex algorithms
// 

#pragma once

#include <math.h>
#include <stdio.h>
#include <complex.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Types
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef _C_POLAR_T
#define _C_POLAR_T
typedef struct _C_double_polar {
    double _Val[2];
} _C_double_polar;

typedef struct _C_float_polar {
    float _Val[2];
} _C_float_polar;

typedef struct _C_ldouble_polar {
    long double _Val[2];
} _C_ldouble_polar;
#endif

typedef _C_double_polar  _Dpolar;
typedef _C_float_polar   _Fpolar;
typedef _C_ldouble_polar _Lpolar;

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Macros
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef PI
#define PI (3.1415926535897932384626433832795028841971693993751058209749445923)
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Functions
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

/* self-defined complex function */

void cprint(_Dcomplex a);

void cswap(_Dcomplex* a, _Dcomplex* b);

_Dcomplex cadd(_Dcomplex a, _Dcomplex b);

_Dcomplex csub(_Dcomplex a, _Dcomplex b);

_Dcomplex cmul(_Dcomplex a, _Dcomplex b);

_Dcomplex cdiv(_Dcomplex a, _Dcomplex b);

// polar to complex
_Dcomplex p2c(_Dpolar p);

/* polar coordinate system function */

_Dpolar _Pbuild(double _Radial, double _Angular);

_Dpolar padd(_Dpolar a, _Dpolar b);

_Dpolar psub(_Dpolar a, _Dpolar b);

_Dpolar pmul(_Dpolar a, _Dpolar b);

_Dpolar pdiv(_Dpolar a, _Dpolar b);

_Dpolar pexp(_Dpolar p, double n);

// complex to polar
_Dpolar c2p(_Dcomplex z);
