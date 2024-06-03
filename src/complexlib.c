//
// complexlib.c
//
//      complex library
//

#if defined(ARDUINO) // and other C compilers with no <complex.h>

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "complexlib.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#if defined(M_PI) && !defined(PI)
#define PI M_PI
#endif

#ifndef PI
#define PI (3.1415926535897932384626433832795028841971693993751058209749445923)
#endif

/* double complex type ------------------------------ */

int cisequal(dcomplex a, dcomplex b) {
    if (a.re ==  b.re && a.im ==  b.im) return  1;
    if (a.re == -b.re && a.im == -b.im) return -1;
    return 0;
}

double creal(dcomplex z) {
    return z.re;
}

double cimag(dcomplex z) {
    return z.im;
}

double carg(dcomplex z) {
    return atan2(z.im, z.re);
}

double cabs(dcomplex z) {
    return sqrt(z.re * z.re + z.im * z.im);
}

double norm(dcomplex z) {
    return z.re * z.re + z.im * z.im;
}

dcomplex cminus(dcomplex z) {
    return cbuild(-z.re, -z.im);
}

dcomplex conj(dcomplex z) {
    dcomplex r = {z.re, -z.im};
    return r;
}

dcomplex cacos(dcomplex z) {
    dcomplex tmp = _Csubrc(1., _Cpowcr(z, 2.)); // 1 - square(z)
    dcomplex isqrt_abs_tmp = _Cmulcr(cbuild(0., 1.), sqrt(cabs(tmp)));
    double theta = carg(tmp) / 2.;
    dcomplex e_iphi = cbuild(cos(theta), sin(theta));
    dcomplex ln_longform = clog(_Caddcc(z, _Cmulcc(isqrt_abs_tmp, e_iphi)));
    return _Cmulcc(cbuild(0., -1.), ln_longform);
}

dcomplex casin(dcomplex z) {
    dcomplex iz = cbuild(-z.im, z.re);
    dcomplex tmp = _Csubrc(1., _Cpowcr(z, 2.)); // 1 - square(z)
    double sqrt_abs_tmp = sqrt(cabs(tmp));
    double theta = carg(tmp) / 2.;
    dcomplex e_iphi = cbuild(cos(theta), sin(theta));
    dcomplex ln_longform = clog(_Caddcc(iz, _Cmulrc(sqrt_abs_tmp, e_iphi)));
    return _Cmulcc(cbuild(0., -1.), ln_longform);
}

dcomplex catan(dcomplex z) {
    dcomplex up = cbuild(-z.re, 1. - z.im);
    dcomplex down = cbuild(z.re, 1. + z.im);
    dcomplex post = clog(_Cdivcc(up, down));
    dcomplex ans = _Cmulcc(cbuild(0., -0.5), post);
    if (z.re == 0) {
        if (z.im == 1) return cbuild(0, INFINITY);
        if (z.im == -1) return cbuild(0, -INFINITY);
        if (z.im > 0) return _Caddcr(ans, PI);
        return ans;
    }
}

dcomplex ccos(dcomplex z) {
    if (z.re == 0) return cbuild(cosh(z.im), 0);
    if (z.im == 0) return cbuild(cos(z.re), 0);
    dcomplex epiz = _Ceulerc(z);
    dcomplex eniz = _Ceulerc(cminus(z));
    dcomplex up = _Caddcc(epiz, eniz);
    return cbuild(up.re * .5, up.im * .5);
}

dcomplex csin(dcomplex z) {
    if (z.re == 0) return cbuild(0, sinh(z.im));
    if (z.im == 0) return cbuild(sin(z.re), 0);
    dcomplex epiz = _Ceulerc(z);
    dcomplex eniz = _Ceulerc(cminus(z));
    dcomplex up = _Csubcc(epiz, eniz);
    return cbuild(up.im * .5, -up.re * .5);
}

dcomplex ctan(dcomplex z) {
    if (z.re == 0) return cbuild(0, tanh(z.im));
    if (z.im == 0) return cbuild(tan(z.re), 0);
    dcomplex epiz = _Ceulerc(z);
    dcomplex eniz = _Ceulerc(cminus(z));
    dcomplex up = _Csubcc(epiz, eniz);
    dcomplex down = _Caddcc(epiz, eniz);
    dcomplex niup = cbuild(up.im, -up.re);
    return _Cdivcc(niup, down);
}

dcomplex cacosh(dcomplex z) {
    dcomplex tmp = _Csubcr(_Cpowcr(z, 2.), 1.); // square(z) - 1
    double sqrt_abs_tmp = sqrt(cabs(tmp));
    double theta = carg(tmp) / 2.;
    dcomplex e_iphi = cbuild(cos(theta), sin(theta));
    return clog(_Caddcc(z, _Cmulrc(sqrt_abs_tmp, e_iphi)));
}

dcomplex casinh(dcomplex z) {
    dcomplex tmp = _Caddrc(1., _Cpowcr(z, 2.)); // 1 + square(z)
    double sqrt_abs_tmp = sqrt(cabs(tmp));
    double theta = carg(tmp) / 2.;
    dcomplex e_iphi = cbuild(cos(theta), sin(theta));
    return clog(_Caddcc(z, _Cmulrc(sqrt_abs_tmp, e_iphi)));
}

dcomplex catanh(dcomplex z) {
    dcomplex up = cbuild(1. + z.re, z.im);
    dcomplex down = cbuild(1. - z.re, z.im);
    dcomplex form = clog(_Cdivcc(up, down));
    return _Cmulrc(0.5, form);
}

dcomplex ccosh(dcomplex z) {
    if (z.re == 0) return cbuild(0, cos(z.im));
    if (z.im == 0) return cbuild(cosh(z.re), 0);
    double cosb = cos(z.im);
    double sinb = sin(z.im);
    double prevconst = exp(z.re) / 2.;
    double postconst = exp(-z.re) / 2.;
    dcomplex prev = _Cmulcr(cbuild(cosb, sinb), prevconst);
    dcomplex post = _Cmulcr(cbuild(cosb, -sinb), postconst);
    return _Caddcc(prev, post);
}

dcomplex csinh(dcomplex z) {
    if (z.re == 0) return cbuild(0, sin(z.im));
    if (z.im == 0) return cbuild(sinh(z.re), 0);
    double cosb = cos(z.im);
    double sinb = sin(z.im);
    double prevconst = exp(z.re) / 2.;
    double postconst = exp(-z.re) / 2.;
    dcomplex prev = _Cmulcr(cbuild(cosb, sinb), prevconst);
    dcomplex post = _Cmulcr(cbuild(cosb, -sinb), postconst);
    return _Csubcc(prev, post);
}

dcomplex ctanh(dcomplex z) {
    if (z.re == 0) return cbuild(0, tan(z.im));
    if (z.im == 0) return cbuild(tanh(z.re), 0);
    double cosb = cos(z.im);
    double sinb = sin(z.im);
    double prevconst = exp(z.re);
    double postconst = exp(-z.re);
    dcomplex prev = _Cmulcr(cbuild(cosb, sinb), prevconst);
    dcomplex post = _Cmulcr(cbuild(cosb, -sinb), postconst);
    dcomplex up = _Csubcc(prev, post);
    dcomplex down = _Caddcc(prev, post);
    return _Cdivcc(up, down);
}

dcomplex cexp(dcomplex z) {
    return _Cmulrc(exp(z.re), cbuild(cos(z.im), sin(z.im)));
}

dcomplex clog(dcomplex z) {
    return cbuild(log(norm(z)) / 2, carg(z));
}

dcomplex clog10(dcomplex z) {
    return _Cdivcr(clog(z), log(10.));
}

dcomplex csqrt(dcomplex z) {
    double phi = carg(z) / 2.;
    return _Cmulrc(sqrt(cabs(z)), cbuild(cos(phi), sin(phi)));
}

dcomplex cproj(dcomplex z) {
    if (isinf(z.re) || isinf(z.im))
        return cbuild(INFINITY, 0.);
    if (isnan(z.im))
        return cbuild(NAN, NAN); /* undifined */
    return z;
}

dcomplex csgn(dcomplex z) {
    if (z.im == 0) {
        // sgn function
        double sgn = z.re < -DBL_EPSILON ? -1 : z.re > DBL_EPSILON;
        return cbuild(sgn, 0.);
    }
    return _Cdivcr(z, cabs(z));
}

dcomplex _Ceulerr(double x) {
    return cbuild(cos(x), sin(x));
}

dcomplex _Ceulerc(dcomplex z) {
    return _Cmulcr(_Ceulerr(z.re), exp(-z.im));
}

void _Ceulererr(void) {
    fprintf(stderr, "TypeError: ceuler( <dcomplex, double> z ).\n");
    exit(EXIT_FAILURE);
}

dcomplex _Clognc(dcomplex z, dcomplex n) {
    return _Cdivcc(clog(z), clog(n));
}

dcomplex _Clognr(dcomplex z, double n) {
    return _Cdivcr(clog(z), log(n));
}

dcomplex _Clognz(double z, dcomplex n) {
    return _Cdivrc(log(z), clog(n));
}

void _Clognerr(void) {
    fprintf(stderr, "TypeError: clogn( <dcomplex, double> z, <dcomplex, double> n).\n");
    exit(EXIT_FAILURE);
}

dcomplex _Cpowcc(dcomplex z, dcomplex n) {
    double a = n.re;
    double b = n.im;
    double r = cabs(z);
    double phi = carg(z);
    double theta1 = a * phi;
    double theta2 = b * log(r);
    dcomplex z1 = cbuild(cos(theta1), sin(theta1));
    dcomplex z2 = cbuild(cos(theta2), sin(theta2));
    return _Cmulrc(pow(r, a) * exp(-b * phi), _Cmulcc(z1, z2));
}

dcomplex _Cpowcr(dcomplex z, double n) {
    double theta = n * carg(z);
    return _Cmulrc(pow(cabs(z), n), cbuild(cos(theta), sin(theta)));
}

dcomplex _Cpowrc(double z, dcomplex n) {
    double tmp = pow(z, n.re);
    double theta = n.im * log(z > 0 ? z : -z);
    dcomplex z0 = cbuild(cos(theta), sin(theta));
    if (z > 0) return _Cmulrc(tmp, z0);
    if (z < 0) return _Cmulrc(tmp * exp(-PI * n.im), z0);
}

void _Cpowerr(void) {
    fprintf(stderr, "TypeError: cpow( <dcomplex, double> z, <dcomplex, double> n).\n");
    exit(EXIT_FAILURE);
}

dcomplex cbuild(double re, double im) {
    dcomplex z = { re, im };
    return z;
}

dcomplex _Caddcc(dcomplex a, dcomplex b) {
    return cbuild(a.re + b.re, a.im + b.im);
}

dcomplex _Caddcr(dcomplex a, double b) {
    return cbuild(a.re + b, a.im);
}

dcomplex _Caddrc(double a, dcomplex b) {
    return cbuild(a + b.re, b.im);
}

void _Cadderr(void) {
    fprintf(stderr, "TypeError: cadd( <dcomplex, double> a, <dcomplex, double> b).\n");
    exit(EXIT_FAILURE);
}

dcomplex _Csubcc(dcomplex a, dcomplex b) {
    return cbuild(a.re - b.re, a.im - b.im);
}

dcomplex _Csubcr(dcomplex a, double b) {
    return cbuild(a.re - b, a.im);
}

dcomplex _Csubrc(double a, dcomplex b) {
    return cbuild(a - b.re, -b.im);
}

void _Csuberr(void) {
    fprintf(stderr, "TypeError: csub( <dcomplex, double> a, <dcomplex, double> b).\n");
    exit(EXIT_FAILURE);
}

dcomplex _Cmulcc(dcomplex a, dcomplex b) {
    return cbuild(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

dcomplex _Cmulcr(dcomplex a, double b) {
    return cbuild(a.re * b, a.im * b);
}

dcomplex _Cmulrc(double a, dcomplex b) {
    return cbuild(a * b.re, a * b.im);
}

void _Cmulerr(void) {
    fprintf(stderr, "TypeError: cmul( <dcomplex, double> a, <dcomplex, double> b).\n");
    exit(EXIT_FAILURE);
}

dcomplex _Cdivcc(dcomplex a, dcomplex b) {
    double arg = carg(a) - carg(b);
    double abs = cabs(a) / cabs(b);
    return cbuild(abs * cos(arg), abs * sin(arg));
}

dcomplex _Cdivcr(dcomplex a, double b) {
    return cbuild(a.re / b, a.im / b);
}

dcomplex _Cdivrc(double a, dcomplex b) {
    return _Cmulcr(conj(b), a / norm(b));
}

void _Cdiverr(void) {
    fprintf(stderr, "TypeError: cdiv( <dcomplex, double> a, <dcomplex, double> b).\n");
    exit(EXIT_FAILURE);
}


/* float complex type ------------------------------ */

int cisequalf(fcomplex a, fcomplex b) {
    if (a.re ==  b.re && a.im ==  b.im) return  1;
    if (a.re == -b.re && a.im == -b.im) return -1;
    return 0;
}

float crealf(fcomplex z) {
    return z.re;
}

float cimagf(fcomplex z) {
    return z.im;
}

float cargf(fcomplex z) {
    return atan2f(z.im, z.re);
}

float cabsf(fcomplex z) {
    return sqrtf(z.re * z.re + z.im * z.im);
}

float normf(fcomplex z) {
    return z.re * z.re + z.im * z.im;
}

fcomplex cminusf(fcomplex z) {
    return cbuildf(-z.re, -z.im);
}

fcomplex conjf(fcomplex z) {
    fcomplex r = {z.re, -z.im};
    return r;
}

fcomplex cacosf(fcomplex z) {
    fcomplex tmp = _CFsubrc(1., _CFpowcr(z, 2.)); // 1 - square(z)
    fcomplex isqrt_abs_tmp = _CFmulcr(cbuildf(0., 1.), sqrtf(cabsf(tmp)));
    float theta = cargf(tmp) / 2.;
    fcomplex e_iphi = cbuildf(cosf(theta), sinf(theta));
    fcomplex ln_longform = clogf(_CFaddcc(z, _CFmulcc(isqrt_abs_tmp, e_iphi)));
    return _CFmulcc(cbuildf(0., -1.), ln_longform);
}

fcomplex casinf(fcomplex z) {
    fcomplex iz = cbuildf(-z.im, z.re);
    fcomplex tmp = _CFsubrc(1., _CFpowcr(z, 2.)); // 1 - square(z)
    float sqrt_abs_tmp = sqrtf(cabsf(tmp));
    float theta = cargf(tmp) / 2.;
    fcomplex e_iphi = cbuildf(cosf(theta), sinf(theta));
    fcomplex ln_longform = clogf(_CFaddcc(iz, _CFmulrc(sqrt_abs_tmp, e_iphi)));
    return _CFmulcc(cbuildf(0., -1.), ln_longform);
}

fcomplex catanf(fcomplex z) {
    fcomplex up = cbuildf(-z.re, 1. - z.im);
    fcomplex down = cbuildf(z.re, 1. + z.im);
    fcomplex post = clogf(_CFdivcc(up, down));
    fcomplex ans = _CFmulcc(cbuildf(0., -0.5), post);
    if (z.re == 0) {
        if (z.im == 1) return cbuildf(0, INFINITY);
        if (z.im == -1) return cbuildf(0, -INFINITY);
        if (z.im > 0) return _CFaddcr(ans, PI);
        return ans;
    }
}

fcomplex ccosf(fcomplex z) {
    if (z.re == 0) return cbuildf(coshf(z.im), 0);
    if (z.im == 0) return cbuildf(cosf(z.re), 0);
    fcomplex epiz = _CFeulerc(z);
    fcomplex eniz = _CFeulerc(cminusf(z));
    fcomplex up = _CFaddcc(epiz, eniz);
    return cbuildf(up.re * .5, up.im * .5);
}

fcomplex csinf(fcomplex z) {
    if (z.re == 0) return cbuildf(0, sinhf(z.im));
    if (z.im == 0) return cbuildf(sinf(z.re), 0);
    fcomplex epiz = _CFeulerc(z);
    fcomplex eniz = _CFeulerc(cminusf(z));
    fcomplex up = _CFsubcc(epiz, eniz);
    return cbuildf(up.im * .5, -up.re * .5);
}

fcomplex ctanf(fcomplex z) {
    if (z.re == 0) return cbuildf(0, tanhf(z.im));
    if (z.im == 0) return cbuildf(tanf(z.re), 0);
    fcomplex epiz = _CFeulerc(z);
    fcomplex eniz = _CFeulerc(cminusf(z));
    fcomplex up = _CFsubcc(epiz, eniz);
    fcomplex down = _CFaddcc(epiz, eniz);
    fcomplex niup = cbuildf(up.im, -up.re);
    return _CFdivcc(niup, down);
}

fcomplex cacoshf(fcomplex z) {
    fcomplex tmp = _CFsubcr(_CFpowcr(z, 2.), 1.); // square(z) - 1
    float sqrt_abs_tmp = sqrtf(cabsf(tmp));
    float theta = cargf(tmp) / 2.;
    fcomplex e_iphi = cbuildf(cosf(theta), sinf(theta));
    return clogf(_CFaddcc(z, _CFmulrc(sqrt_abs_tmp, e_iphi)));
}

fcomplex casinhf(fcomplex z) {
    fcomplex tmp = _CFaddrc(1., _CFpowcr(z, 2.)); // 1 + square(z)
    float sqrt_abs_tmp = sqrtf(cabsf(tmp));
    float theta = cargf(tmp) / 2.;
    fcomplex e_iphi = cbuildf(cosf(theta), sinf(theta));
    return clogf(_CFaddcc(z, _CFmulrc(sqrt_abs_tmp, e_iphi)));
}

fcomplex catanhf(fcomplex z) {
    fcomplex up = cbuildf(1. + z.re, z.im);
    fcomplex down = cbuildf(1. - z.re, z.im);
    fcomplex form = clogf(_CFdivcc(up, down));
    return _CFmulrc(0.5, form);
}

fcomplex ccoshf(fcomplex z) {
    if (z.re == 0) return cbuildf(0, cosf(z.im));
    if (z.im == 0) return cbuildf(coshf(z.re), 0);
    float cosb = cosf(z.im);
    float sinb = sinf(z.im);
    float prevconst = expf(z.re) / 2.;
    float postconst = expf(-z.re) / 2.;
    fcomplex prev = _CFmulcr(cbuildf(cosb, sinb), prevconst);
    fcomplex post = _CFmulcr(cbuildf(cosb, -sinb), postconst);
    return _CFaddcc(prev, post);
}

fcomplex csinhf(fcomplex z) {
    if (z.re == 0) return cbuildf(0, sinf(z.im));
    if (z.im == 0) return cbuildf(sinhf(z.re), 0);
    float cosb = cosf(z.im);
    float sinb = sinf(z.im);
    float prevconst = expf(z.re) / 2.;
    float postconst = expf(-z.re) / 2.;
    fcomplex prev = _CFmulcr(cbuildf(cosb, sinb), prevconst);
    fcomplex post = _CFmulcr(cbuildf(cosb, -sinb), postconst);
    return _CFsubcc(prev, post);
}

fcomplex ctanhf(fcomplex z) {
    if (z.re == 0) return cbuildf(0, tanf(z.im));
    if (z.im == 0) return cbuildf(tanhf(z.re), 0);
    float cosb = cosf(z.im);
    float sinb = sinf(z.im);
    float prevconst = expf(z.re);
    float postconst = expf(-z.re);
    fcomplex prev = _CFmulcr(cbuildf(cosb, sinb), prevconst);
    fcomplex post = _CFmulcr(cbuildf(cosb, -sinb), postconst);
    fcomplex up = _CFsubcc(prev, post);
    fcomplex down = _CFaddcc(prev, post);
    return _CFdivcc(up, down);
}

fcomplex cexpf(fcomplex z) {
    return _CFmulrc(expf(z.re), cbuildf(cosf(z.im), sinf(z.im)));
}

fcomplex clogf(fcomplex z) {
    return cbuildf(logf(normf(z)) / 2, cargf(z));
}

fcomplex clog10f(fcomplex z) {
    return _CFdivcr(clogf(z), logf(10.));
}

fcomplex csqrtf(fcomplex z) {
    float phi = cargf(z) / 2.;
    return _CFmulrc(sqrtf(cabsf(z)), cbuildf(cosf(phi), sinf(phi)));
}

fcomplex cprojf(fcomplex z) {
    if (isinf(z.re) || isinf(z.im))
        return cbuildf(INFINITY, 0.);
    if (isnan(z.im))
        return cbuildf(NAN, NAN); /* undifined */
    return z;
}

fcomplex csgnf(fcomplex z) {
    if (z.im == 0) {
        // sgn function
        float sgn = z.re < -FLT_EPSILON ? -1 : z.re > FLT_EPSILON;
        return cbuildf(sgn, 0.);
    }
    return _CFdivcr(z, cabsf(z));
}

fcomplex _CFeulerr(float x) {
    return cbuildf(cosf(x), sinf(x));
}

fcomplex _CFeulerc(fcomplex z) {
    return _CFmulcr(_CFeulerr(z.re), expf(-z.im));
}

void _CFeulererr(void) {
    fprintf(stderr, "TypeError: ceulerf( <fcomplex, float> z ).\n");
    exit(EXIT_FAILURE);
}

fcomplex _CFlognc(fcomplex z, fcomplex n) {
    return _CFdivcc(clogf(z), clogf(n));
}

fcomplex _CFlognr(fcomplex z, float n) {
    return _CFdivcr(clogf(z), logf(n));
}

fcomplex _CFlognz(float z, fcomplex n) {
    return _CFdivrc(logf(z), clogf(n));
}

void _CFlognerr(void) {
    fprintf(stderr, "TypeError: clognf( <fcomplex, float> z, <fcomplex, float> n).\n");
    exit(EXIT_FAILURE);
}

fcomplex _CFpowcc(fcomplex z, fcomplex n) {
    float a = n.re;
    float b = n.im;
    float r = cabsf(z);
    float phi = cargf(z);
    float theta1 = a * phi;
    float theta2 = b * logf(r);
    fcomplex z1 = cbuildf(cosf(theta1), sinf(theta1));
    fcomplex z2 = cbuildf(cosf(theta2), sinf(theta2));
    return _CFmulrc(powf(r, a) * expf(-b * phi), _CFmulcc(z1, z2));
}

fcomplex _CFpowcr(fcomplex z, float n) {
    float theta = n * cargf(z);
    return _CFmulrc(powf(cabsf(z), n), cbuildf(cosf(theta), sinf(theta)));
}

fcomplex _CFpowrc(float z, fcomplex n) {
    float tmp = powf(z, n.re);
    float theta = n.im * logf(z > 0 ? z : -z);
    fcomplex z0 = cbuildf(cosf(theta), sinf(theta));
    if (z > 0) return _CFmulrc(tmp, z0);
    if (z < 0) return _CFmulrc(tmp * expf(-PI * n.im), z0);
}

void _CFpowerr(void) {
    fprintf(stderr, "TypeError: cpowf( <fcomplex, float> z, <fcomplex, float> n).\n");
    exit(EXIT_FAILURE);
}

fcomplex cbuildf(float re, float im) {
    fcomplex z = { re, im };
    return z;
}

fcomplex _CFaddcc(fcomplex a, fcomplex b) {
    return cbuildf(a.re + b.re, a.im + b.im);
}

fcomplex _CFaddcr(fcomplex a, float b) {
    return cbuildf(a.re + b, a.im);
}

fcomplex _CFaddrc(float a, fcomplex b) {
    return cbuildf(a + b.re, b.im);
}

void _CFadderr(void) {
    fprintf(stderr, "TypeError: caddf( <fcomplex, float> a, <fcomplex, float> b).\n");
    exit(EXIT_FAILURE);
}

fcomplex _CFsubcc(fcomplex a, fcomplex b) {
    return cbuildf(a.re - b.re, a.im - b.im);
}

fcomplex _CFsubcr(fcomplex a, float b) {
    return cbuildf(a.re - b, a.im);
}

fcomplex _CFsubrc(float a, fcomplex b) {
    return cbuildf(a - b.re, -b.im);
}

void _CFsuberr(void) {
    fprintf(stderr, "TypeError: csubf( <fcomplex, float> a, <fcomplex, float> b).\n");
    exit(EXIT_FAILURE);
}

fcomplex _CFmulcc(fcomplex a, fcomplex b) {
    return cbuildf(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

fcomplex _CFmulcr(fcomplex a, float b) {
    return cbuildf(a.re * b, a.im * b);
}

fcomplex _CFmulrc(float a, fcomplex b) {
    return cbuildf(a * b.re, a * b.im);
}

void _CFmulerr(void) {
    fprintf(stderr, "TypeError: cmulf( <fcomplex, float> a, <fcomplex, float> b).\n");
    exit(EXIT_FAILURE);
}

fcomplex _CFdivcc(fcomplex a, fcomplex b) {
    float arg = cargf(a) - cargf(b);
    float abs = cabsf(a) / cabsf(b);
    return cbuildf(abs * cosf(arg), abs * sinf(arg));
}

fcomplex _CFdivcr(fcomplex a, float b) {
    return cbuildf(a.re / b, a.im / b);
}

fcomplex _CFdivrc(float a, fcomplex b) {
    return _CFmulcr(conjf(b), a / normf(b));
}

void _CFdiverr(void) {
    fprintf(stderr, "TypeError: cdivf( <fcomplex, float> a, <fcomplex, float> b).\n");
    exit(EXIT_FAILURE);
}


/* long double complex type ------------------------------ */

int cisequall(lcomplex a, lcomplex b) {
    if (a.re ==  b.re && a.im ==  b.im) return  1;
    if (a.re == -b.re && a.im == -b.im) return -1;
    return 0;
}

long double creall(lcomplex z) {
    return z.re;
}

long double cimagl(lcomplex z) {
    return z.im;
}

long double cargl(lcomplex z) {
    return atan2(z.im, z.re);
}

long double cabsl(lcomplex z) {
    return sqrt(z.re * z.re + z.im * z.im);
}

long double norml(lcomplex z) {
    return z.re * z.re + z.im * z.im;
}

lcomplex cminusl(lcomplex z) {
    return cbuildl(-z.re, -z.im);
}

lcomplex conjl(lcomplex z) {
    lcomplex r = {z.re, -z.im};
    return r;
}

lcomplex cacosl(lcomplex z) {
    lcomplex tmp = _CLsubrc(1., _CLpowcr(z, 2.)); // 1 - square(z)
    lcomplex isqrt_abs_tmp = _CLmulcr(cbuildl(0., 1.), sqrt(cabsl(tmp)));
    long double theta = cargl(tmp) / 2.;
    lcomplex e_iphi = cbuildl(cos(theta), sin(theta));
    lcomplex ln_longform = clogl(_CLaddcc(z, _CLmulcc(isqrt_abs_tmp, e_iphi)));
    return _CLmulcc(cbuildl(0., -1.), ln_longform);
}

lcomplex casinl(lcomplex z) {
    lcomplex iz = cbuildl(-z.im, z.re);
    lcomplex tmp = _CLsubrc(1., _CLpowcr(z, 2.)); // 1 - square(z)
    long double sqrt_abs_tmp = sqrt(cabsl(tmp));
    long double theta = cargl(tmp) / 2.;
    lcomplex e_iphi = cbuildl(cos(theta), sin(theta));
    lcomplex ln_longform = clogl(_CLaddcc(iz, _CLmulrc(sqrt_abs_tmp, e_iphi)));
    return _CLmulcc(cbuildl(0., -1.), ln_longform);
}

lcomplex catanl(lcomplex z) {
    lcomplex up = cbuildl(-z.re, 1. - z.im);
    lcomplex down = cbuildl(z.re, 1. + z.im);
    lcomplex post = clogl(_CLdivcc(up, down));
    lcomplex ans = _CLmulcc(cbuildl(0., -0.5), post);
    if (z.re == 0) {
        if (z.im == 1) return cbuildl(0, INFINITY);
        if (z.im == -1) return cbuildl(0, -INFINITY);
        if (z.im > 0) return _CLaddcr(ans, PI);
        return ans;
    }
}

lcomplex ccosl(lcomplex z) {
    if (z.re == 0) return cbuildl(cosh(z.im), 0);
    if (z.im == 0) return cbuildl(cos(z.re), 0);
    lcomplex epiz = _CLeulerc(z);
    lcomplex eniz = _CLeulerc(cminusl(z));
    lcomplex up = _CLaddcc(epiz, eniz);
    return cbuildl(up.re * .5, up.im * .5);
}

lcomplex csinl(lcomplex z) {
    if (z.re == 0) return cbuildl(0, sinh(z.im));
    if (z.im == 0) return cbuildl(sin(z.re), 0);
    lcomplex epiz = _CLeulerc(z);
    lcomplex eniz = _CLeulerc(cminusl(z));
    lcomplex up = _CLsubcc(epiz, eniz);
    return cbuildl(up.im * .5, -up.re * .5);
}

lcomplex ctanl(lcomplex z) {
    if (z.re == 0) return cbuildl(0, tanh(z.im));
    if (z.im == 0) return cbuildl(tan(z.re), 0);
    lcomplex epiz = _CLeulerc(z);
    lcomplex eniz = _CLeulerc(cminusl(z));
    lcomplex up = _CLsubcc(epiz, eniz);
    lcomplex down = _CLaddcc(epiz, eniz);
    lcomplex niup = cbuildl(up.im, -up.re);
    return _CLdivcc(niup, down);
}

lcomplex cacoshl(lcomplex z) {
    lcomplex tmp = _CLsubcr(_CLpowcr(z, 2.), 1.); // square(z) - 1
    long double sqrt_abs_tmp = sqrt(cabsl(tmp));
    long double theta = cargl(tmp) / 2.;
    lcomplex e_iphi = cbuildl(cos(theta), sin(theta));
    return clogl(_CLaddcc(z, _CLmulrc(sqrt_abs_tmp, e_iphi)));
}

lcomplex casinhl(lcomplex z) {
    lcomplex tmp = _CLaddrc(1., _CLpowcr(z, 2.)); // 1 + square(z)
    long double sqrt_abs_tmp = sqrt(cabsl(tmp));
    long double theta = cargl(tmp) / 2.;
    lcomplex e_iphi = cbuildl(cos(theta), sin(theta));
    return clogl(_CLaddcc(z, _CLmulrc(sqrt_abs_tmp, e_iphi)));
}

lcomplex catanhl(lcomplex z) {
    lcomplex up = cbuildl(1. + z.re, z.im);
    lcomplex down = cbuildl(1. - z.re, z.im);
    lcomplex form = clogl(_CLdivcc(up, down));
    return _CLmulrc(0.5, form);
}

lcomplex ccoshl(lcomplex z) {
    if (z.re == 0) return cbuildl(0, cos(z.im));
    if (z.im == 0) return cbuildl(cosh(z.re), 0);
    long double cosb = cos(z.im);
    long double sinb = sin(z.im);
    long double prevconst = exp(z.re) / 2.;
    long double postconst = exp(-z.re) / 2.;
    lcomplex prev = _CLmulcr(cbuildl(cosb, sinb), prevconst);
    lcomplex post = _CLmulcr(cbuildl(cosb, -sinb), postconst);
    return _CLaddcc(prev, post);
}

lcomplex csinhl(lcomplex z) {
    if (z.re == 0) return cbuildl(0, sin(z.im));
    if (z.im == 0) return cbuildl(sinh(z.re), 0);
    long double cosb = cos(z.im);
    long double sinb = sin(z.im);
    long double prevconst = exp(z.re) / 2.;
    long double postconst = exp(-z.re) / 2.;
    lcomplex prev = _CLmulcr(cbuildl(cosb, sinb), prevconst);
    lcomplex post = _CLmulcr(cbuildl(cosb, -sinb), postconst);
    return _CLsubcc(prev, post);
}

lcomplex ctanhl(lcomplex z) {
    if (z.re == 0) return cbuildl(0, tan(z.im));
    if (z.im == 0) return cbuildl(tanh(z.re), 0);
    long double cosb = cos(z.im);
    long double sinb = sin(z.im);
    long double prevconst = exp(z.re);
    long double postconst = exp(-z.re);
    lcomplex prev = _CLmulcr(cbuildl(cosb, sinb), prevconst);
    lcomplex post = _CLmulcr(cbuildl(cosb, -sinb), postconst);
    lcomplex up = _CLsubcc(prev, post);
    lcomplex down = _CLaddcc(prev, post);
    return _CLdivcc(up, down);
}

lcomplex cexpl(lcomplex z) {
    return _CLmulrc(exp(z.re), cbuildl(cos(z.im), sin(z.im)));
}

lcomplex clogL(lcomplex z) {
    return cbuildl(log(norml(z)) / 2, cargl(z));
}

lcomplex clog10l(lcomplex z) {
    return _CLdivcr(clogl(z), log(10.));
}

lcomplex csqrtl(lcomplex z) {
    long double phi = cargl(z) / 2.;
    return _CLmulrc(sqrt(cabsl(z)), cbuildl(cos(phi), sin(phi)));
}

lcomplex cprojl(lcomplex z) {
    if (isinf(z.re) || isinf(z.im))
        return cbuildl(INFINITY, 0.);
    if (isnan(z.im))
        return cbuildl(NAN, NAN); /* undifined */
    return z;
}

lcomplex csgnl(lcomplex z) {
    if (z.im == 0) {
        // sgn function
        long double sgn = z.re < -LDBL_EPSILON ? -1 : z.re > LDBL_EPSILON;
        return cbuildl(sgn, 0.);
    }
    return _CLdivcr(z, cabsl(z));
}

lcomplex _CLeulerr(long double x) {
    return cbuildl(cos(x), sin(x));
}

lcomplex _CLeulerc(lcomplex z) {
    return _CLmulcr(_CLeulerr(z.re), exp(-z.im));
}

void _CLeulererr(void) {
    fprintf(stderr, "TypeError: ceulerl( <lcomplex, long double> z ).\n");
    exit(EXIT_FAILURE);
}

lcomplex _CLlognc(lcomplex z, lcomplex n) {
    return _CLdivcc(clogl(z), clogl(n));
}

lcomplex _CLlognr(lcomplex z, long double n) {
    return _CLdivcr(clogl(z), log(n));
}

lcomplex _CLlognz(long double z, lcomplex n) {
    return _CLdivrc(log(z), clogl(n));
}

void _CLlognerr(void) {
    fprintf(stderr, "TypeError: clognl( <lcomplex, long double> z, <lcomplex, long double> n).\n");
    exit(EXIT_FAILURE);
}

lcomplex _CLpowcc(lcomplex z, lcomplex n) {
    long double a = n.re;
    long double b = n.im;
    long double r = cabsl(z);
    long double phi = cargl(z);
    long double theta1 = a * phi;
    long double theta2 = b * log(r);
    lcomplex z1 = cbuildl(cos(theta1), sin(theta1));
    lcomplex z2 = cbuildl(cos(theta2), sin(theta2));
    return _CLmulrc(pow(r, a) * exp(-b * phi), _CLmulcc(z1, z2));
}

lcomplex _CLpowcr(lcomplex z, long double n) {
    long double theta = n * cargl(z);
    return _CLmulrc(pow(cabsl(z), n), cbuildl(cos(theta), sin(theta)));
}

lcomplex _CLpowrc(long double z, lcomplex n) {
    long double tmp = pow(z, n.re);
    long double theta = n.im * log(z > 0 ? z : -z);
    lcomplex z0 = cbuildl(cos(theta), sin(theta));
    if (z > 0) return _CLmulrc(tmp, z0);
    if (z < 0) return _CLmulrc(tmp * exp(-PI * n.im), z0);
}

void _CLpowerr(void) {
    fprintf(stderr, "TypeError: cpowl( <lcomplex, long double> z, <lcomplex, long double> n).\n");
    exit(EXIT_FAILURE);
}

lcomplex cbuildl(long double re, long double im) {
    lcomplex z = { re, im };
    return z;
}

lcomplex _CLaddcc(lcomplex a, lcomplex b) {
    return cbuildl(a.re + b.re, a.im + b.im);
}

lcomplex _CLaddcr(lcomplex a, long double b) {
    return cbuildl(a.re + b, a.im);
}

lcomplex _CLaddrc(long double a, lcomplex b) {
    return cbuildl(a + b.re, b.im);
}

void _CLadderr(void) {
    fprintf(stderr, "TypeError: caddl( <lcomplex, long double> a, <lcomplex, long double> b).\n");
    exit(EXIT_FAILURE);
}

lcomplex _CLsubcc(lcomplex a, lcomplex b) {
    return cbuildl(a.re - b.re, a.im - b.im);
}

lcomplex _CLsubcr(lcomplex a, long double b) {
    return cbuildl(a.re - b, a.im);
}

lcomplex _CLsubrc(long double a, lcomplex b) {
    return cbuildl(a - b.re, -b.im);
}

void _CLsuberr(void) {
    fprintf(stderr, "TypeError: csubl( <lcomplex, long double> a, <lcomplex, long double> b).\n");
    exit(EXIT_FAILURE);
}

lcomplex _CLmulcc(lcomplex a, lcomplex b) {
    return cbuildl(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

lcomplex _CLmulcr(lcomplex a, long double b) {
    return cbuildl(a.re * b, a.im * b);
}

lcomplex _CLmulrc(long double a, lcomplex b) {
    return cbuildl(a * b.re, a * b.im);
}

void _CLmulerr(void) {
    fprintf(stderr, "TypeError: cmull( <lcomplex, long double> a, <lcomplex, long double> b).\n");
    exit(EXIT_FAILURE);
}

lcomplex _CLdivcc(lcomplex a, lcomplex b) {
    long double arg = cargl(a) - cargl(b);
    long double abs = cabsl(a) / cabsl(b);
    return cbuildl(abs * cos(arg), abs * sin(arg));
}

lcomplex _CLdivcr(lcomplex a, long double b) {
    return cbuildl(a.re / b, a.im / b);
}

lcomplex _CLdivrc(long double a, lcomplex b) {
    return _CLmulcr(conjl(b), a / norml(b));
}

void _CLdiverr(void) {
    fprintf(stderr, "TypeError: cdivl( <lcomplex, long double> a, <lcomplex, long double> b).\n");
    exit(EXIT_FAILURE);
}

#endif
