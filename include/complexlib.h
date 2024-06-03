//
// conmplexlib.h
//
//      header of complex library
//

#pragma once

#if defined(ARDUINO)

#include <math.h>

typedef struct double_complex_t {
    double re;
    double im;
} dcomplex;

typedef struct float_complex_t {
    float re;
    float im;
} fcomplex;

typedef struct long_double_complex_t {
    long double re;
    long double im;
} lcomplex;

/* double complex type ------------------------------ */

int cisequal(dcomplex a, dcomplex b);

double creal(dcomplex z);
double cimag(dcomplex z);
double carg(dcomplex z);
double cabs(dcomplex z);
double norm(dcomplex z);

dcomplex cminus(dcomplex z);
dcomplex conj(dcomplex z);
dcomplex cacos(dcomplex z);
dcomplex casin(dcomplex z);
dcomplex catan(dcomplex z);
dcomplex ccos(dcomplex z);
dcomplex csin(dcomplex z);
dcomplex ctan(dcomplex z);
dcomplex cacosh(dcomplex z);
dcomplex casinh(dcomplex z);
dcomplex catanh(dcomplex z);
dcomplex ccosh(dcomplex z);
dcomplex csinh(dcomplex z);
dcomplex ctanh(dcomplex z);
dcomplex cexp(dcomplex z);
dcomplex clog(dcomplex z);
dcomplex clog10(dcomplex z);
dcomplex csqrt(dcomplex z);
dcomplex cproj(dcomplex z);
dcomplex csgn(dcomplex z);

dcomplex _Ceulerr(double x);
dcomplex _Ceulerc(dcomplex z);
void _Ceulererr(void);

#define ceuler(z) _Generic((a),     \
    dcomplex:_Ceulerc,              \
    double: _Ceulerr,               \
    default: _Ceulererr             \
)(z)

dcomplex _Clognc(dcomplex z, dcomplex n);
dcomplex _Clognr(dcomplex z, double n);
dcomplex _Clognz(double z, dcomplex n);
void _Clognerr(void);

#define clogn(a, b) _Generic((a),   \
    dcomplex: _Generic((b),         \
        dcomplex: _Clognc,          \
        double: _Clognr,            \
        default: _Clognerr          \
    ),                              \
    double: _Generic((b),           \
        dcomplex: _Clognz,          \
        default: _Clognerr          \
    ),                              \
    default: _Clognerr              \
)(a, b)

dcomplex _Cpowcc(dcomplex z, dcomplex n);
dcomplex _Cpowcr(dcomplex z, double n);
dcomplex _Cpowrc(double z, dcomplex n);
void _Cpowerr(void);

#define cpow(a, b) _Generic((a),    \
    dcomplex: _Generic((b),         \
        dcomplex: _Cpowcc,          \
        double: _Cpowcr,            \
        default: _Cpowerr           \
    ),                              \
    double: _Generic((b),           \
        dcomplex: _Cpowrc,          \
        default: _Cpowerr           \
    ),                              \
    default: _Cpowerr               \
)(a, b)

dcomplex cbuild(double real, double imag);

dcomplex _Caddcc(dcomplex a, dcomplex b);
dcomplex _Caddcr(dcomplex a, double b);
dcomplex _Caddrc(double a, dcomplex b);
void _Cadderr(void);

#define cadd(a, b) _Generic((a),    \
    dcomplex: _Generic((b),         \
        dcomplex: _Caddcc,          \
        double: _Caddcr,            \
        default: _Cadderr           \
    ),                              \
    double: _Generic((b),           \
        dcomplex: _Caddrc,          \
        default: _Cadderr           \
    ),                              \
    default: _Cadderr               \
)(a, b)

dcomplex _Csubcc(dcomplex a, dcomplex b);
dcomplex _Csubcr(dcomplex a, double b);
dcomplex _Csubrc(double a, dcomplex b);
void _Csuberr(void);

#define csub(a, b) _Generic((a),    \
    dcomplex: _Generic((b),         \
        dcomplex: _Csubcc,          \
        double: _Csubcr,            \
        default: _Csuberr           \
    ),                              \
    double: _Generic((b),           \
        dcomplex: _Csubcr,          \
        default: _Csuberr           \
    ),                              \
    default: _Csuberr               \
)(a, b)

dcomplex _Cmulcc(dcomplex a, dcomplex b);
dcomplex _Cmulcr(dcomplex a, double b);
dcomplex _Cmulrc(double a, dcomplex b);
void _Cmulerr(void);

#define cmul(a, b) _Generic((a),    \
    dcomplex: _Generic((b),         \
        dcomplex: _Cmulcc,          \
        double  : _Cmulcr,          \
        default : _Cmulerr          \
    ),                              \
    double: _Generic((b),           \
        dcomplex: _Cmulrc,          \
        default : _Cmulerr          \
    ),                              \
    default: _Cmulerr               \
)(a, b)

dcomplex _Cdivcc(dcomplex a, dcomplex b);
dcomplex _Cdivcr(dcomplex a, double b);
dcomplex _Cdivrc(double a, dcomplex b);
void _Cdiverr(void);

#define cdiv(a, b) _Generic((a),    \
    dcomplex: _Generic((b),         \
        dcomplex: _Cdivcc,          \
        double: _Cdivcr,            \
        default: _Cdiverr           \
    ),                              \
    double: _Generic((b),           \
        dcomplex: _Cdivcr,          \
        default: _Cdiverr           \
    ),                              \
    default: _Cdiverr               \
)(a, b)


/* float complex type ------------------------------ */

int cisequalf(fcomplex a, fcomplex b);

float crealf(fcomplex z);
float cimagf(fcomplex z);
float cargf(fcomplex z);
float cabsf(fcomplex z);
float normf(fcomplex z);

fcomplex cminusf(fcomplex z);
fcomplex conjf(fcomplex z);
fcomplex cacosf(fcomplex z);
fcomplex casinf(fcomplex z);
fcomplex catanf(fcomplex z);
fcomplex ccosf(fcomplex z);
fcomplex csinf(fcomplex z);
fcomplex ctanf(fcomplex z);
fcomplex cacoshf(fcomplex z);
fcomplex casinhf(fcomplex z);
fcomplex catanhf(fcomplex z);
fcomplex ccoshf(fcomplex z);
fcomplex csinhf(fcomplex z);
fcomplex ctanhf(fcomplex z);
fcomplex cexpf(fcomplex z);
fcomplex clogf(fcomplex z);
fcomplex clog10f(fcomplex z);
fcomplex csqrtf(fcomplex z);
fcomplex cprojf(fcomplex z);
fcomplex csgnf(fcomplex z);

fcomplex _CFeulerr(float x);
fcomplex _CFeulerc(fcomplex z);
void _CFeulererr(void);

#define ceulerf(z) _Generic((a),     \
    fcomplex:_CFeulerc,              \
    float: _CFeulerr,                \
    default: _CFeulererr             \
)(z)

fcomplex _CFlognc(fcomplex z, fcomplex n);
fcomplex _CFlognr(fcomplex z, float n);
fcomplex _CFlognz(float z, fcomplex n);
void _CFlognerr(void);

#define clognf(a, b) _Generic((a),   \
    fcomplex: _Generic((b),          \
        fcomplex: _CFlognc,          \
        float: _CFlognr,             \
        default: _CFlognerr          \
    ),                               \
    float: _Generic((b),             \
        fcomplex: _CFlognz,          \
        default: _CFlognerr          \
    ),                               \
    default: _CFlognerr              \
)(a, b)

fcomplex _CFpowcc(fcomplex z, fcomplex n);
fcomplex _CFpowcr(fcomplex z, float n);
fcomplex _CFpowrc(float z, fcomplex n);
void _CFpowerr(void);

#define cpowf(a, b) _Generic((a),    \
    fcomplex: _Generic((b),          \
        fcomplex: _CFpowcc,          \
        float: _CFpowcr,             \
        default: _CFpowerr           \
    ),                               \
    float: _Generic((b),             \
        fcomplex: _CFpowrc,          \
        default: _CFpowerr           \
    ),                               \
    default: _CFpowerr               \
)(a, b)

fcomplex cbuildf(float real, float imag);

fcomplex _CFaddcc(fcomplex a, fcomplex b);
fcomplex _CFaddcr(fcomplex a, float b);
fcomplex _CFaddrc(float a, fcomplex b);
void _CFadderr(void);

#define caddf(a, b) _Generic((a),    \
    fcomplex: _Generic((b),          \
        fcomplex: _CFaddcc,          \
        float: _CFaddcr,             \
        default: _CFadderr           \
    ),                               \
    float: _Generic((b),             \
        fcomplex: _CFaddrc,          \
        default: _CFadderr           \
    ),                               \
    default: _CFadderr               \
)(a, b)

fcomplex _CFsubcc(fcomplex a, fcomplex b);
fcomplex _CFsubcr(fcomplex a, float b);
fcomplex _CFsubrc(float a, fcomplex b);
void _CFsuberr(void);

#define csubf(a, b) _Generic((a),    \
    fcomplex: _Generic((b),          \
        fcomplex: _CFsubcc,          \
        float: _CFsubcr,             \
        default: _CFsuberr           \
    ),                               \
    float: _Generic((b),             \
        fcomplex: _CFsubcr,          \
        default: _CFsuberr           \
    ),                               \
    default: _CFsuberr               \
)(a, b)

fcomplex _CFmulcc(fcomplex a, fcomplex b);
fcomplex _CFmulcr(fcomplex a, float b);
fcomplex _CFmulrc(float a, fcomplex b);
void _CFmulerr(void);

#define cmulf(a, b) _Generic((a),    \
    fcomplex: _Generic((b),          \
        fcomplex: _CFmulcc,          \
        float  : _CFmulcr,           \
        default : _CFmulerr          \
    ),                               \
    float: _Generic((b),             \
        fcomplex: _CFmulrc,          \
        default : _CFmulerr          \
    ),                               \
    default: _CFmulerr               \
)(a, b)

fcomplex _CFdivcc(fcomplex a, fcomplex b);
fcomplex _CFdivcr(fcomplex a, float b);
fcomplex _CFdivrc(float a, fcomplex b);
void _CFdiverr(void);

#define cdivf(a, b) _Generic((a),    \
    fcomplex: _Generic((b),          \
        fcomplex: _CFdivcc,          \
        float: _CFdivcr,             \
        default: _CFdiverr           \
    ),                               \
    float: _Generic((b),             \
        fcomplex: _CFdivcr,          \
        default: _CFdiverr           \
    ),                               \
    default: _CFdiverr               \
)(a, b)


/* long double complex type ------------------------------ */

int cisequall(lcomplex a, lcomplex b);

long double creall(lcomplex z);
long double cimagl(lcomplex z);
long double cargl(lcomplex z);
long double cabsl(lcomplex z);
long double norml(lcomplex z);

lcomplex cminusl(lcomplex z);
lcomplex conjl(lcomplex z);
lcomplex cacosl(lcomplex z);
lcomplex casinl(lcomplex z);
lcomplex catanl(lcomplex z);
lcomplex ccosl(lcomplex z);
lcomplex csinl(lcomplex z);
lcomplex ctanl(lcomplex z);
lcomplex cacoshl(lcomplex z);
lcomplex casinhl(lcomplex z);
lcomplex catanhl(lcomplex z);
lcomplex ccoshl(lcomplex z);
lcomplex csinhl(lcomplex z);
lcomplex ctanhl(lcomplex z);
lcomplex cexpl(lcomplex z);
lcomplex clogl(lcomplex z);
lcomplex clog10l(lcomplex z);
lcomplex csqrtl(lcomplex z);
lcomplex cprojl(lcomplex z);
lcomplex csgnl(lcomplex z);

lcomplex _CLeulerr(long double x);
lcomplex _CLeulerc(lcomplex z);
void _CLeulererr(void);

#define ceulerl(z) _Generic((a),     \
    lcomplex:_CLeulerc,              \
    long double: _CLeulerr,          \
    default: _CLeulererr             \
)(z)

lcomplex _CLlognc(lcomplex z, lcomplex n);
lcomplex _CLlognr(lcomplex z, long double n);
lcomplex _CLlognz(long double z, lcomplex n);
void _CLlognerr(void);

#define clognl(a, b) _Generic((a),   \
    lcomplex: _Generic((b),          \
        lcomplex: _CLlognc,          \
        long double: _CLlognr,       \
        default: _CLlognerr          \
    ),                               \
    long double: _Generic((b),       \
        lcomplex: _CLlognz,          \
        default: _CLlognerr          \
    ),                               \
    default: _CLlognerr              \
)(a, b)

lcomplex _CLpowcc(lcomplex z, lcomplex n);
lcomplex _CLpowcr(lcomplex z, long double n);
lcomplex _CLpowrc(long double z, lcomplex n);
void _CLpowerr(void);

#define cpowl(a, b) _Generic((a),    \
    lcomplex: _Generic((b),          \
        lcomplex: _CLpowcc,          \
        long double: _CLpowcr,       \
        default: _CLpowerr           \
    ),                               \
    long double: _Generic((b),       \
        lcomplex: _CLpowrc,          \
        default: _CLpowerr           \
    ),                               \
    default: _CLpowerr               \
)(a, b)

lcomplex cbuildl(long double real, long double imag);

lcomplex _CLaddcc(lcomplex a, lcomplex b);
lcomplex _CLaddcr(lcomplex a, long double b);
lcomplex _CLaddrc(long double a, lcomplex b);
void _CLadderr(void);

#define caddl(a, b) _Generic((a),    \
    lcomplex: _Generic((b),          \
        lcomplex: _CLaddcc,          \
        long double: _CLaddcr,       \
        default: _CLadderr           \
    ),                               \
    long double: _Generic((b),       \
        lcomplex: _CLaddrc,          \
        default: _CLadderr           \
    ),                               \
    default: _CLadderr               \
)(a, b)

lcomplex _CLsubcc(lcomplex a, lcomplex b);
lcomplex _CLsubcr(lcomplex a, long double b);
lcomplex _CLsubrc(long double a, lcomplex b);
void _CLsuberr(void);

#define csubl(a, b) _Generic((a),    \
    lcomplex: _Generic((b),          \
        lcomplex: _CLsubcc,          \
        long double: _CLsubcr,       \
        default: _CLsuberr           \
    ),                               \
    long double: _Generic((b),       \
        lcomplex: _CLsubcr,          \
        default: _CLsuberr           \
    ),                               \
    default: _CLsuberr               \
)(a, b)

lcomplex _CLmulcc(lcomplex a, lcomplex b);
lcomplex _CLmulcr(lcomplex a, long double b);
lcomplex _CLmulrc(long double a, lcomplex b);
void _CLmulerr(void);

#define cmull(a, b) _Generic((a),    \
    lcomplex: _Generic((b),          \
        lcomplex: _CLmulcc,          \
        long double  : _CLmulcr,     \
        default : _CLmulerr          \
    ),                               \
    long double: _Generic((b),       \
        lcomplex: _CLmulrc,          \
        default : _CLmulerr          \
    ),                               \
    default: _CLmulerr               \
)(a, b)

lcomplex _CLdivcc(lcomplex a, lcomplex b);
lcomplex _CLdivcr(lcomplex a, long double b);
lcomplex _CLdivrc(long double a, lcomplex b);
void _CLdiverr(void);

#define cdivl(a, b) _Generic((a),    \
    lcomplex: _Generic((b),          \
        lcomplex: _CLdivcc,          \
        long double: _CLdivcr,       \
        default: _CLdiverr           \
    ),                               \
    long double: _Generic((b),       \
        lcomplex: _CLdivcr,          \
        default: _CLdiverr           \
    ),                               \
    default: _CLdiverr               \
)(a, b)

#endif
