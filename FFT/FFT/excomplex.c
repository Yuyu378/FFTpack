// 
// excomplex.c
// 
//		Extend complex implementation
// 

#include "excomplex.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Functions
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

void cprint(_Dcomplex a) {
	if (a._Val[1] < 0) printf("%lf - %lfi\n", a._Val[0], -a._Val[1]);
	else printf("%lf + %lfi\n", a._Val[0], a._Val[1]);
	return;
}

void cswap(_Dcomplex* a, _Dcomplex* b) {
	_Dcomplex tmp = _Cbuild(0., 0.);
	tmp = *a;
	*a = *b;
	*b = tmp;
	return;
}

_Dcomplex cadd(_Dcomplex a, _Dcomplex b) {
	return _Cbuild(a._Val[0] + b._Val[0], a._Val[1] + b._Val[1]);
}

_Dcomplex csub(_Dcomplex a, _Dcomplex b) {
	return _Cbuild(a._Val[0] - b._Val[0], a._Val[1] - b._Val[1]);
}

_Dcomplex cmul(_Dcomplex a, _Dcomplex b) {
	return _Cmulcc(a, b);
}

_Dcomplex cdiv(_Dcomplex a, _Dcomplex b) {
	return p2c(pdiv(c2p(a), c2p(b)));
}

_Dcomplex p2c(_Dpolar p) {
	return _Cbuild(p._Val[0] * cos(p._Val[1]), p._Val[0] * sin(p._Val[1]));
}

_Dpolar _Pbuild(double _Radial, double _Angular) {
	_Dpolar p = { _Radial, _Angular };
	return p;
}

_Dpolar padd(_Dpolar a, _Dpolar b) {
	return c2p(cadd(p2c(a), p2c(b)));
}

_Dpolar psub(_Dpolar a, _Dpolar b) {
	return c2p(csub(p2c(a), p2c(b)));
}

_Dpolar pmul(_Dpolar a, _Dpolar b) {
	return _Pbuild(a._Val[0] * b._Val[0], a._Val[1] + b._Val[1]);
}

_Dpolar pdiv(_Dpolar a, _Dpolar b) {
	return _Pbuild(a._Val[0] / b._Val[0], a._Val[1] - b._Val[1]);
}

_Dpolar pexp(_Dpolar p, double n) {
	return _Pbuild(pow(p._Val[0], n), p._Val[1] * n);
}

_Dpolar c2p(_Dcomplex z) {
	return _Pbuild(cabs(z), carg(z));
}
