//	
// Jan 28, 2023
// Feb 01, 2023
// 
// mymath.c
// 

#include "mymath.h"

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
	_Dcomplex c = _Cbuild(0., 0.);
	c._Val[0] = a._Val[0] + b._Val[0];
	c._Val[1] = a._Val[1] + b._Val[1];
	return c;
}

_Dcomplex csub(_Dcomplex a, _Dcomplex b) {
	_Dcomplex c = _Cbuild(0., 0.);
	c._Val[0] = a._Val[0] - b._Val[0];
	c._Val[1] = a._Val[1] - b._Val[1];
	return c;
}

_Dcomplex cmul(_Dcomplex a, _Dcomplex b) {
	_Dcomplex c = _Cbuild(0., 0.);
	c._Val[0] = cabs(a) * cabs(b) * cos(carg(a) + carg(b));
	c._Val[1] = cabs(a) * cabs(b) * sin(carg(a) + carg(b));
	return c;
}

_Dcomplex cdiv(_Dcomplex a, _Dcomplex b) {
	_Dcomplex c = _Cbuild(0., 0.);
	c._Val[0] = cabs(a) / cabs(b) * cos(carg(a) - carg(b));
	c._Val[1] = cabs(a) / cabs(b) * sin(carg(a) - carg(b));
	return c;
}

// Numeric --------------------------------------------------------------------

bool isPowerOf2(unsigned num) {
	if (num && !(num & (num - 1))) return true;
	return false;
}

unsigned nextPowerOf2(unsigned num) {
	unsigned val = 1;
	if (isPowerOf2(num)) return num;
	while (val < num) val <<= 1;
	return val;
}

bool isprime(unsigned num) {
	if (num <= 1) return false;
	if (num == 2) return true;
	if (num % 2 == 0) return false;
	for (unsigned i = 3; i * i <= num; i += 2) {
		if (num % i == 0) return false;
	}
	return true;
}

unsigned* factor(unsigned num) {
	unsigned* result = 0;
	unsigned* tmp = 0;
	unsigned n = 2;
	if ((result = (unsigned*)calloc(n, sizeof(unsigned))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	if (num == 1) {
		result[0] = 0;
		result[n - 1] = 1;
		return result;
	}

	for (unsigned i = 2; i <= num; i++) {
		if (i != 2 && i % 2 == 0) i++;
		if (!isprime(i)) continue;
		if (isprime(num)) i = num;
		while (num != i) {
			if (num % i == 0) {
				tmp = (unsigned*)realloc(result, sizeof(unsigned) * n);
				if (tmp == NULL) return 0;
				result = tmp;
				tmp = NULL;
				result[n - 1] = i;
				num /= i;
				n += 1;
			}
			else break;
		}
	}

	tmp = (unsigned*)realloc(result, sizeof(unsigned) * n);
	if (tmp == NULL) return 0;
	result = tmp;
	tmp = NULL;
	result[0] = n;
	result[n - 1] = num;
	return result;
}

int extendedEuclidean(unsigned a, unsigned b, int* x, int* y) {
	if (b == 0) {
		*x = 1, * y = 0;
		return a;
	}
	int d = extendedEuclidean(b, a % b, y, x);
	*y -= a / b * *x;
	return d;
}

int* exgcd(unsigned a, unsigned b) {
	int x = 0, y = 0, g = 0;
	g = extendedEuclidean(a, b, &x, &y);

	int* r = 0;
	if ((r = (int*)calloc(4, sizeof(int))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	r[0] = 4, r[1] = x, r[2] = y, r[3] = g;
	return r;
}

unsigned expModuloN(unsigned g, unsigned k, unsigned n) {
	unsigned a = 1;
	while (k) {
		if (k & 1) a = (a * g) % n;
		k >>= 1;
		g = (g * g) % n;
	}
	return a;
}

unsigned expModuloNInverse(unsigned a, unsigned n) {
	int* arr = exgcd(a, n);
	int x = 0;
	if (arr[3] == 1 && arr[1] > 0) {
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

bool isPrimitiveRoot(unsigned g, unsigned p) {
	unsigned* r = 0;
	r = factor(p - 1);
	for (unsigned i = 1; i < r[0]; i++) {
		if (expModuloN(g, (p - 1) / r[i], p) == 1) {
			free(r);
			return false;
		}
	}
	free(r);
	return true;
}

unsigned findPrimitiveRoot(unsigned n) {
	if (!isprime(n)) {
		printf("n should be prime!\n");
		return 0;
	}
	for (unsigned rootCandidate = 2; rootCandidate < n; rootCandidate++) {
		if (isPrimitiveRoot(rootCandidate, n)) return rootCandidate;
	}
	printf("Primitive root not found");
	return 0;
}
