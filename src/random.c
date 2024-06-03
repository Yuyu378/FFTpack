// 
// random.c
// 
//      random functions
// 

#include "random.h"

#include <math.h>
#include <stdint.h>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#if defined(M_PI) && !defined(PI)
#define PI M_PI
#endif

#ifndef PI
#define PI (3.1415926535897932384626433832795028841971693993751058209749445923)
#endif

static uint64_t state = 0x4d595df4d0f33173;		            // Or something seed-dependent
static uint64_t const multiplier = 6364136223846793005u;
static uint64_t const increment  = 1442695040888963407u;	// Or an arbitrary odd constant

/**
 * @brief Rotates a 32-bit integer right by a specified number of bits.
 *
 * @param x The 32-bit integer to rotate.
 * @param r The number of bits to rotate by.
 * @return The rotated 32-bit integer.
 */
static uint32_t rotr32(uint32_t x, unsigned r) {
	return x >> r | x << (-r & 31);
}

void pcg32_init(uint64_t seed) {
	state = seed + increment;
	(void)pcg32();
}

uint32_t pcg32(void) {
	uint64_t x = state;
	unsigned count = (unsigned)(x >> 59);		// 59 = 64 - 5
	state = x * multiplier + increment;
	x ^= x >> 18;								// 18 = (64 - 27)/2
	return rotr32((uint32_t)(x >> 27), count);	// 27 = 32 - 5
}

int randint(int low, int high) {
    return low + (int)(randreal()  * (high - low + 1));
}

_Bool randbool(void) {
    return (_Bool)randint(0, 1);
}

double randreal(void) {
    return pcg32() * (1.0 / 0xffffffff); 
    /* divided by 2^32-1 */
}

double uniform(double low, double high) {
    return low + randreal() * (high - low + 1.0);
}

double normal(double mu, double sigma) {
    double u1, u2;
    do {
        u1 = randreal();
    } while (u1 == 0);
    u2 = randreal();
    double mag = sigma * sqrt(-2. * log(u1));
    double omg = 2 * PI * u2;
    return mag * (randbool() ? cos(omg) : sin(omg)) + mu;
}

double randn(void) {
    return normal(0.0, 1.0);
}
