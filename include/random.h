// 
// random.h
// 
//      Header of random functions
// 

#pragma once

#include <stdint.h>

/**
 * @brief Initializes the PCG random number generator with a seed.
 * @param seed The seed to initialize the random number generator. 
 *             Can use `(uint64_t)time(NULL)` as the seed.
 * @cite 
 * Wikipedia contributors. (2024, February 8). Permuted congruential generator. 
 * In Wikipedia, The Free Encyclopedia. Retrieved 06:20, June 2, 2024, from 
 * https://en.wikipedia.org/w/index.php?title=Permuted_congruential_generator&oldid=1204900506
 * @return A pseudo-random 32-bit integer.
 */
void pcg32_init(uint64_t seed);

/**
 * @brief Generates a pseudo-random 32-bit integer using the PCG algorithm.
 *        This function updates the internal state and produces a random number.
 *
 * @return A pseudo-random 32-bit integer.
 */
uint32_t pcg32(void);

/**
 * @brief Generates a random integer between the specified range [low, high].
 *
 * @param low The lower bound of the range (inclusive).
 * @param high The upper bound of the range (inclusive).
 * @return A random integer between low and high.
 */
int randint(int low, int high);

/**
 * @brief Generates a random boolean value.
 *
 * @return A random boolean value (true or false).
 */
bool randbool(void);

/**
 * @brief Generates a random real number in the range [0, 1].
 *
 * @return A random real number in the range [0, 1].
 */
double randreal(void);

/**
 * @brief Generates a random real number in the specified range [low, high].
 *
 * @param low The lower bound of the range.
 * @param high The upper bound of the range.
 * @return A random real number between low and high.
 */
double uniform(double low, double high);

/**
 * @brief Generates a random number following a normal (Gaussian) distribution.
 *
 * @param mu The mean of the distribution.
 * @param sigma The standard deviation of the distribution.
 * @return A random number following a normal distribution with mean `mu` and 
 *         standard deviation `sigma`.
 */
double normal(double mu, double sigma);

/**
 * @brief Generates a random number following a standard normal distribution 
 *        (mean = 0, sigma = 1).
 *
 * @return A random number following a standard normal distribution.
 */
double randn(void);
