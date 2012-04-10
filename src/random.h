#ifndef RANDOM_H_
#define RANDOM_H_
#include <R.h>

/**
 * This header contains wrapper methods to generate random numbers.
 */

/**
 *  Returns a random double in [0,1)
 */
static inline double doublerand_1()
{
	return unif_rand();
}


/**
 * Returns a random integer value in [0,maxVal-1]
 */
static inline unsigned int intrand(unsigned int maxVal)
{
	return (unsigned int)(unif_rand() * maxVal);
}

#endif /*RANDOM_H_*/
