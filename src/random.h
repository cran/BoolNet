#ifndef RANDOM_H_
#define RANDOM_H_

/**
 * This header contains wrapper methods to generate random numbers.
 * The implementation uses a mersenne twister.
 */

#include "mersennetwister.h"
#include <stdlib.h>

/**
 * Sets the seed of the random number generator to <seed>
 */

static inline void setrandomseed(unsigned long seed)
{
	init_genrand(seed);
}

/**
 *  Returns a random double in [0,1)
 */
static inline double doublerand_1()
{
	return genrand_real2();
}


/**
 *  Returns a random double in [0,<maxVal>)
 */
static inline double doublerand(double maxVal)
{
	return doublerand_1() * maxVal;
}



/**
 * Returns a random integer value in [0,maxVal-1]
 */
static inline unsigned int intrand(unsigned int maxVal)
{
	return genrand_int32() % maxVal;
}

/**
 * Returns a random permutation of the numbers 0..<length>-1
 */
static inline unsigned int* randomPermutation(unsigned int length)
{
	unsigned int * res = (unsigned int *)malloc(length*sizeof(unsigned int));
	unsigned int i;
	for (i = 0; i < length; i++)
		res[i] = i;
	if (length <= 1)
		return res;
	for (i = length; i >= 0; i--)
	{
		unsigned int idx1 = intrand(length);
		unsigned int idx2 = intrand(length);
		unsigned int tmp = res[idx1];
		res[idx1] = res[idx2];
		res[idx2] = tmp;
	}
	return res;
}

#endif /*RANDOM_H_*/
