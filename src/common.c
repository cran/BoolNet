#include "common.h"
#include <string.h>

/**
 * Common utilities for the BoolNet package
 *
 * This is part of the BooleanNetwork R package.
 *
 * Copyright 2009/2010 by Christoph Müssel and Zhou Dao
 *
 * Contact christoph.muessel@uni-ulm.de
 */

/**
 * Encode a vector of binary values in an integer.
 * The rightmost element in <bin> is the leftmost bit in <dec>
 * <dec> is an array of <num> elements, and <bin> points
 * to an integer to which the result is written.
 */
void bin2dec(int *dec, int *bin, int *numBits)
{
	// clear output first
	unsigned int numElts;
	if (*numBits % BITS_PER_BLOCK_32 == 0)
		numElts = *numBits / BITS_PER_BLOCK_32;
	else
		numElts = *numBits / BITS_PER_BLOCK_32 + 1;

	memset(dec,0,numElts*sizeof(int));

	// decode input and write binary integers
	unsigned int * unsigned_dec = (unsigned int *) dec;
	unsigned int i;

	for(i = 0; i < *numBits; ++i)
	{
		unsigned_dec[i / BITS_PER_BLOCK_32] |= (bin[i] << (i % BITS_PER_BLOCK_32));
	}
}

/**
 * Decode an integer to a vector of binary values.
 * The rightmost element in <bin> is the leftmost bit in <dec>
 * <bin> points to the result vector, <dec> is a number
 * to be decoded, and <num> is the number of bits/elements in bin
 */
void dec2bin(int *bin, int *dec, int *numBits)
{
	unsigned int i;
	unsigned int * unsigned_dec = (unsigned int *) dec;

	for(i = 0; i < *numBits; ++i)
		if( (unsigned_dec[i / BITS_PER_BLOCK_32] & (1 << (i % BITS_PER_BLOCK_32))) != 0)
			bin[i] = 1;
		else
			bin[i] = 0;
}

/**
 * Inserts values of fixed genes into states - this is required as
 * fixed genes are not encoded in the internal state representations.
 * <value> is a pointer to a state to be corrected.
 * <fixedGenes> is an array specifying which genes are fixed, as contained in
 * the BooleanNetwork structure.
 * <numGenes> is the length of <fixedGenes>.
 * The function changes the state pointed to by <value> and has no return value.
 */
void insertFixedGenes(unsigned int * value, int* fixedGenes, unsigned int numGenes)
{
	unsigned int tmp[numGenes];
	unsigned int i, j = 0;

	// build an array of Boolean values for the genes
	for (i = 0; i < numGenes; ++i)
	{
		if (fixedGenes[i] != -1)
		// this gene is fixed
		{
			tmp[i] = fixedGenes[i];
		}
		else
		// not a fixed gene => take value from original state
		{
			tmp[i] = ((value[j / BITS_PER_BLOCK_32] & (1 << (j % BITS_PER_BLOCK_32))) != 0) ? 1 : 0;
			++j;
		}
	}

	// re-encode Boolean array to integer value
	bin2dec((int *)value,(int*)tmp,(int*)&numGenes);
}

/**
 * Removes values of fixed genes from states - this is required as
 * fixed genes are not encoded in the internal state representations.
 * <value> is a pointer to a state to be corrected.
 * <fixedGenes> is an array specifying which genes are fixed, as contained in
 * the BooleanNetwork structure.
 * <numGenes> is the length of <fixedGenes>.
 * The function changes the state pointed to by <value> and has no return value.
 */
void removeFixedGenes(unsigned int * value, int* fixedGenes, unsigned int numGenes)
{	
	unsigned int tmp[numGenes];
	memset(tmp,0,sizeof(unsigned int) * numGenes);
	unsigned int i, j = 0;
	
	// build an array of Boolean values for the genes
	for (i = 0; i < numGenes; ++i)
	{
		if (fixedGenes[i] == -1)
		{
			tmp[j] = ((value[i / BITS_PER_BLOCK_32] & (1 << (i % BITS_PER_BLOCK_32))) != 0) ? 1 : 0;
			++j;
		}
	}

	// re-encode Boolean array to integer value
	bin2dec((int *)value,(int*)tmp,(int*)&numGenes);
}
