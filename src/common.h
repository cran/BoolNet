#ifndef COMMON_H
#define COMMON_H
/**
 * Common structures and definitions
 */

// the number of bits that can be stored in one component of a 32-bit state array
#define BITS_PER_BLOCK_32 (sizeof(unsigned int) * 8)

// Retrieve the <i>-th bit in <x>
#define GET_BIT(x,i) (((x) & (1 << (i))) != 0)

// Set the <i>-th bit in <x> to 1
#define SET_BIT(x,i) ((x) | (1 << (i)))

// Set the <i>-th bit in <x> to 0
#define CLEAR_BIT(x,i) ((x) & (~(1 << (i))))

// Set the <i>-th bit in <x> to <v>
#define SET_BIT_TO_VAL(x,i,v) (((x) & (~(1 << (i)))) | ((v) << (i)))

/**
 * Encode a vector of binary values in an integer.
 * The rightmost element in <bin> is the leftmost bit in <dec>
 * <dec> is an array of <num> elements, and <bin> points
 * to an integer to which the result is written.
 */
extern void bin2dec(int *dec, int *bin, int *numBits);

/**
 * Decode an integer to a vector of binary values.
 * The rightmost element in <bin> is the leftmost bit in <dec>
 * <bin> points to the result vector, <dec> is a number
 * to be decoded, and <num> is the number of bits/elements in bin
 */
extern void dec2bin(int *bin, int *dec, int *numBits);

/**
 * Inserts values of fixed genes into states - this is required as
 * fixed genes are not encoded in the internal state representations.
 * <value> is a pointer to a state to be corrected.
 * <fixedGenes> is an array specifying which genes are fixed, as contained in
 * the BooleanNetwork structure.
 * <numGenes> is the length of <fixedGenes>.
 * The function changes the state pointed to by <value> and has no return value.
 */
extern void insertFixedGenes(unsigned int * value, int* fixedGenes, unsigned int numGenes);

#endif
