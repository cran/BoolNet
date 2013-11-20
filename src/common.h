#ifndef COMMON_H
#define COMMON_H

#include <R.h>
#include "uthash.h"

/**
 * Common structures and definitions
 */

// the number of bits that can be stored in one component of a 32-bit state array
#define BITS_PER_BLOCK_32 (sizeof(unsigned int) * 8)

// Retrieve the <i>-th bit in <x>
#define GET_BIT(x,i) (((x) & ((unsigned int)1 << (i))) != 0)

// Set the <i>-th bit in <x> to 1
#define SET_BIT(x,i) ((x) | ((unsigned int)1 << (i)))

// Set the <i>-th bit in <x> to 0
#define CLEAR_BIT(x,i) ((x) & (~((unsigned int)1 << (i))))

// Set the <i>-th bit in <x> to <v>
#define SET_BIT_TO_VAL(x,i,v) (((x) & (~((unsigned int)1 << (i)))) | ((v) << (i)))

typedef struct
{
  // a pointer to the allocated memory
	void * ptr;

	// used by the hash table
	UT_hash_handle hh;
} AllocatedMemory;

// map that stores all allocated memory pointers
// to free them on interrupt
extern AllocatedMemory * memoryMap;

/**
 * Custom function to allocate memory that stores
 * the pointers in the global map.
 */
static inline void* CALLOC(unsigned int n, unsigned int sz) 
{
  void * ptr = calloc(n, sz); 
  AllocatedMemory * m = calloc(1, sizeof(AllocatedMemory)); 
  m->ptr = ptr; 
  HASH_ADD_PTR(memoryMap, ptr, m);
  return ptr;
}

/**
 * Custom function to free memory that was
 * allocated using CALLOC().
 */
static inline void FREE(void * ptr) 
{
  AllocatedMemory * m; 
  HASH_FIND_PTR(memoryMap, &ptr, m); 
  HASH_DEL(memoryMap, m); 
  free(m); 
  free(ptr);
}

/**
 * Free all remaining memory blocks stored in
 * the global map.
 */
extern void freeAllMemory();

/**
 * A structure that provides flexible
 * large array blocks as a list of arrays
 */
typedef struct ALE
{
  void * array;
    
  struct ALE * next;
} ArrayListElement;


static inline void allocNewArray(ArrayListElement ** head, unsigned int numElements, unsigned int elementSize)
{
  ArrayListElement * el = CALLOC(1, sizeof(ArrayListElement));
  el->array = CALLOC(numElements, elementSize);
  el->next = *head;
  *head = el;
}

/*
 * Free an array list <head>.
 */
static inline void freeArrayList(ArrayListElement * head)
{
  ArrayListElement * current = head;
  while (current != NULL)
	{
		ArrayListElement * next = current->next;
		FREE(current->array);
		FREE(current);
		current = next;
	}
}


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


/**
 * Removes values of fixed genes from states - this is required as
 * fixed genes are not encoded in the internal state representations.
 * <value> is a pointer to a state to be corrected.
 * <fixedGenes> is an array specifying which genes are fixed, as contained in
 * the BooleanNetwork structure.
 * <numGenes> is the length of <fixedGenes>.
 * The function changes the state pointed to by <value> and has no return value.
 */
extern void removeFixedGenes(unsigned int * value, int* fixedGenes, unsigned int numGenes);
#endif
