#ifndef COMMON_H
#define COMMON_H

#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include "uthash.h"

/**
 * Common structures and definitions
 */

/**
 * Macros for bit manipulation
 */

// the number of bits that can be stored in one component of a 32-bit state array
#define BITS_PER_BLOCK_32 (sizeof(unsigned int) * 8)

// Retrieve the <i>-th bit in <x>
#define GET_BIT(x,i) (((x) & ((unsigned long long)1 << (i))) != 0)

// Set the <i>-th bit in <x> to 1
#define SET_BIT(x,i) ((x) | ((unsigned long long)1 << (i)))

// Set the <i>-th bit in <x> to 0
#define CLEAR_BIT(x,i) ((x) & (~((unsigned long long)1 << (i))))

// Set the <i>-th bit in <x> to <v>
#define SET_BIT_TO_VAL(x,i,v) (((x) & (~((unsigned long long)1 << (i)))) | ((v) << (i)))

// Retrieve the <i>-th bit in an array of 32-bit integers <x>
#define GET_BIT_ARRAY(x,i) (((*(&(x) + i / BITS_PER_BLOCK_32)) & ((unsigned int)1 << (i % BITS_PER_BLOCK_32))) != 0)

// Set the <i>-th bit in an array of 32-bit integers <x> to 1
#define SET_BIT_ARRAY(x,i) (*(&(x) + i / BITS_PER_BLOCK_32) |= ((unsigned int)1 << (i % BITS_PER_BLOCK_32)))

// Set the <i>-th bit in an array of 32-bit integers <x> to 0
#define CLEAR_BIT_ARRAY(x,i) (*(&(x) + i / BITS_PER_BLOCK_32) &= (~((unsigned int)1 << (i % BITS_PER_BLOCK_32))))

/**
 * Macros for the management of vector-like buffers with automatic resizing
 */

// Allocate a buffer <buf> of type <type>[] and set the capacity counter <sz> and the number of elements <count>
#define ALLOC_BUFFER(buf, type, sz, count) do{sz = 4; count = 0; buf = calloc(sz, sizeof(type)); } while(0);

// Add a new element <el> of type <type> to the buffer <buf>, and update the capacity counter <sz> and the number of elements <count>
#define PUT_BUFFER(buf, type, sz, count, el) do{if (sz == count) {sz *= 2; buf = realloc(buf, sz * sizeof(type));} buf[count++] = (type) el; } while(0);

// Set the capacity <sz> of the buffer <buf> of type <type>[] to the current number of elements <count>
#define FINALIZE_BUFFER(buf, type, sz, count) do{sz = count; buf = realloc(buf, sz * sizeof(type)); } while(0);

/**
 * A hash structure for the allocated memory map
 */
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
static inline void* CALLOC(size_t n, size_t sz) 
{
  void * ptr = calloc(n, sz); 
  
  if (ptr == NULL)
    error("Out of memory!");
    
  AllocatedMemory * m = calloc(1, sizeof(AllocatedMemory));  
  m->ptr = ptr; 
  HASH_ADD_PTR(memoryMap, ptr, m);
  return ptr;
}

static inline void* REALLOC(void* ptr, size_t new_sz) 
{
  if (ptr == NULL)
    return CALLOC(new_sz, 1);
    
  void * newptr = realloc(ptr, new_sz); 
  
  if (newptr == NULL)
    error("Out of memory!");
  
  if (newptr != ptr)
  { 
    AllocatedMemory * m;
    HASH_FIND_PTR(memoryMap, &ptr, m); 
    HASH_DEL(memoryMap, m); 
    m->ptr = newptr; 
    HASH_ADD_PTR(memoryMap, ptr, m);
  }
  return newptr;
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
 * Add one to a binary number represented by an array of characters
 * with 0/1 elements. 
 * <fixed> is a vector of fixed positions that are not touched or NULL.
 * Return true if there is a next state, or false in case of an overflow.
 */
static inline bool getNextState(unsigned char * state, int * fixed, unsigned int numBits)
{
  if (numBits == 0)
    return false;
    
  int i = numBits - 1;
  
  while(true)
  {
    while (fixed != NULL && fixed[i] != -1)
    {     
      --i;
    }
    if (i < 0)
      return false;
      
    if (state[i] == 0)
    {
      state[i] = 1;
      return true;
    }
    else
    {
      if (i == 0)
        return false;
        
      state[i] = 0;
      --i;
    }
  }
}

/**
 * Free all remaining memory blocks stored in
 * the global map.
 * moved to init.c
 */
//extern void freeAllMemory();

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
 * Moved to init.c
 */
extern void bin2decC(int *dec, int *bin, int *numBits);

/**
 * Decode an integer to a vector of binary values.
 * The rightmost element in <bin> is the leftmost bit in <dec>
 * <bin> points to the result vector, <dec> is a number
 * to be decoded, and <num> is the number of bits/elements in bin
 * moved to init.c
 */
extern void dec2binC(int *bin, int *dec, int *numBits);

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

extern SEXP getListElement(SEXP list, char *str);

#endif
