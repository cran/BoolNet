/**
 * C code for the reconstruction of Boolean networks
 *
 * This is part of the BoolNet R package.
 *
 * Copyright 2009/2010 by Christoph Müssel
 *
 * Contact christoph.muessel@uni-ulm.de
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "common.h"

/**
 * Structure for lists of
 * transition functions for
 * the return value of reconstruction functions
 */
typedef struct FLE
{
	// the number of input genes
	unsigned int k;

	// the indices of the input genes
	int * inputGenes;

  // The transition function.
  // ATTENTION:
  // If all functions are built (returnPBN = true), this
	// is a bit vector containing the transition function
	// If "don't care" values are allowed (returnPBN = false),
	// each integer element encodes one bit
	unsigned int * transitionFunction;

	// the next element in the list
	struct FLE * next;
} FunctionListElement;

/**
 * Structure to build an internal stack
 * of functions. This is used to determine
 * all equally-rated functions recursively
 */
typedef struct FSE
{
	// the bit position to be processed next
	unsigned int pos;

	// the (incomplete) transition function
	unsigned int * transitionFunction;

	// the next element on the stack
	struct FSE * next;
} FunctionStackElement;


typedef struct 
{
  unsigned int pos;
  unsigned int numFixed;
  unsigned int numAvailable;
  unsigned int k;
  unsigned int n;
  unsigned int * indexMapping;
  int * comb;
  int * intComb;
} InputCombination;


/**
 * Add a new element to the function list specified by <root>.
 * Set the corresponding values <k>,<inputGenes> and <transitionFunction>
 * by copying the supplied values.
 * <transitionFunctionSize> is the number of elements needed for the transition function vector.
 * Returns the new list element.
 */
static inline FunctionListElement * addFunctionListElement(FunctionListElement ** root,
											 unsigned int k,
											 unsigned int transitionFunctionSize,
											 int * inputGenes,
											 unsigned int * transitionFunction)
{
	FunctionListElement * el = CALLOC(1,sizeof(FunctionListElement));
	el->k = k;

	el->inputGenes = CALLOC(k,sizeof(unsigned int));
	memcpy(el->inputGenes,inputGenes,sizeof(unsigned int) * k);

	el->transitionFunction = CALLOC(transitionFunctionSize,sizeof(unsigned int));
	memcpy(el->transitionFunction,transitionFunction,sizeof(unsigned int) * transitionFunctionSize);

	el->next = *root;
	*root = el;
	return el;
}

/**
 * Free all elements of the list <root>
 */
void freeFunctionList(FunctionListElement ** root)
{
	if (*root == NULL)
		return;
	FunctionListElement * current = *root;
	do
	{
		FunctionListElement * next = current->next;
		FREE(current->inputGenes);
		FREE(current->transitionFunction);
		FREE(current);
		current = next;
	}
	while(current != NULL);
	*root = NULL;
}

/**
 * Push a new element on top of the stack <stack>.
 * <transitionFunction> and <pos> describe the data to be pushed.
 * These elements are copied.
 * <transitionFunctionSize> is the number of elements of <transitionFunction>.
 * Returns the new stack element.
 */
static inline FunctionStackElement * pushFunctionStackElement(FunctionStackElement ** stack,
											 unsigned int * transitionFunction,
											 unsigned int transitionFunctionSize,
											 unsigned int pos)
{
	FunctionStackElement * el = CALLOC(1,sizeof(FunctionStackElement));

	el->pos = pos;

	el->transitionFunction = CALLOC(transitionFunctionSize,sizeof(unsigned int));
	memcpy(el->transitionFunction,transitionFunction,sizeof(unsigned int) * transitionFunctionSize);

	el->next = *stack;
	*stack = el;
	return el;
}

/**
 * Remove the top-level element from <stack>.
 */
static inline void deleteFunctionStackElement(FunctionStackElement ** stack)
{
	FunctionStackElement * el = *stack;
	*stack = (*stack)->next;
	FREE(el->transitionFunction);
	FREE(el);
}


/**
 * A function that iteratively returns all <n> choose <k> combinations of input genes
 * if called multiple times.
 * The function receives an input array <comb> containing the previous combination.
 * Before the first call, this combination must be initialized to {k-1, k-2, ... , 0},
 * and <pos> must be set to 0.
 * After each call, <comb> contains the next combination. If all combinations have been
 * listed, the function returns false, otherwise true.
 */
static inline bool nextCombination(InputCombination * comb)
{
	bool posChanged = false;

	// find the first position that has not been set
	// to its maximum number
	while(comb->pos < comb->k - comb->numFixed && 
	      comb->intComb[comb->pos] == comb->numAvailable - comb->pos - 1)
	{
		++ (comb->pos);
		posChanged = true;
	}
	if (comb->pos == comb->k - comb->numFixed)
	{
		// all elements have been listed
		return false;
  }
	if (posChanged)
	// reset lower-order positions, and increase
	// the position found previously
	{
		unsigned int i;
		++comb->intComb[comb->pos];
		for (i = comb->pos; i > 0; --i)
			comb->intComb[i-1] = comb->intComb[i] + 1;
		comb->pos = 0;
	}
	else
	// the current position has not been set to its maximum value
	// => increase it
		++comb->intComb[comb->pos];
	
	unsigned int j;	
	for (j = 0; j < comb->k - comb->numFixed; ++j)
	{
		comb->comb[comb->numFixed + j] = comb->indexMapping[comb->intComb[j]];
  }
  
  //for (j = 0; j < comb->k; ++j)
  //  Rprintf("%d ",comb->comb[j]);
  //Rprintf("\n");  
    
	
	return true;
}

static inline InputCombination * initCombination(int * requiredDependencies, 
                                                 int * excludedDependencies, 
                                                 unsigned int k, unsigned int n)
{
  InputCombination * res = CALLOC(1,sizeof(InputCombination));
  res->comb = CALLOC(k,sizeof(int));
  res->indexMapping = CALLOC(n,sizeof(int));
  
  res->k = k;
  res->n = n;
  res->numFixed = 0;
  res->numAvailable = 0;
  
  unsigned int j;
  
  for (j = 0; j < n; ++j)
  {
    if (requiredDependencies != NULL && requiredDependencies[j] != 0)
    {
      res->comb[res->numFixed] = j;
      ++(res->numFixed);
    }
    else
    if (excludedDependencies == NULL || excludedDependencies[j] == 0)
    {
      res->indexMapping[res->numAvailable] = j;
      ++(res->numAvailable);
    }
  }
  
  res->intComb = CALLOC(res->numAvailable,sizeof(unsigned int));
	for (j = 0; j < res->k - res->numFixed; ++j)
	{
	  res->intComb[j] = k - res->numFixed - j - 1;
		res->comb[res->numFixed + j] = res->indexMapping[res->intComb[j]];
  }
  return res;
}

static inline void freeInputCombination(InputCombination * comb)
{
  FREE(comb->comb);
  FREE(comb->intComb);  
  FREE(comb->indexMapping);
  FREE(comb);
}

/**
 * Lähdesmäki's best-fit extension algorithm to infer Boolean networks
 * from time series.
 * <inputStates> and <outputStates> represent state transitions. 
 * Here, <ceil(numGenes / 32)> consecutive array entries describe one state, 
 * thus the array size is <ceil(numGenes / 32) * numberOfStates>.
 * <perturbations> is an array of the same length as <inputStates> specifying
 * whether inputs were normally regulated (-1), overexpressed (1), or knocked out (0).
 * <numGenes> is the number of genes in the time series.
 * <maxK> is the maximum number of inputs a function can have.
 * <requiredDependencies> is a flattened Boolean matrix specifying dependencies
 * that *must* occur in the reconstructed networks.
 * Similarly, <excludedDependencies> specifies dependencies
 * that *must not* occur in the reconstructed networks. 
 * <result> receives one list of possible functions for each gene, and
 * <errors> stores the error of these functions on the time series.
 * <allSolutions> specifies whether the algorithm stops as soon as a perfect solution
 * is found (false) or not (true).
 * <returnPBN> specifies whether incomplete functions (with "don't care" values)
 * are expanded recursively (true) or returned as-is (false).
 */
void bestFitExtension(unsigned int * inputStates, 
            unsigned int * outputStates,
            unsigned int * perturbations,
					  unsigned int numStates, 
					  unsigned int numGenes, 
            int * requiredDependencies, 
            int * excludedDependencies,					  
					  unsigned int maxK,
					  FunctionListElement ** result, 
					  unsigned int * errors, 
					  bool allSolutions,
					  bool returnPBN)
{
	unsigned int i, j;
	unsigned int numElts;

	// calculate block size in time series array
	if (numGenes % BITS_PER_BLOCK_32 == 0)
		numElts = numGenes / BITS_PER_BLOCK_32;
	else
		numElts = numGenes / BITS_PER_BLOCK_32 + 1;

	// store the lengths of the current best solutions for
	// each gene - initialized with maximum value.
	unsigned int bestLength[numGenes];
	memset(bestLength,0xFF,sizeof(unsigned int) * numGenes);

	for (i = 0; i < numGenes; ++i)
	// iterate over genes
	{
		unsigned int currentMaxK = maxK, minK = 0, excludedCount = 0;
    
    if (requiredDependencies != NULL)
      for (j = 0; j < numGenes; ++j)
        if (requiredDependencies[i*numGenes + j] != 0)
          ++minK;
    
   if (excludedDependencies != NULL)
      for (j = 0; j < numGenes; ++j)
        if (excludedDependencies[i*numGenes + j] != 0)
          ++excludedCount;

    if (currentMaxK < minK)
      currentMaxK = minK;
      
    if (currentMaxK > numGenes - excludedCount)
      currentMaxK =  numGenes - excludedCount;
	
		// set error to maximum
		errors[i] = ~0;
		unsigned int k, s;

    if (minK == 0)
    {
		  // check for constant genes
		  unsigned int geneVal = GET_BIT(outputStates[i/BITS_PER_BLOCK_32],i % BITS_PER_BLOCK_32);
		  bool isConst = true;

		  unsigned int const0Err = (geneVal > 0), const1Err = (geneVal == 0);

		  for (s = 1; s < numStates; ++s)
		  {
    		if (perturbations != NULL && perturbations[s * numGenes + i] != -1)
    		// ignore this state if the current gene is perturbed
    		  continue;
    		  
			  unsigned int nextBit = GET_BIT(outputStates[s*numElts + i/BITS_PER_BLOCK_32],i % BITS_PER_BLOCK_32);
			  if (nextBit	!= geneVal)
			  {
				  isConst = false;
			  }
			  if (nextBit)
				  ++const0Err;
			  else
				  ++const1Err;
		  }

		  if (isConst)
		  // gene is constant => simplest function already found!
		  {
			  int inputGenes = -1;
			  addFunctionListElement(&result[i],1,1,&inputGenes,&geneVal);
			  errors[i] = 0;
			  bestLength[i] = 0;
		  }
		  else
		  {
			  if (const0Err <= const1Err)
			  {
				  int inputGenes = -1;
				  unsigned int val = 0;
				  addFunctionListElement(&result[i],1,1,&inputGenes,&val);
				  errors[i] = const0Err;
				  bestLength[i] = 0;
			  }

			  if (const1Err <= const0Err)
			  {
				  int inputGenes = -1;
				  unsigned int val = 1;
				  addFunctionListElement(&result[i],1,1,&inputGenes,&val);
				  errors[i] = const1Err;
				  bestLength[i] = 0;
			  }
		  }
    }
    
		for (k = (minK > 1? minK : 1); k <= currentMaxK; ++k)
		// iterate over possible numbers of inputs
		{
			if (errors[i] == 0 && !allSolutions)
				break;

			// initialize gene combination vector
			
			InputCombination * comb = initCombination((requiredDependencies == NULL? 
			                                            NULL : 
			                                            &requiredDependencies[i*numGenes]), 
                                                (excludedDependencies == NULL?
                                                  NULL :
                                                  &excludedDependencies[i*numGenes]), 
                                                k, numGenes);
			
			// calculate the number of array elements needed
			// to represent the output of a function
			unsigned int array_sz = (unsigned int)1 << k;
			unsigned int numEltsFunc;
			if (array_sz % BITS_PER_BLOCK_32 == 0)
				numEltsFunc = array_sz / BITS_PER_BLOCK_32;
			else
				numEltsFunc = array_sz / BITS_PER_BLOCK_32 + 1;

			unsigned int c_0[array_sz];
			unsigned int c_1[array_sz];

			do
			{
  			R_CheckUserInterrupt();
			
				memset(c_0,0,sizeof(unsigned int) * array_sz);
				memset(c_1,0,sizeof(unsigned int) * array_sz);

				for (s = 0; s < numStates; ++s)
				// iterate over states and count errors
				{
					if (perturbations != NULL && perturbations[s * numGenes + i] != -1)
				  // the current gene is perturbed in this state 
				  // => ignore the state for inference of this gene					
    		    continue;
				
					unsigned int input = 0, bit;
					for (bit = 0; bit <  k; ++bit)
					// build input by transferring the bits of the input genes
					// in the current state into a condensed input value
					{

						input |= (GET_BIT(inputStates[s * numElts
										   + comb->comb[bit]/BITS_PER_BLOCK_32],
										   comb->comb[bit] % BITS_PER_BLOCK_32))
								  << bit;
					}


					// determine the response value after a state transition
					unsigned int response = GET_BIT(outputStates[s * numElts
													 + i/BITS_PER_BLOCK_32],i % BITS_PER_BLOCK_32);

					if (response == 0)
						++c_1[input];
					else
						++c_0[input];
				}

				unsigned int f[numEltsFunc];
				memset(f,0,sizeof(unsigned int) * numEltsFunc);

				unsigned int error = 0, c;
				for (c = 0; c < array_sz; ++c)
				{
					if (c_0[c] > c_1[c])
						error += c_1[c];
					else
						error += c_0[c];

				}

				if (error < errors[i])
				// the new solution is better than all previously found solutions
				// => reset the solution list and the best error
				{
					errors[i] = error;
					bestLength[i] = k;
					freeFunctionList(&result[i]);
				}

				if (error <= errors[i] && (bestLength[i] >= k || allSolutions))
				{

          if (!returnPBN)
          // encode a function with "don't care" value
          // note: in contrast to the code below, 
          // each entry of f corresponds to a single bit here          
          {
   					unsigned int f[array_sz];
  					memset(f,0,sizeof(unsigned int) * array_sz);
					  for (unsigned int l = 0; l < array_sz; ++l)
					  {
      				if (c_1[l] < c_0[l])
      				{
      				  f[l] = 1;
      				}
      				else
      				if (c_1[l] > c_0[l])
      				{
      				  f[l] = 0;
      				}
      				else
      				{
      				  f[l] = -1;
      				}
					  }					
					  addFunctionListElement(&result[i],k,array_sz,comb->comb,f);
          }
          else
          {
            // recursively determine all functions by filling "don't care" values.
            // note: here, the elements of f are used as bit vectors
            // to save memory

					  // start with an element initialized to zero,
					  // and first examine the 0-th position

	  				unsigned int f[numEltsFunc];
  					memset(f,0,sizeof(unsigned int) * numEltsFunc);
   					FunctionStackElement * stack =  NULL;
					  pushFunctionStackElement(&stack,f,numEltsFunc,0);
					  do
					  {
					    R_CheckUserInterrupt();
					    
						  // get top-level stack element
						  FunctionStackElement * el = stack;

						  if (el->pos == array_sz)
						  // the top-level element on the stack is a complete function
						  // => clean it up
						  {
							  // create a new element in the result function list containing the
							  // completed function on the stack
							  addFunctionListElement(&result[i],k,numEltsFunc,
							                         comb->comb,el->transitionFunction);
							  // remove the completed function from the stack
							  deleteFunctionStackElement(&stack);
						  }
						  else
						  if (c_1[el->pos] == c_0[el->pos])
						  // the 0-error and the 1-error are equal
						  // => create two solution branches, one with the <pos>-th bit set to 1
						  // and one with the <pos>-th bit set to 0
						  {
							  // create a new element on the stack with the <pos>-th bit set to 1
							  FunctionStackElement * new = pushFunctionStackElement(&stack,el->transitionFunction,
																				    numEltsFunc,el->pos+1);
							  new->transitionFunction[el->pos / BITS_PER_BLOCK_32] |= (unsigned int)1 << (el->pos % BITS_PER_BLOCK_32);

							  // increment the position of the old element, which remains on the stack
							  // with the <pos>-th bit set to 0 (due to initialization)
							  ++el->pos;
						  }
						  else
						  if (c_1[el->pos] < c_0[el->pos])
						  // the error is lower if the <pos>-th bit is set to 1
						  {
							  // set the <pos>-th bit of the top-level stack element to 1,
							  // and increment the position to be examined
							  el->transitionFunction[el->pos / BITS_PER_BLOCK_32] |= (unsigned int)1 << (el->pos % BITS_PER_BLOCK_32);
							  ++el->pos;
						  }
						  else
						  // the error is lower if the <pos>-th bit is set to 0
						  {
							  // due to initialization, the <pos>-th bit is already set to 0
							  // => increment the position to be examined
							  ++el->pos;
						  }
					  }
					  // terminate if all stack elements have been completed
					  while (stack != NULL);
					}
				}

			}
			// get next input gene combination vector
			while(nextCombination(comb));
      
      freeInputCombination(comb);
		}
	}
}

/**
 * Calculate the entropy of a set of <numChosenIndices> genes (gene_1,...,gene_n)
 * specified by <chosenIndices>.
 * All values less than <numGenes> specify genes in the input states, and values greater or
 * equal to <numGenes> specify genes in the output states.
 * <inputStates> and <outputStates> contain <numStates> states, each consisting
 * of <elementsPerEntry> array elements.
 * <perturbations> is an array of the same length as <inputStates> specifying
 * whether inputs were normally regulated (-1), overexpressed (1), or knocked out (0).
 * The table counting the occurrences of gene value combinations is returned in <table>, which
 * must be initialized to a size of 2^<numChosenIndices>.
 * The return value is the entropy H(gene1,...,gene_n).
 */
static inline double entropy(
              unsigned int * inputStates, 
              unsigned int * outputStates,
              unsigned int * perturbations,
              unsigned int currentGene,
				      unsigned int numStates, unsigned int elementsPerEntry,
				      unsigned int numGenes, int * chosenIndices,
				      unsigned int numChosenIndices,
				      unsigned int * table)
{
	unsigned int numEntries = (unsigned int)1 << numChosenIndices;

	// reset table to 0
	memset(table,0,sizeof(unsigned int) * numEntries);
	unsigned int state, validStates = 0;

	for (state = 0; state < numStates; ++state)
	{
		if (perturbations != NULL && perturbations[state * numGenes + currentGene] != -1)
    // the current gene is perturbed in this state 
    // => ignore the state in the entropy calculation
	    continue;
	  
	  ++validStates;
	    
		unsigned int tableIndex = 0, geneIndex;
		for (geneIndex = 0; geneIndex < numChosenIndices; ++geneIndex)
		// count the number of occurrences of gene value combinations
		{
			// encode the index of the table element to be increased
			if (chosenIndices[geneIndex] < numGenes)
			// this is a gene in the input states
			{
				unsigned int chosenIndex = chosenIndices[geneIndex];
				tableIndex |= (GET_BIT(inputStates[state * elementsPerEntry
								  + chosenIndex/BITS_PER_BLOCK_32],chosenIndex % BITS_PER_BLOCK_32))
								  << geneIndex;
			}
			else
			// this is a gene in the output states
			{

				unsigned int chosenIndex = chosenIndices[geneIndex]  - numGenes;
				tableIndex |= (GET_BIT(outputStates[state * elementsPerEntry
							  + chosenIndex/BITS_PER_BLOCK_32],chosenIndex % BITS_PER_BLOCK_32))
							  << geneIndex;
			}
		}
		// increase the corresponding table entry
		++table[tableIndex];
	}

	// calculate entropy
	double result = 0.0;
	unsigned int tableIndex;
	for (tableIndex = 0; tableIndex < numEntries; ++tableIndex)
	{
		if (table[tableIndex] != 0)
			result += ((double)table[tableIndex])/validStates*log2(((double)table[tableIndex])/validStates);
	}
	return -result;
}

/**
 * Liang's REVEAL algorithm for reconstruction of Boolean networks.
 * This version is enhanced for dealing with inconsistent samples by
 * taking the solutions that produce the minimum error.
 * <inputStates> and <outputStates> represent state transitions. 
 * Here, <ceil(numGenes / 32)> consecutive array entries describe one state, 
 * thus the array size is <ceil(numGenes / 32) * numberOfStates>.
 * <perturbations> is an array of the same length as <inputStates> specifying
 * whether inputs were normally regulated (-1), overexpressed (1), or knocked out (0).
 * <numGenes> is the number of genes in the time series.
 * <maxK> is the maximum number of inputs a function can have.
 * <requiredDependencies> is a flattened Boolean matrix specifying dependencies
 * that *must* occur in the reconstructed networks.
 * Similarly, <excludedDependencies> specifies dependencies
 * that *must not* occur in the reconstructed networks. 
 * <result> receives one list of possible functions for each gene, and
 * <errors> stores the error of these functions on the time series.
 * <allSolutions> specifies whether the algorithm stops as soon as a perfect solution
 * is found (false) or not (true).
 * <returnPBN> specifies whether incomplete functions (with "don't care" values)
 * are expanded recursively (true) or returned as-is (false).
 */
void reveal(unsigned int * inputStates, 
            unsigned int * outputStates,
			      unsigned int * perturbations,            
			      unsigned int numStates, 
			      unsigned int numGenes, 
            int * requiredDependencies, 
            int * excludedDependencies,	
			      unsigned int maxK,
			      FunctionListElement ** result, 
			      unsigned int * errors, 
			      bool allSolutions,
			      bool returnPBN)
{
	unsigned int i, j;
	unsigned int numElts;

	// calculate block size in time series array
	if (numGenes % BITS_PER_BLOCK_32 == 0)
		numElts = numGenes / BITS_PER_BLOCK_32;
	else
		numElts = numGenes / BITS_PER_BLOCK_32 + 1;

	for (i = 0; i < numGenes; ++i)
	// iterate over genes
	{
		unsigned int currentMaxK = maxK, minK = 0, excludedCount = 0;
    
    if (requiredDependencies != NULL)
      for (j = 0; j < numGenes; ++j)
        if (requiredDependencies[i*numGenes + j] != 0)
          ++minK;
          
   if (excludedDependencies != NULL)
      for (j = 0; j < numGenes; ++j)
        if (excludedDependencies[i*numGenes + j] != 0)
          ++excludedCount;      
    
    if (currentMaxK < minK)
      currentMaxK = minK;
      
    if (currentMaxK > numGenes - excludedCount)
      currentMaxK =  numGenes - excludedCount;
	
		// set error to maximum
		// Note: in REVEAL, <errors> is only used as an indicator
		// whether a solution has been found, i.e. it is 0 if
		// a solution was found, and anything else otherwise.
		errors[i] = ~0;
		unsigned int k, s;

		// check for constant genes
		unsigned int geneVal = GET_BIT(outputStates[i/BITS_PER_BLOCK_32],i % BITS_PER_BLOCK_32);
		bool isConst = true;

		for (s = 1; s < numStates; ++s)
		{
			if ((perturbations == NULL || perturbations[s * numGenes + i] == -1) &&
			    GET_BIT(outputStates[s*numElts + i/BITS_PER_BLOCK_32],i % BITS_PER_BLOCK_32)
          				!= geneVal)
  		// ignore this state if the current gene is perturbed          				
			{
				isConst = false;
				break;
			}
		}

		if (isConst)
		// gene is constant => simplest function already found!
		{
			int inputGenes = -1;
			addFunctionListElement(&result[i],1,1,&inputGenes,&geneVal);
			errors[i] = 0;
		}

		for (k = (minK > 1? minK : 1); k <= currentMaxK; ++k)
		// iterate over possible numbers of inputs
		{
			if (errors[i] == 0 && !allSolutions)
				break;

			// initialize gene combination vector
			InputCombination * comb = initCombination((requiredDependencies == NULL? 
			                                            NULL : 
			                                            &requiredDependencies[i*numGenes]), 
                                                (excludedDependencies == NULL?
                                                  NULL :
                                                  &excludedDependencies[i*numGenes]), 
                                                k, numGenes);

			// calculate the number of array elements needed
			// to represent the output of a function
			unsigned int array_sz = (unsigned int)1 << k;
			unsigned int numEltsFunc;
			if (array_sz % BITS_PER_BLOCK_32 == 0)
				numEltsFunc = array_sz / BITS_PER_BLOCK_32;
			else
				numEltsFunc = array_sz / BITS_PER_BLOCK_32 + 1;


			do
			{
        R_CheckUserInterrupt();
        
				// calculate entropy of input genes
				unsigned int table_input[array_sz];
				double entropy_input = entropy(inputStates, outputStates,  
				                               perturbations, i, numStates,
              											   numElts, numGenes, comb->comb, k, table_input);

				// calculate entropy of the combination of input genes and output value
				int comb_output[k+1];
				memcpy(comb_output,comb->comb,sizeof(unsigned int) * k);
				comb_output[k] = numGenes + i;
				unsigned int table_output[(unsigned int)1 << (k+1)];

				double entropy_all = entropy(inputStates, outputStates, 
  		                               perturbations, i, numStates,				
	             										   numElts, numGenes, comb_output, k+1, table_output);
 
				if (fabs(entropy_input - entropy_all) < 1E-6)
				// the two entropies are the same => output is fully explained by input
				{

					// a solution has been found => set errors to 0 so that
					// the algorithm does not search for higher k's
					errors[i] = 0;

          if (!returnPBN)
          // encode a function with "don't care" value
          // note: in contrast to the code below, 
          // each entry of f corresponds to a single bit here
          {
   					unsigned int f[array_sz];
  					memset(f,0,sizeof(unsigned int) * array_sz);

					  for (unsigned int l = 0; l < array_sz; ++l)
					  {
					    if (table_input[l] == 0 || table_output[l | ((unsigned int)1 << k)] == table_output[l])
      				{
      				  f[l] = -1;
      				}
      				else
      				if (table_output[l | ((unsigned int)1 << k)] > table_output[l])
      				{
      				  f[l] = 1;
      				}
      				else
      				{
      				  f[l] = 0;
      				}
					  }					
					  addFunctionListElement(&result[i],k,array_sz,comb->comb,f);
          }
          else
          {
            // recursively determine all functions by filling "don't care" values.
            // note: here, the elements of f are used as bit vectors
            // to save memory
   					// start with an element initialized to zero,
  					// and first examine the 0-th position
  					unsigned int f[numEltsFunc];
  					memset(f,0,sizeof(unsigned int) * numEltsFunc);
   					FunctionStackElement * stack =  NULL;
					  pushFunctionStackElement(&stack,f,numEltsFunc,0);
					  do
					  {
					    R_CheckUserInterrupt();
					
						  // get top-level stack element
						  FunctionStackElement * el = stack;

						  if (el->pos == array_sz)
						  // the top-level element on the stack is a complete function
						  // => clean it up
						  {
							  // create a new element in the result function list containing the
							  // completed function on the stack
							  addFunctionListElement(&result[i],k,numEltsFunc,
							                         comb->comb,el->transitionFunction);
							  // remove the completed function from the stack
							  deleteFunctionStackElement(&stack);
						  }
						  else
						  if (table_input[el->pos] == 0 || table_output[el->pos | ((unsigned int)1 << k)] == table_output[el->pos])
						  // no information is available if the <pos>-th bit must be set to one or zero
						  // => create two solution branches, one with the <pos>-th bit set to 1
						  // and one with the <pos>-th bit set to 0
						  {
							  // create a new element on the stack with the <pos>-th bit set to 1
							  FunctionStackElement * new = pushFunctionStackElement(&stack,el->transitionFunction,
																				    numEltsFunc,el->pos+1);
							  new->transitionFunction[el->pos / BITS_PER_BLOCK_32] |= (unsigned int)1 << (el->pos % BITS_PER_BLOCK_32);

							  // increment the position of the old element, which remains on the stack
							  // with the <pos>-th bit set to 0 (due to initialization)
							  ++el->pos;
						  }
						  else
						  if (table_output[el->pos | ((unsigned int)1 << k)] > table_output[el->pos])
						  // the <pos>-th bit must be set to 1
						  {
							  // set the <pos>-th bit of the top-level stack element to 1,
							  // and increment the position to be examined
							  el->transitionFunction[el->pos / BITS_PER_BLOCK_32] |= (unsigned int)1 << (el->pos % BITS_PER_BLOCK_32);
							  ++el->pos;
						  }
						  else
						  // the <pos>-th bit must be set to 0
						  {
							  // due to initialization, the <pos>-th bit is already set to 0
							  // => increment the position to be examined
							  ++el->pos;
						  }
					  }
					  while (stack != NULL);
					}  
				}
			}
			// get next input gene combination vector
			while(nextCombination(comb));
			
			freeInputCombination(comb);

		}
	}
}

/**
 * R Wrapper for bestFitExtension() and reveal()
 * <measurements> contains an array of <number of genes> * <number of states> binary
 * values, where <number of genes> consecutive values form one state.
 * <numberOfStates> specifies the number of states in the array.
 * <maxK> is the highest number of input genes for a function.
 * If <method> is 0, bestFitExtension is called. If <method> is 1, reveal() is called.
 * Returns a list of lists. The outer list has <number of genes> elements, and the inner
 * list has one element for each equally-rated function for the current gene. Each of these
 * elements consists of a vector of input genes, the function as a binary vector, and
 * the error this function makes on the input time series.
 */
SEXP reconstructNetwork_R(SEXP inputStates, 
                          SEXP outputStates, 
                          SEXP perturbations,
                          SEXP numberOfStates, 
                          SEXP requiredDependencies,
                          SEXP excludedDependencies,                          
                          SEXP maxK, 
                          SEXP method, 
                          SEXP allSolutions,
                          SEXP returnPBN)
{
	int * _inputStates = INTEGER(inputStates);
	int * _outputStates = INTEGER(outputStates);
	
  unsigned int * _perturbations = (Rf_isNull(perturbations)? NULL : 
                          (unsigned int *)INTEGER(perturbations));
  
  int * _requiredDependencies = (Rf_isNull(requiredDependencies)? NULL : 
                                 INTEGER(requiredDependencies));
                                 
  int * _excludedDependencies = (Rf_isNull(excludedDependencies)? NULL : 
                                 INTEGER(excludedDependencies)); 
                                 
	int _numberOfStates = *INTEGER(numberOfStates);
	int _maxK = *INTEGER(maxK);
	int _method = *INTEGER(method);
	int _allSolutions = *INTEGER(allSolutions);
	int _returnPBN = *INTEGER(returnPBN);	

	unsigned int numGenes = length(inputStates) /  _numberOfStates;

	unsigned int numElts;
	if (numGenes % BITS_PER_BLOCK_32 == 0)
		numElts = numGenes / BITS_PER_BLOCK_32;
	else
		numElts = numGenes / BITS_PER_BLOCK_32 + 1;

	unsigned int encodedInputStates[numElts * _numberOfStates];
	memset(encodedInputStates,0,numElts * _numberOfStates * sizeof(unsigned int));

	unsigned int encodedOutputStates[numElts * _numberOfStates];
	memset(encodedOutputStates,0,numElts * _numberOfStates * sizeof(unsigned int));

	unsigned int state,gene;
	for (state = 0; state < _numberOfStates; ++state)
	// binary encoding of input states
	{
		for (gene = 0; gene < numGenes; ++gene)
			encodedInputStates[(numElts * state) + gene / BITS_PER_BLOCK_32] |=
				(_inputStates[state*numGenes + gene] << (gene % BITS_PER_BLOCK_32));
	}

	for (state = 0; state < _numberOfStates; ++state)
	// binary encoding of output states
	{
		for (gene = 0; gene < numGenes; ++gene)
			encodedOutputStates[(numElts * state) + gene / BITS_PER_BLOCK_32] |=
				(_outputStates[state*numGenes + gene] << (gene % BITS_PER_BLOCK_32));
	}

	FunctionListElement ** res = CALLOC(numGenes,sizeof(FunctionListElement *));
	unsigned int * errors = CALLOC(numGenes,sizeof(unsigned int));

	if (_method == 0)
		// perform Lähdesmäki's best-fit extension
		bestFitExtension(encodedInputStates,encodedOutputStates,_perturbations,
						 _numberOfStates,numGenes,
						 _requiredDependencies,_excludedDependencies,
						 _maxK,res,errors,_allSolutions,_returnPBN);
	else
		// start REVEAL algorithm
		reveal(encodedInputStates,encodedOutputStates,_perturbations,
								 _numberOfStates,numGenes,
								 _requiredDependencies,_excludedDependencies,
								 _maxK,res,errors,_allSolutions,_returnPBN);

//	for (gene = 0; gene < numGenes; ++gene)
//	{
//		printf("Gene %d\n",gene+1);
//		FunctionListElement * cur = res[gene];
//		while (cur != NULL)
//		{
//			printf("Input: ");
//			unsigned int k;
//			for (k = 0; k < cur->k; ++k)
//				printf("%d ",cur->inputGenes[k] + 1);
//			printf("\nFunction: ");
//			for (k = 0; k < (1 << cur->k); ++k)
//				printf("%d",GET_BIT(cur->transitionFunction[k / BITS_PER_BLOCK_32],k%BITS_PER_BLOCK_32));
//			printf("\nError: %d\n\n",errors[gene]);
//			cur = cur->next;
//		}
//	}

	// pack R objects
	SEXP resSXP;
	PROTECT(resSXP = allocList(numGenes));
	SEXP listPos1 = resSXP;

	for(gene = 0; gene < numGenes; ++gene)
	{
		SEXP geneSXP;
		unsigned int listLength = 0;
		FunctionListElement * cur = res[gene];

		while (cur != NULL)
		{
			++listLength;
			cur = cur->next;
		}

		PROTECT(geneSXP = allocList(listLength));

		cur = res[gene];

		unsigned int i;
		SEXP listPos2 = geneSXP;
		for (i = 0; i < listLength; ++i)
		{
			SEXP entrySXP,inputSXP,funcSXP,errorSXP;
			PROTECT(entrySXP = allocList(3));
			SET_TAG(entrySXP, install("input"));
			SET_TAG(CDR(entrySXP), install("func"));
			SET_TAG(CDR(CDR(entrySXP)), install("error"));

			PROTECT(inputSXP = allocVector(INTSXP,cur->k));
			int * array = INTEGER(inputSXP);
			unsigned int j;
			unsigned int numBits;
			if (cur->k == 1 && cur->inputGenes[0] == -1)
			{
				numBits = 1;
				array[0] = 0;
			}
			else
			{
				numBits = (unsigned int)1 << cur->k;
				for (j = 0; j < cur->k; j++)
					array[j] = cur->inputGenes[cur->k - j - 1] + 1;
			}
			SETCAR(entrySXP,inputSXP);
			UNPROTECT(1);

			PROTECT(funcSXP = allocVector(INTSXP,numBits));
			array = INTEGER(funcSXP);

      if (_returnPBN)
      // for PBN, functions are internally represented as bit vectors to save memory
  			dec2binC(array,(int*)cur->transitionFunction,(int*)&numBits);
  	  else
  	  // for incomplete truth tables, we need a third "don't care" value 
  	  // hence, each bit is an integer in this case (but far fewer entries => less memory)
	  		memcpy(array,(int*)cur->transitionFunction,numBits*sizeof(int));
			
			SETCADR(entrySXP,funcSXP);
			UNPROTECT(1);

			PROTECT(errorSXP = allocVector(INTSXP,1));
			array = INTEGER(errorSXP);
			*array = errors[gene];
			SETCADDR(entrySXP,errorSXP);
			UNPROTECT(1);

			SETCAR(listPos2,entrySXP);
			UNPROTECT(1);

			listPos2 = CDR(listPos2);
			cur = cur->next;
		}
		UNPROTECT(1);
		SETCAR(listPos1,geneSXP);
		listPos1 = CDR(listPos1);
	}
	UNPROTECT(1);

	// free resources
	for (gene = 0; gene < numGenes; ++gene)
		freeFunctionList(&res[gene]);
	FREE(errors);
	FREE(res);

	return resSXP;
}
