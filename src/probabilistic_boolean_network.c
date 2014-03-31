/**
 * C code for Markov simulations of Probabilistic Boolean Networks
 *
 * This is part of the BooleanNetwork R package.
 *
 * Copyright 2009/2010 by Christoph Müssel
 *
 * Contact christoph.muessel@uni-ulm.de
 *
 */
#include "uthash.h"
#include "common.h"

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define MATRIX_POOL_SIZE 1024

typedef struct
{
	int * inputGenes;
	int * transitionFunction;

	unsigned int numGenes;

	double probability;

	unsigned int functionIndex;
} PBNFunction;

/**
 * Internal structure describing a Boolean network
 */
typedef struct
{
	// the number of genes in the network
	unsigned int numGenes;

	// the number of non-fixed genes in the network
	unsigned int numNonFixedGenes;

	// a vector specifying whether the genes are fixed:
	// -1 means the gene is not fixed, 1 and 0 means the
	// genes are fixed to the corresponding values
	int * fixedGenes;

	// an index array with the <i>-th entry
	// specifying the bit position of the <i>-th gene
	// in a state array - this is not always equal to <i>,
	// as fixed genes are not stored
	unsigned int * nonFixedGeneBits;

	// an array containing an array of transition functions for each gene
	PBNFunction ** transitionFunctions;

	// the lengths of the arrays in <transitionFunctions>
	unsigned int * numFunctionsPerGene;

} ProbabilisticBooleanNetwork;

/**
 * An entry in the state probability matrix
 */
typedef struct
{
	// the initial state
	// (next state is the index of the matrix array)
	unsigned int state;

	// the probability for this state transition
	double probability;

	// used by the hash table
	UT_hash_handle hh;
} MatrixEntry;

/*
 * Structure that maintains a sparse matrix
 */
typedef struct
{
  MatrixEntry ** matrix;  
  
  ArrayListElement * entryPool;
  unsigned int poolArraySize;
  unsigned int currentEntry;  
  unsigned int numCols;
} SparseMatrix;

/*
 * Initialize a matrix with <numCols> columns
 * maintaining an entry pool with arrays of size
 * <poolArraySize>
 */
static inline SparseMatrix * allocSparseMatrix(unsigned int numCols, unsigned int poolArraySize)
{
  SparseMatrix * res = CALLOC(1, sizeof(SparseMatrix));
  res->matrix = CALLOC(numCols, sizeof(MatrixEntry *));
  memset(res->matrix,0,sizeof(MatrixEntry *) * numCols);
  res->poolArraySize = poolArraySize;
  res->currentEntry = 0;
  res->entryPool = NULL;
  res->numCols = numCols;
  return res;
}


/**
 * Calculate a transition table with one entry for each state. Each entry consists of
 * one bit per (non-fixed) function in <net>.
 * <tableSize> receives the number of elements in the table.
 * <numElements> receives the number of array elements occupied by one state.
 * Returns an array of states.
 */
unsigned int * probabilisticTransitionTable(ProbabilisticBooleanNetwork * net,
											unsigned int * tableSize, unsigned int * numElements)
{
	// determine number of fixed genes
	unsigned int totalFunctionCount = 0;
	unsigned int i, j, k, numNonFixed = 0;

	for (i = 0; i < net->numGenes; ++i)
	{
		if (net->fixedGenes[i] == -1)
		{
			totalFunctionCount += net->numFunctionsPerGene[i];
			++numNonFixed;
		}
	}

	if (totalFunctionCount % BITS_PER_BLOCK_32 == 0)
		*numElements = totalFunctionCount / BITS_PER_BLOCK_32;
	else
		*numElements = totalFunctionCount / BITS_PER_BLOCK_32 + 1;

	// allocate truth table with 2^(non-fixed genes) elements
	*tableSize = ((unsigned int)1 << numNonFixed);
	unsigned int * table = CALLOC(*tableSize * *numElements,sizeof(unsigned int));
	if (table == 0)
	{
		Rf_error("Too few memory available!");
	}

	unsigned int initialState = 0;

	// calculate state transitions
	for(initialState = 0; initialState < *tableSize; ++initialState)
	{
    R_CheckUserInterrupt();
		//state is simply the binary encoding of the counter

		//calculate transitions
		for (i = 0; i < net->numGenes; ++i)
		{
			if (net->fixedGenes[i] == -1)
			{
				for (j = 0; j < net->numFunctionsPerGene[i]; ++j)
				{
					PBNFunction * current = &net->transitionFunctions[i][j];
					unsigned int inputdec = 0;

					for (k = 0; k < current->numGenes; ++k)
					{
						if (current->inputGenes[k])
						// if the input of the function is not 0 (constant gene), take input bit
						{
							unsigned int gene = current->inputGenes[k] - 1;
							unsigned int bit;

							if (net->fixedGenes[gene] == -1)
								bit = (GET_BIT(initialState,
											   net->nonFixedGeneBits[gene]));
							else
								// fixed genes are not encoded in the states
								// => take them from fixedGenes vector
								bit = net->fixedGenes[gene];
							inputdec |= bit	<< (current->numGenes - k - 1);
						}
					}

					int transition = current->transitionFunction[inputdec];

					if(transition != -1)
						// apply transition function
						table[initialState * *numElements + current->functionIndex / BITS_PER_BLOCK_32]
								  |= (transition << (current->functionIndex % BITS_PER_BLOCK_32));
					else
						// this is a dummy function for a constant gene
						// => value does not change
						table[initialState * *numElements + current->functionIndex / BITS_PER_BLOCK_32]
						      |= (GET_BIT(initialState, net->nonFixedGeneBits[i]) << (current->functionIndex % BITS_PER_BLOCK_32));
				}
			}
		}

	}
	return table;
}

/**
 * Calculate combinations of function indices.
 * <maxVals> is an array containing the numbers of functions for each gene.
 * <combination> is the previous function combination from which the next
 * combination is calculated.
 * <size> is the number of functions.
 * Returns false if there are no more combinations, and true otherwise.
 */
static inline bool nextFunctionCombination(unsigned int * maxVals, unsigned int * combination, unsigned int size)
{
	unsigned int i;
	while(true)
	{
		for (i = 0; i < size; ++i)
		{
			++combination[i];
			if (combination[i] < maxVals[i])
				return true;
			else
			{
				if (i == size - 1)
					return false;
				else
					combination[i] = 0;
			}
		}
	}
}

/**
 * Add <value> to the element of <hash> whose start state is <state>,
 * or add the corresponding element to the hash table
 */
static inline void addToMatrixEntry(SparseMatrix * matrix, unsigned int column, unsigned int state, double value)
{
	MatrixEntry * entry;
	HASH_FIND_INT(matrix->matrix[column], &state, entry);

	if (entry == NULL)
	// add a new hash table entry
	{
	  if (matrix->currentEntry % matrix->poolArraySize == 0)
	    allocNewArray(&matrix->entryPool, matrix->poolArraySize, sizeof(MatrixEntry));
	
		entry = &(((MatrixEntry *)matrix->entryPool->array)[matrix->currentEntry % matrix->poolArraySize]);
		++matrix->currentEntry;
		
		entry->probability = 0.0;
		entry->state = state;
		HASH_ADD_INT(matrix->matrix[column],state,entry);
	}

	// increment probability
	entry->probability += value;
}

/**
 * Extract a state from the combined transition table of all functions <table>
 * depending on the function combination <combination>.
 * <numElements> is the number of array elements occupied by one state.
 * <initialState> is the state for which the successor state is looked up.
 * <net> holds the network information.
 * Returns the extracted state as an integer.
 */
static inline unsigned int extractState(unsigned int * table, unsigned int numElements,
						 unsigned int initialState, unsigned int * combination,
						 ProbabilisticBooleanNetwork * net)
{
	unsigned int i, j, res = 0;
	for (i = 0, j = 0; i < net->numGenes; ++i)
	{
		if (net->fixedGenes[i] == -1)
		// exclude fixed genes
		{
			unsigned int functionIndex = net->transitionFunctions[i][combination[j]].functionIndex;

			// insert the <functionIndex>th bit at the next position
			res |= (GET_BIT(table[initialState * numElements + functionIndex / BITS_PER_BLOCK_32],
							functionIndex % BITS_PER_BLOCK_32))
					<< net->nonFixedGeneBits[i];
			++j;
		}
	}
	return res;
}


/**
 * Free all hash tables in the matrix array <matrix>
 * and the array itself.
 * <size> is the number of columns of the matrix.
 */
void freeMatrix(SparseMatrix * matrix)
{
	
	unsigned int i;
	for (i = 0; i < matrix->numCols; ++i)
	{
		// free hash table
    HASH_CLEAR(hh,matrix->matrix[i]);
	}
	//
	FREE(matrix->matrix);
	freeArrayList(matrix->entryPool);

	// free array
	FREE(matrix);
}

/**
 * Start a Markov chain simulation on network <net> by performing <numIteration> matrix multiplications.
 * <outcomeSize> receives the number of elements in the returned probability vector.
 * <stateProbabilities> is the matrix of state transition probabilities with <outcomeSize> columns.
 * <startStates> is an optional vector with <numStartStates> encoded start states.
 */
double * markovSimulation(ProbabilisticBooleanNetwork * net, unsigned int numIterations,
						  unsigned int * outcomeSize, SparseMatrix ** stateProbabilities,
						  unsigned int * startStates, unsigned int numStartStates)
{
	unsigned int tableSize = 0, numElements = 0, i, k;

	// calculate combined transition table
	unsigned int * table = probabilisticTransitionTable(net,&tableSize,&numElements);

	SparseMatrix * matrix = allocSparseMatrix(tableSize, MATRIX_POOL_SIZE);
	*stateProbabilities = matrix;

	unsigned int combination[net->numNonFixedGenes];
	memset(combination,0,sizeof(unsigned int) * net->numNonFixedGenes);

	unsigned int maxCombination[net->numNonFixedGenes];
	unsigned int geneIndices[net->numNonFixedGenes];

	// write array with numbers of functions
	for (i = 0; i < net->numGenes; ++i)
	{
		if (net->fixedGenes[i] == -1)
		{
			maxCombination[net->nonFixedGeneBits[i]] = net->numFunctionsPerGene[i];
			geneIndices[net->nonFixedGeneBits[i]] = i;
		}
	}

	// set up transition probability matrix
	do
	// iterate over all possible function combinations
	{
		double probability = net->transitionFunctions[geneIndices[0]][combination[0]].probability;

		for (i = 1; i < net->numNonFixedGenes; ++i)
			probability *= net->transitionFunctions[geneIndices[i]][combination[i]].probability;

		for (i = 0; i < tableSize; ++i)
		// add probabilities for all reached states in the transition table
		{
      R_CheckUserInterrupt();
			unsigned int state = extractState(table,numElements,i,combination,net);
			addToMatrixEntry(matrix,state,i,probability);
		}
	}
	while (nextFunctionCombination(maxCombination,combination,net->numNonFixedGenes));

	FREE(table);
	
	double * outcome = CALLOC(tableSize,sizeof(double));
	double * oldOutcome = CALLOC(tableSize,sizeof(double));	
	*outcomeSize = tableSize;

	if (numStartStates == 0)
	// no start states => equal probability for all states
	{
		for (i = 0; i < tableSize; ++i)
		{
			outcome[i] = 1.0/tableSize;
		}
	}
	else
	// start states supplied => equal probability for all supplied states
	{
		for (i = 0; i < numStartStates; ++i)
		{
			outcome[startStates[i]] = 1.0/numStartStates;
		}
	}

	for (k = 0; k < numIterations; ++k)
	// perform <numIterations> matrix multiplications
	{
	  R_CheckUserInterrupt();
	  
		memcpy(oldOutcome,outcome,sizeof(double) * tableSize);

		for (i = 0; i < tableSize; ++i)
		{
			MatrixEntry * entry;
			outcome[i] = 0.0;

			for(entry = matrix->matrix[i]; entry != NULL; entry=entry->hh.next)
			{
				outcome[i] += oldOutcome[entry->state] * entry->probability;
			}

		}
	}
	
	FREE(oldOutcome);

	return outcome;
}

/**
 * R wrapper for the Markov simulation.
 * <inputGenes> is a concatenated integer vector of the input genes for all functions.
 * <inputGenePositions> is used to split up <inputGenes> into the single function inputs.
 * <transitionFunctions> contains the truth table result columns of all functions and is
 * split up according to <transitionFunctionPositions>.
 * <fixedGenes> specifies whether genes are fixed (0,1) or not (-1).
 * <functionAssignment> is a vector specifying which function belongs to which gene.
 * <probabilities> is a vector of probabilities for the functions.
 * <numSteps> is the number of matrix multiplications to be performed.
 * <startStates> is an optional vector of start states (otherwise all states are considered).
 * <cutoff> is the value below which values are considered to be 0.
 * <returnTable> specifies whether the result should contain a transition table
 * Returns a list containing a data frame of important states and their probabilities
 * and, optionally, the transition table.
 */
SEXP markovSimulation_R(SEXP inputGenes,
					 SEXP inputGenePositions,
					 SEXP transitionFunctions,
					 SEXP transitionFunctionPositions,
					 SEXP fixedGenes,
					 SEXP functionAssignment,
					 SEXP probabilities,
					 SEXP numSteps,
					 SEXP startStates,
					 SEXP cutoff,
					 SEXP returnTable)
{
	ProbabilisticBooleanNetwork network;

	int * _inputGenes = INTEGER(inputGenes);

	int * _inputGenePositions = INTEGER(inputGenePositions);

	int * _transitionFunctions = INTEGER(transitionFunctions);

	int * _transitionFunctionPositions = INTEGER(transitionFunctionPositions);

	int * _functionAssignment = INTEGER(functionAssignment);

	double * _probabilities = REAL(probabilities);

	int _numSteps = *INTEGER(numSteps);

	int * _startStates = NULL;
	if (!isNull(startStates))
		_startStates = INTEGER(startStates);

	double _cutoff = *REAL(cutoff);

	bool _returnTable = (bool)(*INTEGER(returnTable));

	unsigned int i, j;

	network.numGenes = length(fixedGenes);

	network.fixedGenes = CALLOC(network.numGenes,sizeof(int));
	memcpy(network.fixedGenes,INTEGER(fixedGenes),sizeof(unsigned int) * network.numGenes);
	network.nonFixedGeneBits = CALLOC(network.numGenes,sizeof(unsigned int));

	// determine which genes are not fixed
	network.numNonFixedGenes = 0;
	for(i = 0; i < network.numGenes; i++)
	{
			if(network.fixedGenes[i] == -1)
			{
				network.nonFixedGeneBits[i] = network.numNonFixedGenes++;
			}
	}

	network.numFunctionsPerGene = CALLOC(network.numGenes,sizeof(unsigned int));
	network.transitionFunctions = CALLOC(network.numGenes,sizeof(PBNFunction *));

	// count number of functions per gene
	for (i = 0; i < length(functionAssignment); ++i)
	{
		++network.numFunctionsPerGene[_functionAssignment[i]];
	}

	// allocate function vectors
	for (i = 0; i < network.numGenes; ++i)
	{
		network.transitionFunctions[i] = CALLOC(network.numFunctionsPerGene[i],sizeof(PBNFunction));
	}

	unsigned int functionCounter[network.numGenes];
	memset(functionCounter,0,sizeof(unsigned int) * network.numGenes);

	unsigned int functionIndex = 0;

	for (i = 0; i < length(functionAssignment); ++i)
	// decode functions
	{
		PBNFunction * current = &network.transitionFunctions[_functionAssignment[i]]
														   [functionCounter[_functionAssignment[i]]++];

		unsigned int numInputs = _inputGenePositions[i+1] - _inputGenePositions[i];

		current->inputGenes = CALLOC(numInputs,sizeof(int));
		memcpy(current->inputGenes,&_inputGenes[_inputGenePositions[i]],numInputs*sizeof(int));

		current->numGenes = numInputs;

		current->transitionFunction = CALLOC((unsigned int)1 << numInputs,sizeof(int));
		memcpy(current->transitionFunction,&_transitionFunctions[_transitionFunctionPositions[i]],
			   ((unsigned int)1 << numInputs)*sizeof(int));

		current->probability = _probabilities[i];

		if (network.fixedGenes[_functionAssignment[i]] == -1)
		{
			current->functionIndex = functionIndex++;
		}
		else
		{
			current->functionIndex = ~0;
		}
	}

	unsigned int startStateElements;
	if (network.numGenes % (sizeof(unsigned int) * 8) == 0)
	{
		startStateElements = network.numGenes / BITS_PER_BLOCK_32;
	}
	else
	{
		startStateElements = network.numGenes / BITS_PER_BLOCK_32 + 1;
	}

	// decode start states (remove fixed genes)
	unsigned int numStartStates;

	if (isNull(startStates))
		numStartStates = 0;
	else
		numStartStates = length(startStates)/startStateElements;

	unsigned int decodedStartStates[numStartStates];

	for (i = 0; i < numStartStates; ++i)
	{
		decodedStartStates[i] = 0;
		for (j = 0; j < network.numGenes; ++j)
		{
			if (network.fixedGenes[j] == -1)
			{
				decodedStartStates[i] |= GET_BIT(_startStates[i * startStateElements
				                                             + j / BITS_PER_BLOCK_32],
				                                 j % BITS_PER_BLOCK_32)
				                          << network.nonFixedGeneBits[j];
			}
		}
	}

	unsigned int outcomeSize = 0;

	SparseMatrix * matrix;

	// perform simulation
	double * outcome = markovSimulation(&network,_numSteps,&outcomeSize,&matrix,
										decodedStartStates,numStartStates);
										
	// determine number of non-zero results
	unsigned int nonZero = 0;
	for (i = 0; i < outcomeSize; ++i)
	{
		if (outcome[i] > _cutoff)
			++nonZero;
	}

	// encode results
	SEXP resSXP, stateSXP, probSXP;

	if (_returnTable)
		PROTECT(resSXP = allocList(5));
	else
		PROTECT(resSXP = allocList(2));

	PROTECT(stateSXP = allocVector(INTSXP,nonZero*startStateElements));
	PROTECT(probSXP = allocVector(REALSXP,nonZero));

	int * states = INTEGER(stateSXP);
	double * prob = REAL(probSXP);
  //return R_NilValue;

	for (i = 0, j = 0; i < outcomeSize; ++i)
	// encode states with non-zero probability
	{
		if (outcome[i] > _cutoff)
		{
			prob[j] = outcome[i];
			states[j*startStateElements] = i;

			// re-encode fixed genes
			insertFixedGenes((unsigned int *)&states[j*startStateElements],network.fixedGenes,network.numGenes);
			++j;
		}
	}

	SET_TAG(resSXP, install("states"));
	SET_TAG(CDR(resSXP), install("probabilities"));

	SETCAR(resSXP,stateSXP);
	SETCADR(resSXP,probSXP);

	if (_returnTable)
	// encode the transition table
	{
		SET_TAG(CDR(CDR(resSXP)), install("initialStates"));
		SET_TAG(CDR(CDR(CDR(resSXP))), install("nextStates"));
		SET_TAG(CDR(CDR(CDR(CDR(resSXP)))), install("transitionProbabilities"));

		//calculate size of transition table
		unsigned int tableSize = 0;

		for (i = 0; i < outcomeSize; ++i)
		{
			MatrixEntry * entry;

			for(entry = matrix->matrix[i]; entry != NULL; entry=entry->hh.next)
			{
				++tableSize;
			}
		}
		
		SEXP initialStateSXP, nextStateSXP, transitionProbSXP;
		PROTECT(initialStateSXP = allocVector(INTSXP,tableSize * startStateElements));
		PROTECT(nextStateSXP = allocVector(INTSXP,tableSize * startStateElements));
		PROTECT(transitionProbSXP = allocVector(REALSXP,tableSize));

		unsigned int * initialStates = (unsigned int *)INTEGER(initialStateSXP);
		unsigned int * nextStates = (unsigned int *)INTEGER(nextStateSXP);
		double * transitionProbs = (double *)REAL(transitionProbSXP);

		// encode non-zero states in transition matrix
		j = 0;
		for (i = 0; i < outcomeSize; ++i)
		{
			MatrixEntry * entry;

			for(entry = matrix->matrix[i]; entry != NULL; entry=entry->hh.next)
			{
				initialStates[j * startStateElements] = entry->state;
				insertFixedGenes(&initialStates[j*startStateElements],network.fixedGenes,network.numGenes);

				nextStates[j * startStateElements] = i;
				insertFixedGenes(&nextStates[j*startStateElements],network.fixedGenes,network.numGenes);

				transitionProbs[j] = entry->probability;
				++j;
			}
		}

		SETCADDR(resSXP,initialStateSXP);
		SETCADDDR(resSXP,nextStateSXP);
		SETCADDDR(CDR(resSXP),transitionProbSXP);
		UNPROTECT(3);
	}

	UNPROTECT(3);

	freeMatrix(matrix);
	FREE(outcome);

	return resSXP;
}
