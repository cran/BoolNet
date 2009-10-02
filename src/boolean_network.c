/**
 * C code to identify attractors in Boolean networks
 *
 * This is part of the BooleanNetwork R package.
 *
 * Copyright 2009 by Christoph Müssel and Zhou Dao
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
 * Identification of attractors
 */


/**
 * Internal structure describing a Boolean network
 */
typedef struct
{
	// the number of genes in the network
	unsigned int numGenes;

	// a vector specifying whether the genes are fixed:
	// -1 means the gene is not fixed, 1 and 0 means the
	// genes are fixed to the corresponding values
	int * fixedGenes;

	// an index array with the <i>-th entry
	// specifying the bit position of the <i>-th gene
	// in a state array - this is not always equal to <i>,
	// as fixed genes are not stored
	unsigned int * nonFixedGeneBits;

	// a vector encoding the input genes for all transition functions.
	int * inputGenes;

	// a vector of indices to split up <inputGenes> for the single
	// gene transition functions.
	int * inputGenePositions;

	// a vector encoding the return values of all transition functions
	int * transitionFunctions;

	// a vector of indices to split up <transitionFunctions> for the single
	// genes.
	int * transitionFunctionPositions;
} BooleanNetwork;

/**
 * Structure that internally describes an attractor
 */
typedef struct Attractor
{
	// an array of states the attractor consists of
	unsigned int *involvedStates;

	// the number of array elements for one entry of <involvedStates>
	// - i.e. for more than 32 genes, <numElementsPerEntry> successive
	// array elements form one entry.
	unsigned int numElementsPerEntry;

	// the number of states in <involvedStates>
	int length;

	// the number of states in the basin of attraction
	unsigned int basinSize;

	// list pointer to the next element
	struct Attractor *next;
} Attractor, *pAttractor;


/**
 * A structure holding all information
 * retrieved by the algorithms in this file
 */
typedef struct
{
	// the number of elements in the following three arrays
	unsigned long long tableSize;

	// the states before the transition - can be NULL
	unsigned int *initialStates;


	// the resulting states of the transition table
	unsigned int *table;

	// the number of array elements for one entry of the table
	// - i.e. for more than 32 genes, <numElementsPerEntry> successive
	// array elements in <table> and <originalStates> form one table entry.
	unsigned int numElementsPerEntry;

	// the attractors the corresponding states belong to
	unsigned int *attractorAssignment;

	// the number of transitions needed to go from the original state
	// (before transition, not stored here as it is defined by the order)
	// to the attractor
	unsigned int *stepsToAttractor;

	// the list of attractors
	pAttractor attractorList;
} AttractorInfo, *pAttractorInfo;

/**
 * A structure saving search trees for states.
 * State trees are employed if a non-exhaustive search is
 * conducted to provide an efficient method to lookup
 * states that have already been visited.
 */
typedef struct STN
{
	// The left child node (containing "smaller" states)
	struct STN * leftChild;

	// The right child node (containing "larger" states)
	struct STN * rightChild;

	// The node holding the next state after a transition has
	// been applied to the current state
	struct STN * successor;

	// The state itself - a set of binary-encoded integers with
	// <ceil(number of Genes / 32)> elements
	unsigned int * data;

	// The basin of attraction to which the state belongs
	unsigned int attractorAssignment;

	// The number of transitions required to enter the attractor
	unsigned int stepsToAttractor;

} StateTreeNode;


void stateTransition(unsigned int * currentState, unsigned int * nextState, BooleanNetwork * net);

unsigned int nodeCount = 0;


/**
 * Allocate a new AttractorInfo structure for <tableSize> states
 */
pAttractorInfo allocAttractorInfo(unsigned long long tableSize, unsigned int numGenes)
{
	pAttractorInfo res = (pAttractorInfo)calloc(1,sizeof(AttractorInfo));
	if ((numGenes % BITS_PER_BLOCK_32) == 0)
		res->numElementsPerEntry = numGenes/BITS_PER_BLOCK_32;
	else
		res->numElementsPerEntry = numGenes/BITS_PER_BLOCK_32 + 1;
	res->tableSize = tableSize;
	res->initialStates = 0;
	res->table = (unsigned int*) calloc(tableSize * res->numElementsPerEntry,sizeof(unsigned int));
	res->attractorAssignment = (unsigned int*) calloc(tableSize,sizeof(unsigned int));
	res->stepsToAttractor = (unsigned int*) calloc(tableSize,sizeof(unsigned int));
	return res;
}

/**
 * Free a list of attractor structures
 */
void freeAttractorList(pAttractor p)
{
	do
	{
		pAttractor next = p->next;
		free(p->involvedStates);
		free(p);
		p = next;
	}
	while(p != NULL);
}

/**
 * Frees an AttractorInfo structure including
 * all sub-elements and the attractor list
 */
void freeAttractorInfo(pAttractorInfo p)
{
	if (p->initialStates != 0)
		free(p->initialStates);
	free(p->table);
	free(p->attractorAssignment);
	free(p->stepsToAttractor);
	freeAttractorList(p->attractorList);
	free(p);
}

/**
 * Allocate a new node of a state tree.
 * <leftChild> and <rightChild> are pointers to the
 * left and right subtrees.
 * <successor> is the state reached after a state transition.
 * <data> is an array of binary-encoded integers of length <numElements>
 * describing the state.
 * <attractorAssignment> is the basin of attraction the state belongs to,
 * and <stepsToAttractor> is the number of transitions required to enter
 * the attractor.
 *
 * Returns a state tree node with the supplied values.
 */
StateTreeNode * allocTreeNode(StateTreeNode * leftChild,
							  StateTreeNode * rightChild,
							  StateTreeNode * successor,
							  unsigned int * data,
							  unsigned int numElements,
							  unsigned int attractorAssignment,
							  unsigned int stepsToAttractor)
{
	StateTreeNode * res = calloc(1,sizeof(StateTreeNode));
	res->leftChild = leftChild;
	res->rightChild = rightChild;
	res->successor = successor;
	res->data = calloc(numElements,sizeof(unsigned int));
	memcpy(res->data,data,numElements*sizeof(unsigned int));
	res->attractorAssignment = attractorAssignment;
	res->stepsToAttractor = stepsToAttractor;
	++nodeCount;
	return res;
}

/**
 * Recursive helper function for findNode()
 */
StateTreeNode * findNodeRec(StateTreeNode * parent, unsigned int * data, unsigned int numElements)
{
	unsigned int direction = 0;
	int i;
	for (i = numElements - 1; i >= 0; --i)
	{
		if (data[i] > parent->data[i])
			direction = 1;
		if (data[i] < parent->data[i])
			direction = 2;
	}
	switch(direction)
	{
		case 0:
			return parent;
		case 1:
			if (parent->rightChild == 0)
			{
				parent->rightChild =  allocTreeNode(0,0,0,data,numElements,0,0);
				return parent->rightChild;
			}
			else
				return findNodeRec(parent->rightChild,data,numElements);
		case 2:
			if (parent->leftChild == 0)
			{
				parent->leftChild =  allocTreeNode(0,0,0,data,numElements,0,0);
				return parent->leftChild;
			}
			else
				return findNodeRec(parent->leftChild,data,numElements);
	}
	// should never be reached
	return 0;
}

/**
 * Recursively find the node corresponding to state <data> in the state tree <root>,
 * or insert the node if it does not exist.
 * <numElements> is the size of the state vector <data>.
 *
 * Returns the (possibly newly created) node belonging to <data>. If the tree is empty,
 * <root> is set to this node.
 */
StateTreeNode * findNode(StateTreeNode ** root, unsigned int * data, unsigned int numElements)
{
	if (*root == 0)
	{
		*root = allocTreeNode(0,0,0,data,numElements,0,0);
		return *root;
	}
	return findNodeRec(*root,data,numElements);
}

/**
 * Returns the successor of the supplied state node <current>.
 * If the state transition has not yet been calculated,
 * a new node is inserted into the tree <root>.
 * <numElementsPerEntry> is the number of array elements required to store one state.
 * <net> describes the network for which a state transition is performed.
 */
StateTreeNode * findSuccessor(StateTreeNode ** root, StateTreeNode * current,
							 unsigned int numElementsPerEntry, BooleanNetwork * net)
{
	if (current->successor == 0)
	// the state does not exist => calculate state transition and insert it
	{
		unsigned int nextState[numElementsPerEntry];
		stateTransition(current->data,nextState,net);
		current->successor = findNode(root,nextState,numElementsPerEntry);
	}
	return current->successor;
}

/**
 * Traverse the tree supplied by <root> in-order, and write the values
 * of the tree nodes to the corresponding arrays <initialStates>,
 * <table>, <attractorAssignment>, and <stepsToAttractor>.
 * <numElements> is the number of array elements allocated by one state.
 * <nodeNo> is the current value of the node counter used for the array indices
 * and increased during recursion. This should be initially set to 0.
 */
void inOrderSerializeTree(StateTreeNode * root,
						  unsigned int * initialStates,
						  unsigned int * table,
						  unsigned int * attractorAssignment,
						  unsigned int * stepsToAttractor,
						  unsigned int numElements,
						  unsigned int * nodeNo)
{

	if (root->leftChild != 0)
	// recursive descent of left subtree
		inOrderSerializeTree(root->leftChild,initialStates,table,attractorAssignment,
							 stepsToAttractor,numElements,nodeNo);

	// write the state itself
	memcpy(&initialStates[numElements* (*nodeNo)],root->data,numElements*sizeof(unsigned int));
	memcpy(&table[numElements* (*nodeNo)],root->successor->data,numElements*sizeof(unsigned int));
	attractorAssignment[*nodeNo] = root->attractorAssignment;
	stepsToAttractor[*nodeNo] = root->stepsToAttractor;
	*nodeNo = *nodeNo + 1;

	if (root->rightChild != 0)
	// recursive descent of right subtree
		inOrderSerializeTree(root->rightChild,initialStates,table,attractorAssignment,
							 stepsToAttractor,numElements,nodeNo);
}

/**
 * Free a state tree supplied by <node> recursively
 */
void freeTreeNode(StateTreeNode * node)
{
	if (node->leftChild != 0)
		freeTreeNode(node->leftChild);
	if (node->rightChild != 0)
			freeTreeNode(node->rightChild);
	free(node->data);
	free(node);
	--nodeCount;
}

/**
 * Make a transition from <currentState> to the next state.
 * <currentState> is a binary-coded integer with <numberOfGenes> used bits.
 * <fixedGenes> is an array of values specifying whether gene <i> is fixed (0 or 1) or not (-1).
 * <inputGenes> provides the input genes for all transition functions and can be split up
 * for a single function according to <inputGenePositions>.
 * <transitionFunctions> provides the truth tables for all transition functions and can be split up
 * for a single function according to <transitionFunctionPositions>.
 *
 * The return value is the next state, encoded in a single integer.
 */
void stateTransition(unsigned int * currentState, unsigned int * nextState, BooleanNetwork * net)
{
	unsigned int i = 0, k = 0, idx = 0;

	unsigned int elementsPerEntry;

	if (net->numGenes % BITS_PER_BLOCK_32 == 0)
	{
		elementsPerEntry = net->numGenes / BITS_PER_BLOCK_32;
	}
	else
	{
		elementsPerEntry = net->numGenes / BITS_PER_BLOCK_32 + 1;
	}

	for (i = 0; i < elementsPerEntry; ++i)
		nextState[i] = 0;

	for (i = 1; i <= net->numGenes; ++i)
	{
		if (net->fixedGenes[i-1] == -1)
		// the gene is not fixed
		{
			unsigned long long inputdec = 0;

			for (k = net->inputGenePositions[i-1]; k < net->inputGenePositions[i]; k++)
			{
				if (net->inputGenes[k])
				// if the input of the function is not 0 (constant gene), take input bit
				{
					unsigned int gene = net->inputGenes[k] - 1;
					unsigned int bit;

					if (net->fixedGenes[gene] == -1)
						bit = (GET_BIT(currentState[net->nonFixedGeneBits[gene] / BITS_PER_BLOCK_32],
							   net->nonFixedGeneBits[gene] % BITS_PER_BLOCK_32));
					else
						// fixed genes are not encoded in the states
						// => take them from fixedGenes vector
						bit = net->fixedGenes[gene];
					inputdec |= bit	<< (net->inputGenePositions[i] - k - 1);
				}
			}
			// determine transition function
			int transition = net->transitionFunctions[net->transitionFunctionPositions[i-1] + inputdec];

			if(transition != -1)
				// apply transition function
				nextState[idx / BITS_PER_BLOCK_32] |= (transition << (idx % BITS_PER_BLOCK_32));
			else
				// this is a dummy function for a constant gene
				// => value does not change
				nextState[idx / BITS_PER_BLOCK_32] |= (GET_BIT(currentState[idx / BITS_PER_BLOCK_32],
															idx % BITS_PER_BLOCK_32) << (idx % BITS_PER_BLOCK_32));

													//(GET_BIT(currentState[(i-1) / BITS_PER_BLOCK_32],
				                                     //   		 (i-1) % BITS_PER_BLOCK_32) << (idx % BITS_PER_BLOCK_32));
			++idx;
		}
	}
}

/**
 * Retrieves the result column of the state transition table.
 * <numberOfGenes> specifies the total number of genes.
 * <fixedGenes> is an array of values specifying whether gene <i> is fixed (0 or 1) or not (-1).
 * <inputGenes> provides the input genes for all transition functions and can be split up
 * for a single function according to <inputGenePositions>.
 * <transitionFunctions> provides the truth tables for all transition functions and can be split up
 * for a single function according to <transitionFunctionPositions>.
 */
unsigned long long * getTransitionTable(BooleanNetwork * net)
{
	unsigned long long i = 0;

	// determine number of fixed genes
	int numFixed = 0;
	for( i = 0; i < net->numGenes; i++)
		if(net->fixedGenes[i] != -1)
			++numFixed;

	int numNonFixed = net->numGenes - numFixed;

	// allocate truth table with 2^(non-fixed genes) elements
	unsigned long long numberOfElements = pow(2,numNonFixed);
	unsigned long long * table = calloc(numberOfElements,sizeof(unsigned long long));
	if (table == 0)
	{
		Rf_error("Too few memory available!");
	}

	unsigned long long initialState = 0;

	// calculate state transitions
	for(initialState = 0; initialState < numberOfElements; ++initialState)
	{
		//state is simply the binary encoding of the counter
		//calculate transition
		table[initialState] = 0;
		stateTransition((unsigned int *)&initialState,
						(unsigned int *)&table[initialState],
						net);
	}
	return table;
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
 * Retrieves attractors from a given transition table <table> with <numberOfStates> entries.
 * <specialInitializations> is a vector with <numSpecialInitializations>*2 entries, where the elements
 * 2*k entries are the gene indices, and the 2*k+1 entries are the corresponding initialization values (0 or 1)
 *
 * Returns a list of attractors - the last element of this list is empty!
 */
pAttractorInfo getAttractors(unsigned long long * table, unsigned int numberOfStates, int *specialInitializations, int numSpecialInitializations,
							unsigned int numberOfGenes)
{
	unsigned long long i;
	unsigned int current_attractor = 0, elementsPerEntry;

	if (numberOfGenes <= 32)
	{
		elementsPerEntry = 1;
	}
	else
	{
		elementsPerEntry = 2;
	}

	pAttractorInfo result = allocAttractorInfo(numberOfStates,numberOfGenes);

	for (i = 0; i < numberOfStates; ++i)
	{
		memcpy(&result->table[i],&table[i],sizeof(unsigned int) * elementsPerEntry);
	}

	pAttractor attractorHead, attractorList,tmpList;
	attractorHead = attractorList = (pAttractor) calloc(1,sizeof(Attractor));
	attractorList->next = NULL;

	for(i = 0; i < numberOfStates; i++)
	{
		if(!result->attractorAssignment[i])
		// the current state has not yet been visited
		{

			// first attractor has number 1
			current_attractor++;

			unsigned long long begin = i;
			unsigned int steps = 0;

			while(!result->attractorAssignment[begin])
			// iterate while no attractor has been assigned
			{
				++steps;

				// first simply count steps until attractor is reached
				// - to get the distance to the attractor, this number is
				// later subtracted from the maximum distance
				result->stepsToAttractor[begin] = steps;
				result->attractorAssignment[begin] = current_attractor;

				// perform a state transition
				begin = table[begin];
			}
			if(current_attractor == result->attractorAssignment[begin])
			//calculate length and basin size of new attractor
			{
				attractorList->basinSize = steps;

				// fix the number of steps to the attractor by calculating
				// the maximum distance and subtracting the current value from it
				int maxstep = result->stepsToAttractor[begin];

				int rec = 0;
				unsigned long long tmp = i;
				while(tmp != begin)
				{
					rec++;
					result->stepsToAttractor[tmp] = maxstep - result->stepsToAttractor[tmp];
					tmp = table[tmp];
				}

				attractorList->length = steps - rec;

				attractorList->involvedStates = (unsigned int *) calloc(attractorList->length * elementsPerEntry,sizeof(unsigned int));
				attractorList->numElementsPerEntry = elementsPerEntry;

				int a=0;
				do
				{
					result->stepsToAttractor[tmp] = 0;
					//attractorList->involvedStates[a++] = tmp;
					memcpy(&attractorList->involvedStates[a],&tmp,elementsPerEntry*sizeof(unsigned int));
					tmp = table[tmp];
					a += elementsPerEntry;
				}
				while(tmp != begin); /* set steps of attractors to 0; add attractor stub information */

				//generate a next attractor space
				attractorList->next = (pAttractor)calloc(1,sizeof(Attractor));
				attractorList = attractorList->next;
				attractorList->next = NULL;
			}
			else
			//update existing attractor
			{
				// reset attractor number
				current_attractor--;

				// assign states to attractor basin, and
				// correct the numbers of steps to the attractor
				unsigned long long tmp = i;
				int maxstp = result->stepsToAttractor[begin] + steps;
				while(tmp != begin)
				{
					result->attractorAssignment[tmp] = result->attractorAssignment[begin];
					result->stepsToAttractor[tmp] = maxstp - result->stepsToAttractor[tmp] + 1;
					tmp = table[tmp];
				}

				// update basin size in attractor record
				tmpList = attractorHead;
				for (tmp = 1; tmp < result->attractorAssignment[begin]; tmp++)
					tmpList = tmpList->next;

				tmpList->basinSize = tmpList->basinSize + steps;
			}
		}
	}

	if(numSpecialInitializations>0)
	// if special initializations are supplied,
	// correct the basin size by omitting states that
	// are not valid according to these restrictions
	{
		unsigned int newBasin[current_attractor];
		memset(newBasin,0,sizeof(int)*current_attractor);
		for( i = 0; i < numberOfStates; i++)
		{
			int tmp = 1, j;
			for(j = 0; (j < numSpecialInitializations) & (tmp >0);j++)
				tmp = tmp & (((i& 1<<specialInitializations[2*j])>0?1:0 )== (specialInitializations[2*j+1]));

			if(tmp)
				newBasin[result->attractorAssignment[i]-1]++;
		}
		tmpList = attractorHead; i = 0;
		while(tmpList->next !=NULL)
		{
			tmpList->basinSize = newBasin[i++];
			tmpList = tmpList->next;
		}
	}

	result->attractorList = attractorHead;

	free(table);

	return result;
}

/**
 * Retrieves attractors only for a given set of input states supplied in <selectedStates>.
 * Here, <ceil(net->numGenes / 32)> consecutive array entries describe one state, thus
 * the array size is <ceil(net->numGenes / 32) * numberOfStates>
 * <net> describes the network structure.
 * <specialInitializations> is a vector with <numSpecialInitializations>*2 entries, where the elements
 * 2*k entries are the gene indices, and the 2*k+1 entries are the corresponding initialization values (0 or 1).
 */
pAttractorInfo getAttractorsForStates(unsigned int * selectedStates, unsigned int numberOfStates,
									  BooleanNetwork * net, int *specialInitializations, int numSpecialInitializations)
{
	unsigned long long i;
	unsigned int j, current_attractor = 0, elementsPerEntry;


	// calculate the number of array elements required for one state
	// (depending on the number of genes)
	if (net->numGenes % (sizeof(unsigned int) * 8) == 0)
	{
		elementsPerEntry = net->numGenes / BITS_PER_BLOCK_32;
	}
	else
	{
		elementsPerEntry = net->numGenes / BITS_PER_BLOCK_32 + 1;
	}

	// all states are stored in a tree for fast search
	StateTreeNode * stateTree = 0;

	pAttractor attractorHead, attractorList,tmpList;
	attractorHead = attractorList = (pAttractor) calloc(1,sizeof(Attractor));
	attractorList->next = NULL;

	for(i = 0; i < numberOfStates; i++)
	{
		// check whether the current state is already in the state tree, otherwise insert it
		StateTreeNode * currentState = findNode(&stateTree,&selectedStates[i*elementsPerEntry],elementsPerEntry);
		if(!currentState->attractorAssignment)
		// the current state has not yet been visited
		{

			// first attractor has number 1
			current_attractor++;

			unsigned int steps = 0;

			while(!currentState->attractorAssignment)
			// iterate while no attractor has been assigned
			{
				++steps;

				// first simply count steps until attractor is reached
				// - to get the distance to the attractor, this number is
				// later subtracted from the maximum distance
				currentState->stepsToAttractor = steps;
				currentState->attractorAssignment = current_attractor;

				// perform a state transition
				currentState = findSuccessor(&stateTree,currentState,elementsPerEntry,net);
			}
			if(current_attractor == currentState->attractorAssignment)
			//calculate length and basin size of new attractor
			{
				attractorList->basinSize = steps;

				// fix the number of steps to the attractor by calculating
				// the maximum distance and subtracting the current value from it
				int maxstep = currentState->stepsToAttractor;

				int rec = 0;
				StateTreeNode * tmp = findNode(&stateTree,&selectedStates[i*elementsPerEntry],elementsPerEntry);

				while(memcmp(tmp->data,currentState->data,elementsPerEntry*sizeof(unsigned int)))
				{
					rec++;
					tmp->stepsToAttractor = maxstep - tmp->stepsToAttractor;

					// perform a state transition
					tmp = findSuccessor(&stateTree,tmp,elementsPerEntry,net);
				}

				attractorList->length = steps - rec;

				attractorList->involvedStates = (unsigned int *) calloc(attractorList->length * elementsPerEntry,sizeof(unsigned int));
				attractorList->numElementsPerEntry = elementsPerEntry;

				int a=0;
				do
				{
					tmp->stepsToAttractor = 0;
					memcpy(&attractorList->involvedStates[a],tmp->data,elementsPerEntry*sizeof(unsigned int));
					tmp = findSuccessor(&stateTree,tmp,elementsPerEntry,net);
					a += elementsPerEntry;
				}
				while(memcmp(tmp->data,currentState->data,elementsPerEntry*sizeof(unsigned int)));

				//generate a next attractor space
				attractorList->next = (pAttractor)calloc(1,sizeof(Attractor));
				attractorList = attractorList->next;
				attractorList->next = NULL;
			}
			else
			//update existing attractor
			{
				// reset attractor number
				current_attractor--;

				// assign states to attractor basin, and
				// correct the numbers of steps to the attractor
				StateTreeNode * tmp = findNode(&stateTree,&selectedStates[i*elementsPerEntry],elementsPerEntry);
				int maxstp = currentState->stepsToAttractor + steps;

				while(memcmp(tmp->data,currentState->data,elementsPerEntry*sizeof(unsigned int)))
				{
					tmp->attractorAssignment = currentState->attractorAssignment;
					tmp->stepsToAttractor = maxstp - tmp->stepsToAttractor + 1;

					//perform a state transition
					tmp = findSuccessor(&stateTree,tmp,elementsPerEntry,net);
				}

				// update basin size in attractor record
				tmpList = attractorHead;

				unsigned int i;
				for (i = 1; i < currentState->attractorAssignment; ++i)
					tmpList = tmpList->next;

				tmpList->basinSize = tmpList->basinSize + steps;
			}
		}
	}

	if(numSpecialInitializations > 0)
	// if special initializations are supplied,
	// correct the basin size by omitting states that
	// are not valid according to these restrictions
	{
		unsigned int newBasin[current_attractor];
		memset(newBasin,0,sizeof(int)*current_attractor);

		for( i = 0; i < numberOfStates; ++i)
		{
			StateTreeNode * currentState = findNode(&stateTree,&selectedStates[i*elementsPerEntry],elementsPerEntry);
			int fitsInitializations = 1;
			for(j = 0; (j < numSpecialInitializations) & (fitsInitializations != 0); ++j)
			{
				unsigned int genePosition = specialInitializations[2*j];
				unsigned int geneValue = specialInitializations[2*j+1];
				fitsInitializations = fitsInitializations &
										((currentState->data[genePosition / BITS_PER_BLOCK_32] &
										 (1 << (geneValue % BITS_PER_BLOCK_32))) == geneValue);
			}

			if(fitsInitializations)
				++newBasin[currentState->attractorAssignment-1];
		}
		tmpList = attractorHead;
		i = 0;
		while(tmpList->next !=NULL)
		{
			tmpList->basinSize = newBasin[i++];
			tmpList = tmpList->next;
		}
	}

	pAttractorInfo result = allocAttractorInfo(nodeCount,net->numGenes);
	result->attractorList = attractorHead;
	result->initialStates = calloc(result->tableSize * elementsPerEntry,sizeof(unsigned int));

	unsigned int nodeNo = 0;

	// build a series of arrays by in-order traversing the tree
	inOrderSerializeTree(stateTree,
						 result->initialStates,
						 result->table,
						 result->attractorAssignment,
						 result->stepsToAttractor,
						 elementsPerEntry,
						 &nodeNo);
	freeTreeNode(stateTree);
	return result;
}

/**
 * The R wrapper function for getAttractors.
 * Arguments:
 * inputGenes			An integer vector containing the concatenated input gene lists
 * 						of *all* transition functions
 * inputGenePositions	A vector of positions to split up <inputGenes> for each transition function
 * transitionFunctions	The concatenated truth table result columns of the transition functions of all genes.
 * transitionFunctionPositions	A vector of positions to split up <transitionFunctions> for each transition function.
 * fixedGenes			A vector that contains -1 for all genes that can be both ON and OFF, 1 for genes that are always ON,
 * 						and 0 for genes that are always OFF.
 * specialInitializations	A matrix of special initializations supplied by the user. The first row contains the genes,
 * 							the second row contains the corresponding initialization values.
 */
SEXP getAttractors_R(SEXP inputGenes,
					 SEXP inputGenePositions,
					 SEXP transitionFunctions,
					 SEXP transitionFunctionPositions,
					 SEXP fixedGenes,
					 SEXP specialInitializations,
					 SEXP startStates)
{

	// decode information in SEXP for use in C

	BooleanNetwork network;
	network.numGenes = length(fixedGenes);

	network.inputGenes = INTEGER(inputGenes);
	network.inputGenePositions = INTEGER(inputGenePositions);
	network.transitionFunctions = INTEGER(transitionFunctions);
	network.transitionFunctionPositions = INTEGER(transitionFunctionPositions);
	network.fixedGenes = INTEGER(fixedGenes);
	network.nonFixedGeneBits = calloc(network.numGenes,sizeof(unsigned int));

	int* _specialInitializations = NULL;
	if (!isNull(specialInitializations))
		_specialInitializations = INTEGER(specialInitializations);

	// count fixed genes, and create an index array for non-fixed genes:
	// <network.nonFixedGeneBits[i]> contains the bit positions in a state
	// at which the <i>-th gene is stored - this is different from <i>
	// as fixed genes are not stored
	unsigned int numNonFixed = 0, i;

	for(i = 0; i < network.numGenes; i++)
		if(network.fixedGenes[i] == -1)
		{
			network.nonFixedGeneBits[i] = numNonFixed++;
		}

	pAttractorInfo res;

	if (isNull(startStates) || length(startStates) == 0)
	// no start states supplied => perform exhaustive search
	{
		// calculate transition table
		unsigned long long * table = getTransitionTable(&network);

		if (table == 0)
			return R_NilValue;

		unsigned long long numStates = pow(2,numNonFixed);

		// find attractors
		res = getAttractors(table, numStates, _specialInitializations, length(specialInitializations), network.numGenes);
	}
	else
	// start states supplied => only identify attractors for these states
	{
		unsigned int* _startStates = (unsigned int*) INTEGER(startStates);
		res = getAttractorsForStates(_startStates, length(startStates),
									  &network,_specialInitializations, length(specialInitializations));
	}

	// pack results in SEXP structure for return value
	SEXP resSXP;
	SEXP stateInfoSXP;

	// create a result list with two elements (attractors and transition table)
	PROTECT(resSXP = allocList(2));
	SET_TAG(resSXP, install("stateInfo"));
	SET_TAG(CDR(resSXP), install("attractors"));

	// create a 3-element list for the transition table
	PROTECT(stateInfoSXP = allocList(4));
	SET_TAG(stateInfoSXP, install("table"));
	SET_TAG(CDR(stateInfoSXP), install("attractorAssignment"));
	SET_TAG(CDR(CDR(stateInfoSXP)), install("stepsToAttractor"));
	SET_TAG(CDR(CDR(CDR(stateInfoSXP))), install("initialStates"));

	// write transition table result column
	SEXP tableSXP;
	PROTECT(tableSXP = allocVector(INTSXP,res->tableSize * res->numElementsPerEntry));
	int* array = INTEGER(tableSXP);
	for (i = 0; i < res->tableSize; ++i)
	{
		// the transition table results do not contain fixed genes => insert them
		insertFixedGenes(&res->table[i*res->numElementsPerEntry],network.fixedGenes,network.numGenes);
		memcpy(&array[i*res->numElementsPerEntry],&res->table[i*res->numElementsPerEntry],
				res->numElementsPerEntry * sizeof(unsigned int));
	}
	SETCAR(stateInfoSXP,tableSXP);
	UNPROTECT(1);

	// write attractor assignment vector for states in transition table
	SEXP assignmentSXP;
	PROTECT(assignmentSXP = allocVector(INTSXP,res->tableSize));
	array = INTEGER(assignmentSXP);
	memcpy(array,res->attractorAssignment,res->tableSize * sizeof(int));
	SETCADR(stateInfoSXP,assignmentSXP);
	UNPROTECT(1);

	// write a vector with number of transitions from a state to an attractor
	SEXP stepSXP;
	PROTECT(stepSXP = allocVector(INTSXP,res->tableSize));
	array = INTEGER(stepSXP);
	memcpy(array,res->stepsToAttractor,res->tableSize * sizeof(int));
	SETCADDR(stateInfoSXP,stepSXP);
	UNPROTECT(1);

	// if available, write the original states
	SEXP initialStateSXP;
	if (res->initialStates == 0)
		initialStateSXP = R_NilValue;
	else
	// if start states are specified, the initial states for each transition have to be saved as well,
	// as they cannot be inferred by enumeration
	{
		PROTECT(initialStateSXP = allocVector(INTSXP,res->tableSize * res->numElementsPerEntry));
		array = INTEGER(initialStateSXP);
		for (i = 0; i < res->tableSize; ++i)
		{
			// the transition table results do not contain fixed genes => insert them
			insertFixedGenes(&res->initialStates[i*res->numElementsPerEntry],network.fixedGenes,network.numGenes);
			memcpy(&array[i*res->numElementsPerEntry],
					&res->initialStates[i*res->numElementsPerEntry],res->numElementsPerEntry * sizeof(unsigned int));
		}
		SETCADDDR(stateInfoSXP,initialStateSXP);
		UNPROTECT(1);
	}

	// assign to result list
	SETCAR(resSXP,stateInfoSXP);
	UNPROTECT(1);

	// count attractors
	unsigned int numAttractors = 0;
	pAttractor el;

	for(el = res->attractorList; el->next != NULL; el = el->next)
		++numAttractors;

	// write attractors
	SEXP attractorsSXP;
	PROTECT(attractorsSXP = allocList(numAttractors));
	SEXP listPos = attractorsSXP;
	for(el = res->attractorList, i=0; el->next != NULL; el = el->next, ++i)
	{
		SEXP attractorSXP;
		// each attractor is a 2-element list with a list of states included
		// in the attractor and the size of the basin
		PROTECT(attractorSXP = allocList(2));
		SET_TAG(attractorSXP, install("involvedStates"));
		SET_TAG(CDR(attractorSXP), install("basinSize"));

		SEXP stateSXP;
		PROTECT(stateSXP = allocVector(INTSXP,el->length * el->numElementsPerEntry));
		array = INTEGER(stateSXP);
		for (i = 0; i < el->length; ++i)
		{
				// insert fixed gene values, as they are missing in the calculated results
				insertFixedGenes(&el->involvedStates[i*el->numElementsPerEntry],network.fixedGenes,network.numGenes);
				memcpy(&array[i*el->numElementsPerEntry],
					  &el->involvedStates[i*el->numElementsPerEntry],el->numElementsPerEntry*sizeof(unsigned int));
		}
		SETCAR(attractorSXP,stateSXP);

		// write basin size
		SEXP basinSXP;
		PROTECT(basinSXP = allocVector(INTSXP,1));
		array = INTEGER(basinSXP);
		array[0] = el->basinSize;
		SETCADR(attractorSXP,basinSXP);
		SETCAR(listPos,attractorSXP);
		if (el->next != NULL)
			listPos = CDR(listPos);
		UNPROTECT(3);
	}
	UNPROTECT(1);
	SETCADR(resSXP,attractorsSXP);

	UNPROTECT(1);

	// free resources
	freeAttractorInfo(res);
	free(network.nonFixedGeneBits);

	return(resSXP);
}

/**
 * Reconstruction of Boolean networks
 */

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
	unsigned int * inputGenes;

	// a bit vector containing the transition function
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

/**
 * Add a new element to the function list specified by <root>.
 * Set the corresponding values <k>,<inputGenes> and <transitionFunction>
 * by copying the supplied values.
 * <transitionFunctionSize> is the number of elements needed for the transition function vector.
 * Returns the new list element.
 */
inline FunctionListElement * addFunctionListElement(FunctionListElement ** root,
											 unsigned int k,
											 unsigned int transitionFunctionSize,
											 unsigned int * inputGenes,
											 unsigned int * transitionFunction)
{
	FunctionListElement * el = calloc(1,sizeof(FunctionListElement));
	el->k = k;

	el->inputGenes = calloc(k,sizeof(unsigned int));
	memcpy(el->inputGenes,inputGenes,sizeof(unsigned int) * k);

	el->transitionFunction = calloc(transitionFunctionSize,sizeof(unsigned int));
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
		free(current->inputGenes);
		free(current->transitionFunction);
		free(current);
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
inline FunctionStackElement * pushFunctionStackElement(FunctionStackElement ** stack,
											 unsigned int * transitionFunction,
											 unsigned int transitionFunctionSize,
											 unsigned int pos)
{
	FunctionStackElement * el = calloc(1,sizeof(FunctionStackElement));

	el->pos = pos;

	el->transitionFunction = calloc(transitionFunctionSize,sizeof(unsigned int));
	memcpy(el->transitionFunction,transitionFunction,sizeof(unsigned int) * transitionFunctionSize);

	el->next = *stack;
	*stack = el;
	return el;
}

/**
 * Remove the top-level element from <stack>.
 */
inline void deleteFunctionStackElement(FunctionStackElement ** stack)
{
	FunctionStackElement * el = *stack;
	*stack = (*stack)->next;
	free(el->transitionFunction);
	free(el);
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
inline bool nextCombination(unsigned int * comb, unsigned int * pos, unsigned int k, unsigned int n)
{
	bool posChanged = false;

	// find the first position that has not been set
	// to its maximum number
	while(comb[*pos] == n - *pos - 1 && *pos < k)
	{
		++ *pos;
		posChanged = true;
	}
	if (*pos == k)
		// all elements have been listed
		return false;

	if (posChanged)
	// reset lower-order positions, and increase
	// the position found previously
	{
		unsigned int i;
		++comb[*pos];
		for (i = *pos; i > 0; --i)
			comb[i-1] = comb[i] + 1;
		*pos = 0;
	}
	else
	// the current position has not been set to its maximum value
	// => increase it
		++comb[*pos];
	return true;
}


/**
 * Lähdesmäki's best-fit extension algorithm to infer Boolean networks
 * from time series.
 * <states> is an array of time series. Here, <ceil(numGenes / 32)>
 * consecutive array entries describe one state, thus
 * the array size is <ceil(numGenes / 32) * numberOfStates>.
 * <numGenes> is the number of genes in the time series.
 * <maxK> is the maximum number of inputs a function can have.
 * <result> receives one list of possible functions for each gene, and
 * <errors> stores the error of these functions on the time series.
 */
void bestFitExtension(unsigned int * inputStates, unsigned int * outputStates,
					  unsigned int numStates, unsigned int numGenes, unsigned int maxK,
					  FunctionListElement ** result, unsigned int * errors)
{
	unsigned int i;
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
		// set error to maximum
		errors[i] = ~0;
		unsigned int k, s;

		// check for constant genes
		unsigned int geneVal = GET_BIT(outputStates[i/BITS_PER_BLOCK_32],i % BITS_PER_BLOCK_32);
		bool isConst = true;

		for (s = 1; s < numStates; ++s)
		{
			if (GET_BIT(outputStates[s*numElts + i/BITS_PER_BLOCK_32],i % BITS_PER_BLOCK_32)
				!= geneVal)
			{
				isConst = false;
				break;
			}
		}

		if (isConst)
		// gene is constant => simplest function already found!
		{
			unsigned int inputGenes = -1;
			addFunctionListElement(&result[i],1,1,&inputGenes,&geneVal);
			errors[i] = 0;
			bestLength[i] = 0;
		}

		for (k = 1; k <= maxK; ++k)
		// iterate over possible numbers of inputs
		{
			if (errors[i] == 0)
				break;

			// initialize gene combination vector
			unsigned int comb[k];
			unsigned int j;
			for (j = 0; j < k; ++j)
				comb[j] = k - j - 1;
			unsigned int pos = 0;

			// calculate the number of array elements needed
			// to represent the output of a function
			unsigned int array_sz = 1 << k;
			unsigned int numEltsFunc;
			if (array_sz % BITS_PER_BLOCK_32 == 0)
				numEltsFunc = array_sz / BITS_PER_BLOCK_32;
			else
				numEltsFunc = array_sz / BITS_PER_BLOCK_32 + 1;

			unsigned int c_0[array_sz];
			unsigned int c_1[array_sz];

			do
			{
				memset(c_0,0,sizeof(unsigned int) * array_sz);
				memset(c_1,0,sizeof(unsigned int) * array_sz);

				for (s = 0; s < numStates; ++s)
				// iterate over states and count errors
				{
					unsigned int input = 0, bit;
					for (bit = 0; bit <  k; ++bit)
					// build input by transferring the bits of the input genes
					// in the current state into a condensed input value
					{

						input |= (GET_BIT(inputStates[s * numElts
										   + comb[bit]/BITS_PER_BLOCK_32],comb[bit] % BITS_PER_BLOCK_32))
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

				if (error <= errors[i] && bestLength[i] >= k)
				{

					// create a stack of possible functions to be able
					// to split up if the 0-error and the 1-error at a certain position
					// are equal
					FunctionStackElement * stack =  NULL;

					// start with an element initialized to zero,
					// and first examine the 0-th position
					unsigned int f[numEltsFunc];
					memset(f,0,sizeof(unsigned int) * numEltsFunc);

					pushFunctionStackElement(&stack,f,numEltsFunc,0);
					do
					{
						// get top-level stack element
						FunctionStackElement * el = stack;

						if (el->pos == array_sz)
						// the top-level element on the stack is a complete function
						// => clean it up
						{
							// create a new element in the result function list containing the
							// completed function on the stack
							addFunctionListElement(&result[i],k,numEltsFunc,comb,el->transitionFunction);
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
							new->transitionFunction[el->pos / BITS_PER_BLOCK_32] |= 1 << (el->pos % BITS_PER_BLOCK_32);

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
							el->transitionFunction[el->pos / BITS_PER_BLOCK_32] |= 1 << (el->pos % BITS_PER_BLOCK_32);
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
			// get next input gene combination vector
			while(nextCombination(comb,&pos,k,numGenes));

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
 * The table counting the occurrences of gene value combinations is returned in <table>, which
 * must be initialized to a size of 2^<numChosenIndices>.
 * The return value is the entropy H(gene1,...,gene_n).
 */
inline double entropy(unsigned int * inputStates, unsigned int * outputStates,
				      unsigned int numStates, unsigned int elementsPerEntry,
				      unsigned int numGenes, unsigned int * chosenIndices,
				      unsigned int numChosenIndices,
				      unsigned int * table)
{
	unsigned int numEntries = 1 << numChosenIndices;

	// reset table to 0
	memset(table,0,sizeof(unsigned int) * numEntries);
	unsigned int state;

	for (state = 0; state < numStates; ++state)
	{
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
			result += ((double)table[tableIndex])/numStates*log2(((double)table[tableIndex])/numStates);
	}
	return -result;
}

/**
 * Liang's REVEAL algorithm for reconstruction of Boolean networks.
 * This version is enhanced for dealing with inconsistent samples by
 * taking the solutions that produce the minimum error.
 * <states> is an array of time series. Here, <ceil(numGenes / 32)>
 * consecutive array entries describe one state, thus
 * the array size is <ceil(numGenes / 32) * numberOfStates>.
 * <numGenes> is the number of genes in the time series.
 * <maxK> is the maximum number of inputs a function can have.
 * <result> receives one list of possible functions for each gene, and
 * <errors> stores the error of these functions on the time series.
 */
void reveal(unsigned int * inputStates, unsigned int * outputStates,
			unsigned int numStates, unsigned int numGenes, unsigned int maxK,
			FunctionListElement ** result, unsigned int * errors)
{
	unsigned int i;
	unsigned int numElts;

	// calculate block size in time series array
	if (numGenes % BITS_PER_BLOCK_32 == 0)
		numElts = numGenes / BITS_PER_BLOCK_32;
	else
		numElts = numGenes / BITS_PER_BLOCK_32 + 1;

	for (i = 0; i < numGenes; ++i)
	// iterate over genes
	{
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
			if (GET_BIT(outputStates[s*numElts + i/BITS_PER_BLOCK_32],i % BITS_PER_BLOCK_32)
				!= geneVal)
			{
				isConst = false;
				break;
			}
		}

		if (isConst)
		// gene is constant => simplest function already found!
		{
			unsigned int inputGenes = -1;
			addFunctionListElement(&result[i],1,1,&inputGenes,&geneVal);
			errors[i] = 0;
		}

		for (k = 1; k <= maxK; ++k)
		// iterate over possible numbers of inputs
		{
			if (errors[i] == 0)
				break;

			// initialize gene combination vector
			unsigned int comb[k];
			unsigned int j;
			for (j = 0; j < k; ++j)
				comb[j] = k - j - 1;
			unsigned int pos = 0;

			// calculate the number of array elements needed
			// to represent the output of a function
			unsigned int array_sz = 1 << k;
			unsigned int numEltsFunc;
			if (array_sz % BITS_PER_BLOCK_32 == 0)
				numEltsFunc = array_sz / BITS_PER_BLOCK_32;
			else
				numEltsFunc = array_sz / BITS_PER_BLOCK_32 + 1;


			do
			{

				// calculate entropy of input genes
				unsigned int table_input[array_sz];
				double entropy_input = entropy(inputStates, outputStates, numStates,
											   numElts, numGenes, comb, k, table_input);

				// calculate entropy of the combination of input genes and output value
				unsigned int comb_output[k+1];
				memcpy(comb_output,comb,sizeof(unsigned int) * k);
				comb_output[k] = numGenes + i;
				unsigned int table_output[1 << (k+1)];

				double entropy_all = entropy(inputStates, outputStates, numStates,
											 numElts, numGenes, comb_output, k+1, table_output);

				if (fabs(entropy_input - entropy_all) < 1E-6)
				// the two entropies are the same => output is fully explained by input
				{

					// a solution has been found => set errors to 0 so that
					// the algorithm does not search for higher k's
					errors[i] = 0;

					FunctionStackElement * stack =  NULL;

					// start with an element initialized to zero,
					// and first examine the 0-th position
					unsigned int f[numEltsFunc];
					memset(f,0,sizeof(unsigned int) * numEltsFunc);

					pushFunctionStackElement(&stack,f,numEltsFunc,0);
					do
					{
						// get top-level stack element
						FunctionStackElement * el = stack;

						if (el->pos == array_sz)
						// the top-level element on the stack is a complete function
						// => clean it up
						{
							// create a new element in the result function list containing the
							// completed function on the stack
							addFunctionListElement(&result[i],k,numEltsFunc,comb,el->transitionFunction);
							// remove the completed function from the stack
							deleteFunctionStackElement(&stack);
						}
						else
						if (table_input[el->pos] == 0 || table_output[el->pos | (1 << k)] == table_output[el->pos])
						// no information is available if the <pos>-th bit must be set to one or zero
						// => create two solution branches, one with the <pos>-th bit set to 1
						// and one with the <pos>-th bit set to 0
						{
							// create a new element on the stack with the <pos>-th bit set to 1
							FunctionStackElement * new = pushFunctionStackElement(&stack,el->transitionFunction,
																				  numEltsFunc,el->pos+1);
							new->transitionFunction[el->pos / BITS_PER_BLOCK_32] |= 1 << (el->pos % BITS_PER_BLOCK_32);

							// increment the position of the old element, which remains on the stack
							// with the <pos>-th bit set to 0 (due to initialization)
							++el->pos;
						}
						else
						if (table_output[el->pos | (1 << k)] > table_output[el->pos])
						// the <pos>-th bit must be set to 1
						{
							// set the <pos>-th bit of the top-level stack element to 1,
							// and increment the position to be examined
							el->transitionFunction[el->pos / BITS_PER_BLOCK_32] |= 1 << (el->pos % BITS_PER_BLOCK_32);
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
			// get next input gene combination vector
			while(nextCombination(comb,&pos,k,numGenes));

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
SEXP reconstructNetwork_R(SEXP inputStates, SEXP outputStates, SEXP numberOfStates, SEXP maxK, SEXP method)
{
	int * _inputStates = INTEGER(inputStates);
	int * _outputStates = INTEGER(outputStates);
	int _numberOfStates = *INTEGER(numberOfStates);
	int _maxK = *INTEGER(maxK);
	int _method = *INTEGER(method);

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

	FunctionListElement ** res = calloc(numGenes,sizeof(FunctionListElement *));
	unsigned int * errors = calloc(numGenes,sizeof(unsigned int));

	if (_method == 0)
		// perform Lähdesmäki's best-fit extension
		bestFitExtension(encodedInputStates,encodedOutputStates,
						 _numberOfStates,numGenes,_maxK,res,errors);
	else
		// start REVEAL algorithm
		reveal(encodedInputStates,encodedOutputStates,
								 _numberOfStates,numGenes,_maxK,res,errors);

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
			for (j = 0; j < cur->k; j++)
				array[j] = cur->inputGenes[cur->k - j - 1] + 1;
			SETCAR(entrySXP,inputSXP);
			UNPROTECT(1);

			unsigned int numBits = 1 << cur->k;
			PROTECT(funcSXP = allocVector(INTSXP,numBits));
			array = INTEGER(funcSXP);

			dec2bin(array,(int*)cur->transitionFunction,(int*)&numBits);
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
	free(errors);
	free(res);

	return resSXP;
}
