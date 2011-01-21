/**
 * C code to identify attractors in Boolean networks
 *
 * This is part of the BooleanNetwork R package.
 *
 * Copyright 2009/2010 by Christoph Müssel and Zhou Dao
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
#include "random.h"
#include "common.h"

/**
 * Identification of attractors
 */

#define SYNC_MODE 0
#define ASYNC_MODE_RANDOM 1

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


typedef struct TTE
{
	unsigned int * initialState;
	unsigned int * nextState;
	struct TTE * next;
} TransitionTableEntry;

/**
 * Structure that internally describes an attractor
 */
typedef struct Attractor
{

	// an array of states the attractor consists of
	unsigned int *involvedStates;

	// if this is a complex attractor,
	// the transitions of the attractor are stored here
	TransitionTableEntry * table;

	// the number of elements in <table>
	unsigned int transitionTableSize;

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

	union
	{
		// In case of synchronous networks:
		struct
		{
			// The node holding the next state after a transition has
			// been applied to the current state
			struct STN * successor;

			// The basin of attraction to which the state belongs
			unsigned int attractorAssignment;

			// The number of transitions required to enter the attractor
			unsigned int stepsToAttractor;
		} sync;

		// In case of asynchronous networks:
		struct
		{
			// An array of all successor nodes
			struct STN ** successors;

			// The number of successor states
			unsigned int numSuccessors;

			// Dummy variable
			unsigned int unused;

		} async;
	} type;

	// The state itself - a set of binary-encoded integers with
	// <ceil(number of Genes / 32)> elements
	unsigned int * data;

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
	res->table = NULL;
	res->tableSize = tableSize;
	res->initialStates = NULL;
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
 * Free an AttractorInfo structure including
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
 * Insert a new transition from <initialState> to <nextState> into <table>.
 * <numElements> is the number of array elements occupied by a state.
 */
TransitionTableEntry * insertNewTransition(TransitionTableEntry ** table,
										   unsigned int * initialState,
										   unsigned int * nextState,
										   unsigned int numElements)
{
	TransitionTableEntry * entry = (TransitionTableEntry *)calloc(1,sizeof(TransitionTableEntry));
	entry->initialState = calloc(numElements,sizeof(unsigned int));
	entry->nextState = calloc(numElements,sizeof(unsigned int));
	memcpy(entry->initialState,initialState,numElements*sizeof(unsigned int));
	memcpy(entry->nextState,nextState,numElements*sizeof(unsigned int));
	entry->next = *table;
	*table = entry;
	return entry;
}

/**
 * Free a list-type transition table as used in complex attractors
 */
void freeTransitionTableEntry(TransitionTableEntry * t)
{
	do
	{
		TransitionTableEntry * next = t->next;
		free(t->initialState);
		free(t->nextState);
		free(t);
		t = next;
	}
	while (t != NULL);
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
	res->type.sync.successor = successor;
	res->data = calloc(numElements,sizeof(unsigned int));
	memcpy(res->data,data,numElements*sizeof(unsigned int));

	res->type.sync.attractorAssignment = attractorAssignment;
	res->type.sync.stepsToAttractor = stepsToAttractor;
	++nodeCount;
	return res;
}

/**
 * Recursive helper function for findNode()
 */
StateTreeNode * findNodeRec(StateTreeNode * parent, unsigned int * data, unsigned int numElements, bool * found)
{
	unsigned int direction = 0;
	int i;
	for (i = numElements - 1; i >= 0; --i)
	{
		if (data[i] > parent->data[i])
		{
			direction = 1;
			break;
		}
		if (data[i] < parent->data[i])
		{
			direction = 2;
		  break;
		}  
	}
	switch(direction)
	{
		case 0:
			*found = true;
			return parent;
		case 1:
			if (parent->rightChild == 0)
			{
				parent->rightChild =  allocTreeNode(0,0,0,data,numElements,0,0);
				*found = false;
				return parent->rightChild;
			}
			else
				return findNodeRec(parent->rightChild,data,numElements,found);
		case 2:
			if (parent->leftChild == 0)
			{
				parent->leftChild =  allocTreeNode(0,0,0,data,numElements,0,0);
				*found = false;
				return parent->leftChild;
			}
			else
				return findNodeRec(parent->leftChild,data,numElements,found);
	}
	// should never be reached
	return 0;
}

/**
 * Recursively find the node corresponding to state <data> in the state tree <root>,
 * or insert the node if it does not exist.
 * <numElements> is the size of the state vector <data>.
 * The return value of <found> indicates whether the element previously existed in the tree
 *
 * Returns the (possibly newly created) node belonging to <data>. If the tree is empty,
 * <root> is set to this node.
 */
StateTreeNode * findNode(StateTreeNode ** root, unsigned int * data, unsigned int numElements, bool * found)
{
	if (*root == 0)
	{
		*root = allocTreeNode(0,0,0,data,numElements,0,0);
		*found = false;
		return *root;
	}
	return findNodeRec(*root,data,numElements, found);
}

/**
 * Returns the successor of the supplied state node <current>.
 * If the state transition has not yet been calculated,
 * a new node is inserted into the tree <root>.
 * <numElementsPerEntry> is the number of array elements required to store one state.
 * <net> describes the network for which a state transition is performed.
 * <basinCounter> is a counter to be increased when a new state is identified
 */
StateTreeNode * findSuccessor(StateTreeNode ** root, StateTreeNode * current,
							 unsigned int numElementsPerEntry, BooleanNetwork * net, unsigned int * basinCounter)
{
	bool found;
	if (current->type.sync.successor == 0)
	// the state does not exist => calculate state transition and insert it
	{
		unsigned int nextState[numElementsPerEntry];
		stateTransition(current->data,nextState,net);
		current->type.sync.successor = findNode(root,nextState,numElementsPerEntry, &found);
 		++ *basinCounter;
	}
	return current->type.sync.successor;
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
	memcpy(&table[numElements* (*nodeNo)],root->type.sync.successor->data,numElements*sizeof(unsigned int));
	attractorAssignment[*nodeNo] = root->type.sync.attractorAssignment;
	stepsToAttractor[*nodeNo] = root->type.sync.stepsToAttractor;
	*nodeNo = *nodeNo + 1;

	if (root->rightChild != 0)
	// recursive descent of right subtree
		inOrderSerializeTree(root->rightChild,initialStates,table,attractorAssignment,
							 stepsToAttractor,numElements,nodeNo);
}

/**
 * Free a state tree supplied by <node> recursively.
 * If <freeSuccessorArray>, <successors> is assumed to be an array and freed
 */
void freeTreeNode(StateTreeNode * node, bool freeSuccessorArray)
{
	if (node->leftChild != 0)
		freeTreeNode(node->leftChild, freeSuccessorArray);
	if (node->rightChild != 0)
		freeTreeNode(node->rightChild, freeSuccessorArray);
	if (freeSuccessorArray)
		free(node->type.async.successors);
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
 * Retrieves attractors from a given transition table <table> with <numberOfStates> entries.
 * 
 * Returns a list of attractors - the last element of this list is empty!
 */
pAttractorInfo getAttractors(unsigned long long * table, unsigned int numberOfStates, unsigned int numberOfGenes)
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

	result->attractorList = attractorHead;

	free(table);

	return result;
}

/**
 * Retrieves attractors only for a given set of input states supplied in <selectedStates>.
 * Here, <ceil(net->numGenes / 32)> consecutive array entries describe one state, thus
 * the array size is <ceil(net->numGenes / 32) * numberOfStates>
 * <net> describes the network structure.
 */
pAttractorInfo getAttractorsForStates(unsigned int * selectedStates, unsigned int numberOfStates,
									  BooleanNetwork * net)
{
  unsigned long long i;
	unsigned int current_attractor = 0, elementsPerEntry;
	bool found;

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
		StateTreeNode * currentState = findNode(&stateTree,&selectedStates[i*elementsPerEntry],elementsPerEntry,&found);
		if(!currentState->type.sync.attractorAssignment)
		// the current state has not yet been visited
		{

			// first attractor has number 1
			current_attractor++;

			unsigned int steps = 0;
			unsigned int basinSize = 0;

			while(!currentState->type.sync.attractorAssignment)
			// iterate while no attractor has been assigned
			{
				++steps;

				// first simply count steps until attractor is reached
				// - to get the distance to the attractor, this number is
				// later subtracted from the maximum distance
				currentState->type.sync.stepsToAttractor = steps;
				currentState->type.sync.attractorAssignment = current_attractor;

				// perform a state transition
				currentState = findSuccessor(&stateTree,currentState,elementsPerEntry,net,&basinSize);
			}
			if(current_attractor == currentState->type.sync.attractorAssignment)
			//calculate length and basin size of new attractor
			{
				attractorList->basinSize = steps;

				// fix the number of steps to the attractor by calculating
				// the maximum distance and subtracting the current value from it
				int maxstep = currentState->type.sync.stepsToAttractor;

				int rec = 0;
				StateTreeNode * tmp = findNode(&stateTree,&selectedStates[i*elementsPerEntry],elementsPerEntry,&found);

				while(memcmp(tmp->data,currentState->data,elementsPerEntry*sizeof(unsigned int)))
				{
					rec++;
					tmp->type.sync.stepsToAttractor = maxstep - tmp->type.sync.stepsToAttractor;

					// perform a state transition
					tmp = findSuccessor(&stateTree,tmp,elementsPerEntry,net,&basinSize);
				}

				attractorList->length = steps - rec;

				attractorList->involvedStates = (unsigned int *) calloc(attractorList->length * elementsPerEntry,sizeof(unsigned int));
				attractorList->numElementsPerEntry = elementsPerEntry;

				int a=0;
				do
				{
					tmp->type.sync.stepsToAttractor = 0;
					memcpy(&attractorList->involvedStates[a],tmp->data,elementsPerEntry*sizeof(unsigned int));
					tmp = findSuccessor(&stateTree,tmp,elementsPerEntry,net,&basinSize);
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
				StateTreeNode * tmp = findNode(&stateTree,&selectedStates[i*elementsPerEntry],elementsPerEntry,&found);
				int maxstp = currentState->type.sync.stepsToAttractor + steps;

				while(memcmp(tmp->data,currentState->data,elementsPerEntry*sizeof(unsigned int)))
				{
					tmp->type.sync.attractorAssignment = currentState->type.sync.attractorAssignment;
					tmp->type.sync.stepsToAttractor = maxstp - tmp->type.sync.stepsToAttractor + 1;

					//perform a state transition
					tmp = findSuccessor(&stateTree,tmp,elementsPerEntry,net,&basinSize);
				}

				// update basin size in attractor record
				tmpList = attractorHead;

				unsigned int i;
				for (i = 1; i < currentState->type.sync.attractorAssignment; ++i)
					tmpList = tmpList->next;

				tmpList->basinSize = tmpList->basinSize + basinSize;
			}
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
	freeTreeNode(stateTree,false);
	return result;
}

/**
 * Identification of attractors in asynchronous networks
 */

typedef struct SSE
{
	unsigned int * state;
	struct SSE * next;
} StateStackElement;

/**
 * Push a new element on top of the stack <stack>.
 * <state> is the state to push onto the stack and has <numElements> elements.
 * Returns the new stack element.
 */
static inline StateStackElement * pushStateStackElement(StateStackElement ** stack,
												 unsigned int * state,
												 unsigned int numElements)
{
	StateStackElement * el = calloc(1,sizeof(StateStackElement));


	el->state = calloc(numElements,sizeof(unsigned int));
	memcpy(el->state,state,sizeof(unsigned int) * numElements);

	el->next = *stack;
	*stack = el;
	return el;
}

/**
 * Remove the top-level element from <stack>.
 */
static inline void deleteStateStackElement(StateStackElement ** stack)
{
	StateStackElement * el = *stack;
	*stack = (*stack)->next;
	free(el->state);
	free(el);
}


/**
 * Applies the transition function belonging to gene <geneIdx> to state <currentState>.
 * <net> holds the network information.
 * The result is returned in <currentState>.
 */
static inline void applySingleFunction(unsigned int * currentState, unsigned int geneIdx, BooleanNetwork * net)
{
	unsigned int k = 0;


	if (net->fixedGenes[geneIdx] == -1)
	// the gene is not fixed
	{
		unsigned long long inputdec = 0;

		for (k = net->inputGenePositions[geneIdx]; k < net->inputGenePositions[geneIdx+1]; k++)
		{
			if (net->inputGenes[k])
			// if the input of the function is not 0 (constant gene), take input bit
			{
				unsigned int gene = net->inputGenes[k] - 1;
				unsigned int bit;

				if (net->fixedGenes[gene] == -1)
					bit = (GET_BIT(currentState[gene / BITS_PER_BLOCK_32],
						   gene % BITS_PER_BLOCK_32));
				else
					// fixed genes are not encoded in the states
					// => take them from fixedGenes vector
					bit = net->fixedGenes[gene];
				inputdec |= bit	<< (net->inputGenePositions[geneIdx+1] - k - 1);
			}
		}
		// determine transition function
		int transition = net->transitionFunctions[net->transitionFunctionPositions[geneIdx] + inputdec];

		currentState[geneIdx / BITS_PER_BLOCK_32] = CLEAR_BIT(currentState[geneIdx / BITS_PER_BLOCK_32],
															  geneIdx % BITS_PER_BLOCK_32);

		if(transition != -1)
			// apply transition function
			currentState[geneIdx / BITS_PER_BLOCK_32] |= (transition << (geneIdx % BITS_PER_BLOCK_32));
		else
			// this is a dummy function for a constant gene
			// => value does not change
			currentState[geneIdx / BITS_PER_BLOCK_32] |= (GET_BIT(currentState[geneIdx / BITS_PER_BLOCK_32],
														geneIdx % BITS_PER_BLOCK_32) << (geneIdx % BITS_PER_BLOCK_32));

	}
}

/**
 * Calculate a random asynchronous state transition for <currentState>.
 * If <probabilities> is not NULL, this is a vector specifying
 * the cumulative distribution function of the probabilities of genes
 * to be chosen for a transition. Otherwise, each gene has equal probability.
 */
static inline void asynchronousStateTransition(unsigned int * currentState, double * probabilities,
										BooleanNetwork * net)
{
	unsigned int numElts, i;
	if (net->numGenes % BITS_PER_BLOCK_32 == 0)
		numElts = net->numGenes / BITS_PER_BLOCK_32;
	else
		numElts = net->numGenes / BITS_PER_BLOCK_32 + 1;

	if (probabilities == NULL)
	// uniform gene selection
	{
		unsigned int r;
		// ensure that no fixed gene is chosen
		do
		{
			r = intrand(net->numGenes);
		}
		while (net->fixedGenes[r] != -1);
		// make a transition with the chosen gene
		applySingleFunction(currentState,r,net);
	}
	else
	{
		double r = doublerand_1();

		// find the last index in the cumulative distribution that
		// is less than <r>
		for (i = 0; i < net->numGenes; ++i)
		{
			if (probabilities[i] < r && probabilities[i+1] >= r)
				break;
		}
		// make a transition with the chosen gene
		applySingleFunction(currentState,i,net);
	}
}

/**
 * Calculate the forward reachable set of <startState>.
 * <numElements> is the number of array elements used to represent one state.
 * If <avoidSelfLoops> is true, self loops are only considered if there are no other possible transitions.
 * <net> holds the network information.
 * <res> points to the root of the resulting state tree (set).
 * Returns the number of states in the set.
 */
unsigned int buildAsynchronousStateSet(unsigned int * startState, unsigned int numElements,
									   bool avoidSelfLoops, BooleanNetwork * net, StateTreeNode ** res)
{
	unsigned int startNodeCount = nodeCount;
	StateStackElement * stack = NULL;
	unsigned int i;
	bool found=false, newNodes=false;

	// push the start state onto the stack
	pushStateStackElement(&stack,startState,numElements);
	do
	// iterate while stack is not empty (depth-first search)
	{
		unsigned int origstate[numElements];

		memcpy(origstate,stack->state,sizeof(unsigned int) * numElements);

		// remove the top-level stack element
		deleteStateStackElement(&stack);

		StateTreeNode * current = findNode(res,origstate,numElements,&found);

		StateTreeNode ** successors;
		unsigned int numSuccessors;

		if (avoidSelfLoops)
		// try to find successor states that do not lead to the initial state
		{
			unsigned int successorStates[net->numGenes*numElements];
			for (i = 0; i < net->numGenes; ++i)
			// first, calculate all successors
			{
				memcpy(&successorStates[i*numElements],origstate,sizeof(unsigned int) * numElements);
				applySingleFunction(&successorStates[i*numElements],i,net);
			}

			unsigned int numNonSelfLoops = 0;
			bool noSelfLoop[net->numGenes];
			for (i = 0; i < net->numGenes; ++i)
			// now, check which of the successor states are the same as the initial state
			{
				if (memcmp(&successorStates[i*numElements],origstate,sizeof(unsigned int) * numElements) != 0)
				{
					++numNonSelfLoops;
					noSelfLoop[i] = true;
				}
				else
					noSelfLoop[i] = false;
			}
			if (numNonSelfLoops == 0)
			// all transitions are self loops
			// => accept self loop
			{
				successors = calloc(1,sizeof(StateTreeNode * ));
				numSuccessors = 1;
				successors[0] = findNode(res,successorStates,numElements,&found);
				if (!found)
					pushStateStackElement(&stack,successorStates,numElements);
			}
			else
			// there is at least one transition that is no self loop
			// => do not accept self loops
			{
				successors = calloc(numNonSelfLoops,sizeof(StateTreeNode * ));
				numSuccessors = numNonSelfLoops;
				unsigned int j;
				for (i = 0, j = 0; i < net->numGenes; ++i)
				{
					if (noSelfLoop[i])
					// create successor in tree
					{
						successors[j++] = findNode(res,&successorStates[i*numElements],numElements,&found);
						newNodes = newNodes | !found;
						if (!found)
							pushStateStackElement(&stack,&successorStates[i*numElements],numElements);
					}
				}
			}
		}
		else
		// self loops are allowed
		{
			unsigned int state[numElements];
			successors = calloc(net->numGenes,sizeof(StateTreeNode * ));
			numSuccessors = net->numGenes;
			for (i = 0; i < net->numGenes; ++i)
			// calculate all successors
			{
				memcpy(state,origstate,sizeof(unsigned int) * numElements);
				applySingleFunction(state,i,net);
				successors[i] = findNode(res,state,numElements,&found);
				newNodes = newNodes | !found;
				if (!found)
					pushStateStackElement(&stack,state,numElements);
			}
		}

		current->type.async.successors = successors;
		current->type.async.numSuccessors = numSuccessors;
	}
	while (stack != NULL);

	// return the number of elements in the state set
	return (nodeCount - startNodeCount);
}

/**
 * Recursively retrieve an array of states from a tree <root>
 * and store it to <states.
 * <numElements> is the number of array elements used to represent one state.
 * <nodeNo> is the current array entry to be written and should be initialized to 0.
 */
void getStateSet(StateTreeNode * root,
				 unsigned int * states,
				 unsigned int numElements,
				 unsigned int * nodeNo)
{

	if (root->leftChild != 0)
	// recursive descent of left subtree
		getStateSet(root->leftChild,states,numElements,nodeNo);

	// write the state itself
	memcpy(&states[numElements* (*nodeNo)],root->data,numElements*sizeof(unsigned int));
	*nodeNo = *nodeNo + 1;

	if (root->rightChild != 0)
	// recursive descent of right subtree
		getStateSet(root->rightChild,states,numElements,nodeNo);
}

/**
 * Recursively extract the transition table from a state set <root>
 * and store it to a list of transitions <table>.
 * <numElements> is the number of array elements used to represent one state.
 * <size> receives the size of the resulting table.
 */
void getLooseAttractorTransitionTable(StateTreeNode * root,
									 TransitionTableEntry ** table,
									unsigned int numElements,
									unsigned int * size)
{
	if (root->leftChild != 0)
	// recursive descent of left subtree
		getLooseAttractorTransitionTable(root->leftChild,table,numElements,size);

	bool duplicate[root->type.async.numSuccessors];
	memset(duplicate,0,sizeof(bool)*root->type.async.numSuccessors);

	unsigned int i, j;

	// check for duplicate transitions
	for (i = 0; i < root->type.async.numSuccessors; ++i)
	{
		for (j = i + 1; j < root->type.async.numSuccessors; ++j)
		{
			if (memcmp(&root->type.async.successors[i * numElements],
					   &root->type.async.successors[j * numElements],
					   sizeof(unsigned int) * numElements) == 0)
				duplicate[j] = true;
		}
	}

	// copy the unique transitions to the table
	for (i = 0; i < root->type.async.numSuccessors; ++i)
	{
		if (!duplicate[i])
		{
			insertNewTransition(table,root->data,root->type.async.successors[i]->data,numElements);
			++ *size;
		}
	}

	if (root->rightChild != 0)
	// recursive descent of right subtree
		getLooseAttractorTransitionTable(root->rightChild,table,numElements, size);
}

/**
 * Validate whether a set <attractor> with <attractorLength> states is a true attractor.
 * If <avoidSelfLoops> is true, self loops are only considered if there are no other possible transitions.
 * <net> holds the network information.
 */
bool validateAttractor(unsigned int * attractor, unsigned int attractorLength,
					   bool avoidSelfLoops,BooleanNetwork * net)
{
	unsigned int numElts, i;
	if (net->numGenes % BITS_PER_BLOCK_32 == 0)
		numElts = net->numGenes / BITS_PER_BLOCK_32;
	else
		numElts = net->numGenes / BITS_PER_BLOCK_32 + 1;

	for (i = 0; i < attractorLength; ++i)
	// iterate over states
	{
		StateTreeNode * set = NULL;

		// calculate forward reachable set of current state
		unsigned int size_set = buildAsynchronousStateSet(&attractor[i*numElts],numElts,avoidSelfLoops,net,&set);
		if (size_set != attractorLength)
		{
			freeTreeNode(set,true);
			return false;
		}

		unsigned int states_set[numElts * size_set];
		unsigned int nodeNo = 0;
		getStateSet(set,states_set,numElts,&nodeNo);
		freeTreeNode(set,true);

		// compare forward reachable set to original set
		if (memcmp(states_set,attractor,sizeof(unsigned int) * numElts * size_set) != 0)
			// no attractor
			return false;

	}
	return true;
}

/**
 * Calculate complex/loose attractors by performing <randomSteps> random transitions from
 * the <numberOfStates> states supplied in <selectedStates>.
 * If <avoidSelfLoops> is true, self loops are only considered if there are no other possible transitions.
 * If <probabilities> is not NULL, this vector holds the probabilities for each gene to be chosen
 * for a transition.
 */
pAttractorInfo getLooseAttractors(unsigned int * selectedStates, unsigned int numberOfStates,
					    BooleanNetwork * net, unsigned int randomSteps,
					    bool avoidSelfLoops, double * probabilities)
{
	// attractor list has empty dummy element at the end
	pAttractor attractorList = calloc(1,sizeof(Attractor));

	unsigned int numElts, i, j;
	if (net->numGenes % BITS_PER_BLOCK_32 == 0)
		numElts = net->numGenes / BITS_PER_BLOCK_32;
	else
		numElts = net->numGenes / BITS_PER_BLOCK_32 + 1;

	// if probabilities for the genes are supplied, exclude fixed genes (if any)
	// and recalculate probabilities
	double * pProbabilities = NULL;
	double convertedProbabilities[net->numGenes + 1];
	if (probabilities != NULL)
	{
		convertedProbabilities[0] = 0.0;
		double probabilitySum = 0.0;
		for (i = 0; i < net->numGenes; ++i)
		{
			if (net->fixedGenes[i] == -1)
				probabilitySum += probabilities[i];
		}
		for (i = 0; i < net->numGenes; ++i)
		{
			if (net->fixedGenes[i] == -1)
				convertedProbabilities[i+1] = convertedProbabilities[i]
				                               + probabilities[i]/probabilitySum;
			else
				convertedProbabilities[i+1] = convertedProbabilities[i];
		}
		pProbabilities = convertedProbabilities;
	}

	for (i = 0; i < numberOfStates; ++i)
	// iterate over supplied start states
	{
		unsigned int currentState[numElts];
		memcpy(currentState,&selectedStates[i*numElts],sizeof(unsigned int) * numElts);

		unsigned int t = 0;
		for (j = randomSteps; j > 0; --j)
		// perform <randomSteps> random state transitions
		// to reach a potential attractor
		{

			asynchronousStateTransition(currentState,pProbabilities,net);

			++t;
		}

		// calculate forward reachable set of end state
		StateTreeNode * set = NULL;
		unsigned int size_set = buildAsynchronousStateSet(currentState,numElts,avoidSelfLoops,net,&set);
		unsigned int states_set[numElts * size_set];

		unsigned int nodeNo = 0;
		getStateSet(set,states_set,numElts,&nodeNo);


		// search for the current potential attractor in the attractor list
		bool found = false;
		pAttractor current = attractorList;

		while (current != NULL && !found)
		{
			found = ((size_set == current->length)
					&& (memcmp(states_set,current->involvedStates,sizeof(unsigned int) * numElts * size_set) == 0));
			current = current->next;
		}


		if (!found)
		// the potential attractor does not yet exist in the result list
		{
			if (validateAttractor(states_set,size_set,avoidSelfLoops,net))
			// check whether this is a true attractor
			{
				pAttractor attractor = calloc(1,sizeof(Attractor));

				attractor->numElementsPerEntry = numElts;
				attractor->length = size_set;
				attractor->involvedStates = calloc(numElts * size_set,sizeof(unsigned int));
				memcpy(attractor->involvedStates,states_set,sizeof(unsigned int) * numElts * size_set);

				attractor->transitionTableSize = 0;

				if (size_set != 1)
				// if this is a steady-state attractor, no need to store transition table!
					getLooseAttractorTransitionTable(set,&attractor->table,numElts,&(attractor->transitionTableSize));

				attractor->next = attractorList;
				attractorList = attractor;
			}
		}
		freeTreeNode(set,true);
	}
	pAttractorInfo res = allocAttractorInfo(0,net->numGenes);
	res->attractorList = attractorList;
	return res;
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
 * startStates			An optional array of encoded states used to initialize the heuristics
 * 						(not needed for exhaustive search)
 * networkType			An integer that determines whether a synchronous or asynchronous search is performed
 * geneProbabilities	For asynchronous search, the probabilities of each gene to be chosen for update
 * randomSteps			For asynchronous search, the number of random transitions performed to reach a potential attractor
 * avoidSelfLoops		If set to true, self loops are only allowed if no other transitions are possible, which reduces the
 * 						number of edges in the attractors
 * returnTable			If set to true and networkType is synchronous, the transition table is included in the return value.
 */
SEXP getAttractors_R(SEXP inputGenes,
					 SEXP inputGenePositions,
					 SEXP transitionFunctions,
					 SEXP transitionFunctionPositions,
					 SEXP fixedGenes,
					 SEXP startStates,
					 SEXP networkType,
					 SEXP geneProbabilities,
					 SEXP randomSteps,
					 SEXP avoidSelfLoops,
					 SEXP returnTable)
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

	int _networkType = *INTEGER(networkType);
	int _randomSteps = *INTEGER(randomSteps);

	bool _returnTable = (bool)(*INTEGER(returnTable));
	bool _avoidSelfLoops = (bool)(*INTEGER(avoidSelfLoops));

	double * _probabilities = NULL;
	if (!isNull(geneProbabilities) && length(geneProbabilities) > 0)
		_probabilities = REAL(geneProbabilities);

	// count fixed genes, and create an index array for non-fixed genes:
	// <network.nonFixedGeneBits[i]> contains the bit positions in a state
	// at which the <i>-th gene is stored - this is different from <i>
	// as fixed genes are not stored
	unsigned int numNonFixed = 0, i;

	for(i = 0; i < network.numGenes; i++)
	{
		if(network.fixedGenes[i] == -1)
		{
			network.nonFixedGeneBits[i] = numNonFixed++;
		}
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
		res = getAttractors(table, numStates, network.numGenes);
	}
	else
	// start states supplied => only identify attractors for these states
	{
		unsigned int numElts;
		if (network.numGenes % BITS_PER_BLOCK_32 == 0)
			numElts = network.numGenes / BITS_PER_BLOCK_32;
		else
			numElts = network.numGenes / BITS_PER_BLOCK_32 + 1;

		unsigned int* _startStates = (unsigned int*) INTEGER(startStates);
	
		if (_networkType == SYNC_MODE)
		{
		  for (unsigned int i = 0; i < length(startStates) / numElts; ++i)
		  {
		    removeFixedGenes(&_startStates[i*numElts],network.fixedGenes,network.numGenes);
		  }
			res = getAttractorsForStates(_startStates, length(startStates) / numElts,
										&network);
		}
		else
		{
			res = getLooseAttractors(_startStates, length(startStates) / numElts,
								     &network,_randomSteps,
								     _avoidSelfLoops,_probabilities);
		}
	}

	// pack results in SEXP structure for return value
	SEXP resSXP;
	SEXP stateInfoSXP;
	int* array;

	// create a result list with two elements (attractors and transition table)
	PROTECT(resSXP = allocList(2));
	SET_TAG(resSXP, install("stateInfo"));
	SET_TAG(CDR(resSXP), install("attractors"));

	if (res->tableSize != 0 && _returnTable)
	{
		// create a 3-element list for the transition table
		PROTECT(stateInfoSXP = allocList(4));
		SET_TAG(stateInfoSXP, install("table"));
		SET_TAG(CDR(stateInfoSXP), install("attractorAssignment"));
		SET_TAG(CDR(CDR(stateInfoSXP)), install("stepsToAttractor"));
		SET_TAG(CDR(CDR(CDR(stateInfoSXP))), install("initialStates"));

		// write transition table result column
		SEXP tableSXP;
		PROTECT(tableSXP = allocVector(INTSXP,res->tableSize * res->numElementsPerEntry));
		array = INTEGER(tableSXP);
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
	}
	else
	{
		stateInfoSXP = R_NilValue;
	}

	// assign to result list
	SETCAR(resSXP,stateInfoSXP);

	if (res->tableSize != 0 && _returnTable)
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
		if (el->transitionTableSize == 0)
			PROTECT(attractorSXP = allocList(2));
		else
			PROTECT(attractorSXP = allocList(4));
		SET_TAG(attractorSXP, install("involvedStates"));
		SET_TAG(CDR(attractorSXP), install("basinSize"));

		if (el->transitionTableSize != 0)
		{
			SET_TAG(CDR(CDR(attractorSXP)), install("initialStates"));
			SET_TAG(CDR(CDR(CDR(attractorSXP))), install("nextStates"));
		}

		SEXP stateSXP;
		PROTECT(stateSXP = allocVector(INTSXP,el->length * el->numElementsPerEntry));
		array = INTEGER(stateSXP);
		for (i = 0; i < el->length; ++i)
		{
				if (_networkType == SYNC_MODE)
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

		if (el->transitionTableSize != 0)
		{
			SEXP attrInitialStateSXP;
			SEXP attrNextStateSXP;
			PROTECT(attrInitialStateSXP = allocVector(INTSXP,el->numElementsPerEntry * el->transitionTableSize));
			PROTECT(attrNextStateSXP = allocVector(INTSXP,el->numElementsPerEntry * el->transitionTableSize));
			unsigned int * initial = (unsigned int*)INTEGER(attrInitialStateSXP);
			unsigned int * next = (unsigned int*)INTEGER(attrNextStateSXP);

			TransitionTableEntry * currentState = el->table;
			for (i = 0; i < el->transitionTableSize; ++i)
			{
				memcpy(&initial[i * el->numElementsPerEntry],currentState->initialState,
					   sizeof(unsigned int) * el->numElementsPerEntry);
				memcpy(&next[i * el->numElementsPerEntry],currentState->nextState,
					   sizeof(unsigned int) * el->numElementsPerEntry);
				currentState = currentState->next;
			}
			SETCADDR(attractorSXP,attrInitialStateSXP);
			SETCADDDR(attractorSXP,attrNextStateSXP);
			UNPROTECT(2);
		}

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
