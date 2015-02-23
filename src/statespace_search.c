/**
 * C code to identify attractors in Boolean networks
 *
 * This is part of the BoolNet R package.
 *
 * Copyright 2009/2010 by Christoph Müssel and Zhou Dao
 *
 * Contact christoph.muessel@uni-ulm.de
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "random.h"
#include "common.h"
#include "boolean_network.h"
#include "statespace_search.h"
#include <time.h>

#define NODE_ARRAY_SIZE 1024

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

typedef struct
{
  StateTreeNode * root;
  
  unsigned int arraySize;
  unsigned int nodeCount;
  unsigned int numElements;
  
  unsigned int successorPos;
  
  ArrayListElement * nodeArrays;  
  ArrayListElement * dataArrays;
  
  ArrayListElement * successorArrays;  
} StateTree;

void stateTransition(unsigned int * currentState, unsigned int * nextState, TruthTableBooleanNetwork * net);

/**
 * Insert a new transition from <initialState> to <nextState> into <table>.
 * <numElements> is the number of array elements occupied by a state.
 */
TransitionTableEntry * insertNewTransition(TransitionTableEntry ** table,
										   unsigned int * initialState,
										   unsigned int * nextState,
										   unsigned int numElements)
{
	TransitionTableEntry * entry = (TransitionTableEntry *)CALLOC(1,sizeof(TransitionTableEntry));
	entry->initialState = CALLOC(numElements,sizeof(unsigned int));
	entry->nextState = CALLOC(numElements,sizeof(unsigned int));
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
		FREE(t->initialState);
		FREE(t->nextState);
		FREE(t);
		t = next;
	}
	while (t != NULL);
}


static inline StateTree * allocStateTree(unsigned int numElements, unsigned int arraySize)
{
  StateTree * tree = CALLOC(1, sizeof(StateTree));
  tree->root = NULL;
  tree->arraySize = arraySize;
  tree->successorPos = 0;
  tree->nodeCount = 0;
  tree->nodeArrays = NULL;
  tree->dataArrays = NULL;
  tree->successorArrays = NULL;
  tree->numElements = numElements;
  
  return(tree);
}

static inline void newNodeArray(StateTree * tree)
{
  allocNewArray(&tree->nodeArrays, tree->arraySize, sizeof(StateTreeNode));
  allocNewArray(&tree->dataArrays, tree->arraySize * tree->numElements, sizeof(unsigned int));
}

static inline void newSuccessorArray(StateTree * tree)
{
  allocNewArray(&tree->successorArrays, tree->arraySize, sizeof(StateTreeNode *));
  tree->successorPos = 0;
}

static inline void freeStateTree(StateTree * tree)
{
  freeArrayList(tree->nodeArrays);
  freeArrayList(tree->dataArrays);
  freeArrayList(tree->successorArrays);	
	
	FREE(tree);
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
static inline StateTreeNode * allocTreeNode(StateTree * tree,
                StateTreeNode * leftChild,
							  StateTreeNode * rightChild,
							  StateTreeNode * successor,
							  unsigned int * data,
							  unsigned int numElements,
							  unsigned int attractorAssignment,
							  unsigned int stepsToAttractor)
{
  if (tree->nodeCount % tree->arraySize == 0)
      newNodeArray(tree);
          
	StateTreeNode * res = &(((StateTreeNode *)tree->nodeArrays->array)[tree->nodeCount % tree->arraySize]);
	
	res->leftChild = leftChild;
	res->rightChild = rightChild;
	res->type.sync.successor = successor;
	res->data = &(((unsigned int *)tree->dataArrays->array)
	             [(tree->nodeCount % tree->arraySize) * tree->numElements]);
	             	
	memcpy(res->data,data,numElements*sizeof(unsigned int));

	res->type.sync.attractorAssignment = attractorAssignment;
	res->type.sync.stepsToAttractor = stepsToAttractor;
	++tree->nodeCount;
	return res;
}

static inline StateTreeNode ** reserveSuccessorArray(StateTree * tree, unsigned int numSuccessors)
{
  if (tree->successorArrays == NULL || tree->successorPos + numSuccessors >= tree->arraySize)
    newSuccessorArray(tree);
  
  StateTreeNode ** res = &(((StateTreeNode **)tree->successorArrays->array)[tree->successorPos]);
  tree->successorPos += numSuccessors;
  return res;
}

/**
 * Recursive helper function for findNode()
 */
StateTreeNode * findNodeRec(StateTree * tree, StateTreeNode * parent, unsigned int * data, unsigned int numElements, bool * found)
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
				parent->rightChild = allocTreeNode(tree, 0,0,0,data,numElements,0,0);
				*found = false;
				return parent->rightChild;
			}
			else
				return findNodeRec(tree, parent->rightChild,data,numElements,found);
		case 2:
			if (parent->leftChild == 0)
			{
				parent->leftChild = allocTreeNode(tree, 0,0,0,data,numElements,0,0);
				*found = false;
				return parent->leftChild;
			}
			else
				return findNodeRec(tree, parent->leftChild,data,numElements,found);
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
StateTreeNode * findNode(StateTree * tree, unsigned int * data, unsigned int numElements, bool * found)
{
	if (tree->root == 0)
	{
		tree->root = allocTreeNode(tree, 0,0,0,data,numElements,0,0);
		*found = false;
		return tree->root;
	}
	return findNodeRec(tree, tree->root,data,numElements, found);
}

/**
 * Returns the successor of the supplied state node <current>.
 * If the state transition has not yet been calculated,
 * a new node is inserted into the tree <root>.
 * <numElementsPerEntry> is the number of array elements required to store one state.
 * <net> describes the network for which a state transition is performed.
 * <basinCounter> is a counter to be increased when a new state is identified
 */
StateTreeNode * findSuccessor(StateTree * tree, StateTreeNode * current,
							 unsigned int numElementsPerEntry, TruthTableBooleanNetwork * net, unsigned int * basinCounter)
{
	bool found;
	if (current->type.sync.successor == 0)
	// the state does not exist => calculate state transition and insert it
	{
		unsigned int nextState[numElementsPerEntry];
		stateTransition(current->data,nextState,net);
		current->type.sync.successor = findNode(tree,nextState,numElementsPerEntry, &found);
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
 	R_CheckUserInterrupt();
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
/*
void freeTreeNode(StateTreeNode * node, bool freeSuccessorArray)
{
	if (node->leftChild != 0)
		freeTreeNode(node->leftChild, freeSuccessorArray);
	if (node->rightChild != 0)
		freeTreeNode(node->rightChild, freeSuccessorArray);
	if (freeSuccessorArray)
		FREE(node->type.async.successors);
	FREE(node->data);
	FREE(node);
	--nodeCount;
}
*/

/**
 * Make a transition from <currentState> to the next state.
 * States are encoded as arrays of ints and can thus contain an arbitrary number of genes.
 * <currentState> is a binary-coded integer with <numberOfGenes> used bits.
 * <fixedGenes> is an array of values specifying whether gene <i> is fixed (0 or 1) or not (-1).
 * <inputGenes> provides the input genes for all transition functions and can be split up
 * for a single function according to <inputGenePositions>.
 * <transitionFunctions> provides the truth tables for all transition functions and can be split up
 * for a single function according to <transitionFunctionPositions>.
 *
 * The return value is the next state, encoded in a single integer.
 */
void stateTransition(unsigned int * currentState, unsigned int * nextState, TruthTableBooleanNetwork * net)
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
				nextState[idx / BITS_PER_BLOCK_32] |= ((unsigned int)transition << (idx % BITS_PER_BLOCK_32));
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
 * Make a transition from <currentState> to the next state for a bit vector encoded
 * as a single integer (max. 64 genes).
 * <currentState> is a binary-coded integer with <numberOfGenes> used bits.
 * <fixedGenes> is an array of values specifying whether gene <i> is fixed (0 or 1) or not (-1).
 * <inputGenes> provides the input genes for all transition functions and can be split up
 * for a single function according to <inputGenePositions>.
 * <transitionFunctions> provides the truth tables for all transition functions and can be split up
 * for a single function according to <transitionFunctionPositions>.
 *
 * The return value is the next state, encoded in a single integer.
 */
unsigned long long stateTransition_singleInt(unsigned long long currentState, TruthTableBooleanNetwork * net)
{
	unsigned int i = 0, k = 0, idx = 0;

  unsigned long long nextState = 0;

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
						bit = (GET_BIT(currentState,
							     net->nonFixedGeneBits[gene]));
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
				nextState |= (transition << idx);
			else
				// this is a dummy function for a constant gene
				// => value does not change
				nextState |= (GET_BIT(currentState,
															idx) << idx);
			++idx;
		}
	}
	return nextState;
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
unsigned long long * getTransitionTable(TruthTableBooleanNetwork * net)
{
	unsigned long long i = 0;

	// determine number of fixed genes
	int numFixed = 0;
	for( i = 0; i < net->numGenes; i++)
		if(net->fixedGenes[i] != -1)
			++numFixed;

	int numNonFixed = net->numGenes - numFixed;

	// allocate truth table with 2^(non-fixed genes) elements
	unsigned long long numberOfElements = (unsigned long long)1 << numNonFixed;//pow(2,numNonFixed);
	unsigned long long * table = CALLOC(numberOfElements,sizeof(unsigned long long));
	if (table == 0)
	{
		Rf_error("Too few memory available!");
	}

	unsigned long long initialState = 0;

	// calculate state transitions
	for(initialState = 0; initialState < numberOfElements; ++initialState)
	{
    R_CheckUserInterrupt();
		//state is simply the binary encoding of the counter
		//calculate transition
	  table[initialState] = stateTransition_singleInt(initialState, net);			
	}
	return table;
}

/**
 * Retrieves attractors from a given transition table <table> with <numberOfStates> entries.
 * 
 * Returns a list of attractors - the last element of this list is empty!
 */
pAttractorInfo getAttractors(unsigned long long * table, unsigned long long numberOfStates, unsigned int numberOfGenes)
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
		//memcpy(&result->table[i],&table[i],sizeof(unsigned int) * elementsPerEntry);
		result->table[i*elementsPerEntry] = table[i] & 0xFFFFFFFF;
		
  	if (elementsPerEntry == 2)
  					result->table[i*elementsPerEntry + 1] = (table[i] & 0xFFFFFFFF00000000) >> 32;
	}

	pAttractor attractorHead, attractorList,tmpList;
	attractorHead = attractorList = (pAttractor) CALLOC(1,sizeof(Attractor));
	attractorList->next = NULL;

	for(i = 0; i < numberOfStates; i++)
	{
  	R_CheckUserInterrupt();
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

				attractorList->involvedStates = (unsigned int *) CALLOC(attractorList->length * elementsPerEntry,sizeof(unsigned int));
				attractorList->numElementsPerEntry = elementsPerEntry;

				int a=0;
				do
				{
					result->stepsToAttractor[tmp] = 0;
					//attractorList->involvedStates[a++] = tmp;
					//memcpy(&attractorList->involvedStates[a],&tmp,elementsPerEntry*sizeof(unsigned int));
				
					// get low-order and high-order longs in a platform-independent manner
 					attractorList->involvedStates[a] = tmp & 0xFFFFFFFF;
				
					if (elementsPerEntry == 2)
  					attractorList->involvedStates[a+1] = (tmp & 0xFFFFFFFF00000000) >> 32;
					
					tmp = table[tmp];
					a += elementsPerEntry;
				}
				while(tmp != begin); /* set steps of attractors to 0; add attractor stub information */

				//generate a next attractor space

				attractorList->next = (pAttractor)CALLOC(1,sizeof(Attractor));
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
	result->numAttractors = current_attractor - 1;

	FREE(table);

	return result;
}

/**
 * Retrieves attractors only for a given set of input states supplied in <selectedStates>.
 * Here, <ceil(net->numGenes / 32)> consecutive array entries describe one state, thus
 * the array size is <ceil(net->numGenes / 32) * numberOfStates>
 * <net> describes the network structure.
 */
pAttractorInfo getAttractorsForStates(unsigned int * selectedStates, unsigned int numberOfStates,
									  TruthTableBooleanNetwork * net)
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
	StateTree * stateTree = allocStateTree(elementsPerEntry, NODE_ARRAY_SIZE);

	pAttractor attractorHead, attractorList,tmpList;
	attractorHead = attractorList = (pAttractor) CALLOC(1,sizeof(Attractor));
	attractorList->next = NULL;
	
	for(i = 0; i < numberOfStates; i++)
	{
	  // check whether the current state is already in the state tree, otherwise insert it
		StateTreeNode * currentState = findNode(stateTree,&selectedStates[i*elementsPerEntry],elementsPerEntry,&found);
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
		  	R_CheckUserInterrupt();
				++steps;

				// first simply count steps until attractor is reached
				// - to get the distance to the attractor, this number is
				// later subtracted from the maximum distance
				currentState->type.sync.stepsToAttractor = steps;
				currentState->type.sync.attractorAssignment = current_attractor;

				// perform a state transition
				currentState = findSuccessor(stateTree,currentState,elementsPerEntry,net,&basinSize);
			}
			if(current_attractor == currentState->type.sync.attractorAssignment)
			//calculate length and basin size of new attractor
			{
				attractorList->basinSize = steps;

				// fix the number of steps to the attractor by calculating
				// the maximum distance and subtracting the current value from it
				int maxstep = currentState->type.sync.stepsToAttractor;

				int rec = 0;
				StateTreeNode * tmp = findNode(stateTree,&selectedStates[i*elementsPerEntry],elementsPerEntry,&found);

				while(memcmp(tmp->data,currentState->data,elementsPerEntry*sizeof(unsigned int)))
				{
			  	R_CheckUserInterrupt();
					rec++;
					tmp->type.sync.stepsToAttractor = maxstep - tmp->type.sync.stepsToAttractor;

					// perform a state transition
					tmp = findSuccessor(stateTree,tmp,elementsPerEntry,net,&basinSize);
				}

				attractorList->length = steps - rec;

				attractorList->involvedStates = (unsigned int *) CALLOC(attractorList->length * elementsPerEntry,sizeof(unsigned int));
				attractorList->numElementsPerEntry = elementsPerEntry;

				int a=0;
				do
				{
			  	R_CheckUserInterrupt();
					tmp->type.sync.stepsToAttractor = 0;
					memcpy(&attractorList->involvedStates[a],tmp->data,elementsPerEntry*sizeof(unsigned int));
					tmp = findSuccessor(stateTree,tmp,elementsPerEntry,net,&basinSize);
					a += elementsPerEntry;
				}
				while(memcmp(tmp->data,currentState->data,elementsPerEntry*sizeof(unsigned int)));

				//generate a next attractor space
				attractorList->next = (pAttractor)CALLOC(1,sizeof(Attractor));
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
				StateTreeNode * tmp = findNode(stateTree,&selectedStates[i*elementsPerEntry],elementsPerEntry,&found);
				int maxstp = currentState->type.sync.stepsToAttractor + steps;

				while(memcmp(tmp->data,currentState->data,elementsPerEntry*sizeof(unsigned int)))
				{
			  	R_CheckUserInterrupt();
					tmp->type.sync.attractorAssignment = currentState->type.sync.attractorAssignment;
					tmp->type.sync.stepsToAttractor = maxstp - tmp->type.sync.stepsToAttractor + 1;

					//perform a state transition
					tmp = findSuccessor(stateTree,tmp,elementsPerEntry,net,&basinSize);
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

	pAttractorInfo result = allocAttractorInfo(stateTree->nodeCount,net->numGenes);
	result->attractorList = attractorHead;
	result->initialStates = CALLOC(result->tableSize * elementsPerEntry,sizeof(unsigned int));
	result->numAttractors = current_attractor - 1;

	unsigned int nodeNo = 0;

	// build a series of arrays by in-order traversing the tree
	inOrderSerializeTree(stateTree->root,
						 result->initialStates,
						 result->table,
						 result->attractorAssignment,
						 result->stepsToAttractor,
						 elementsPerEntry,
						 &nodeNo);
	freeStateTree(stateTree);
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
	StateStackElement * el = CALLOC(1,sizeof(StateStackElement));


	el->state = CALLOC(numElements,sizeof(unsigned int));
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
  
	FREE(el->state);
	FREE(el);
}


/**
 * Applies the transition function belonging to gene <geneIdx> to state <currentState>.
 * <net> holds the network information.
 * The result is returned in <currentState>.
 */
static inline void applySingleFunction(unsigned int * currentState, unsigned int geneIdx, TruthTableBooleanNetwork * net)
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
										TruthTableBooleanNetwork * net)
{
	unsigned int i;

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
 * Returns the resulting state tree (set)..
 */
StateTree * buildAsynchronousStateSet(unsigned int * startState, unsigned int numElements,
									   bool avoidSelfLoops, TruthTableBooleanNetwork * net)
{
  StateTree * res = allocStateTree(numElements, NODE_ARRAY_SIZE);
	StateStackElement * stack = NULL;
	unsigned int i;
	bool found=false, newNodes=false;

	// push the start state onto the stack
	pushStateStackElement(&stack,startState,numElements);
	do
	// iterate while stack is not empty (depth-first search)
	{
  	R_CheckUserInterrupt();
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
				successors = reserveSuccessorArray(res,1);
				numSuccessors = 1;
				successors[0] = findNode(res,successorStates,numElements,&found);
				if (!found)
					pushStateStackElement(&stack,successorStates,numElements);
			}
			else
			// there is at least one transition that is no self loop
			// => do not accept self loops
			{
				successors = reserveSuccessorArray(res, numNonSelfLoops);
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
			successors = reserveSuccessorArray(res,net->numGenes);
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
	return res;
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
					   bool avoidSelfLoops,TruthTableBooleanNetwork * net)
{
	unsigned int numElts, i;
	if (net->numGenes % BITS_PER_BLOCK_32 == 0)
		numElts = net->numGenes / BITS_PER_BLOCK_32;
	else
		numElts = net->numGenes / BITS_PER_BLOCK_32 + 1;

	for (i = 0; i < attractorLength; ++i)
	// iterate over states
	{
		StateTree * set = buildAsynchronousStateSet(&attractor[i*numElts],numElts,avoidSelfLoops,net);
		if (set->nodeCount != attractorLength)
		{
			freeStateTree(set);
			return false;
		}

		unsigned int size_set = set->nodeCount;
		unsigned int states_set[numElts * size_set];
		unsigned int nodeNo = 0;
		getStateSet(set->root,states_set,numElts,&nodeNo);
		freeStateTree(set);

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
					    TruthTableBooleanNetwork * net, unsigned int randomSteps,
					    bool avoidSelfLoops, double * probabilities)
{
	// attractor list has empty dummy element at the end
	pAttractor attractorList = CALLOC(1,sizeof(Attractor));

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
	unsigned int attractorCount = 0;

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
		StateTree * set = buildAsynchronousStateSet(currentState,numElts,avoidSelfLoops,net);
		unsigned int states_set[numElts * set->nodeCount];

		unsigned int nodeNo = 0;
		getStateSet(set->root,states_set,numElts,&nodeNo);

		// search for the current potential attractor in the attractor list
		bool found = false;
		pAttractor current = attractorList;

		while (current != NULL && !found)
		{
			found = ((set->nodeCount == current->length)
					&& (memcmp(states_set,current->involvedStates,sizeof(unsigned int) * numElts * set->nodeCount) == 0));
			current = current->next;
		}

		if (!found)
		// the potential attractor does not yet exist in the result list
		{
			if (validateAttractor(states_set,set->nodeCount,avoidSelfLoops,net))
			// check whether this is a true attractor
			{
				pAttractor attractor = CALLOC(1,sizeof(Attractor));

				attractor->numElementsPerEntry = numElts;
				attractor->length = set->nodeCount;
				attractor->involvedStates = CALLOC(numElts * set->nodeCount,sizeof(unsigned int));
				memcpy(attractor->involvedStates,states_set,sizeof(unsigned int) * numElts * set->nodeCount);

				attractor->transitionTableSize = 0;

				if (set->nodeCount != 1)
				// if this is a steady-state attractor, no need to store transition table!
					getLooseAttractorTransitionTable(set->root,&attractor->table,numElts,&(attractor->transitionTableSize));

				++attractorCount;
				attractor->next = attractorList;
				attractorList = attractor;
			}
		}
		freeStateTree(set);
	}
	pAttractorInfo res = allocAttractorInfo(0,net->numGenes);
	res->numAttractors = attractorCount;
	res->attractorList = attractorList;
	return res;
}
