#ifndef ATTRACTOR_INFO_H
#define ATTRACTOR_INFO_H

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

	// the number of list elements in attractorList
	unsigned int numAttractors;
} AttractorInfo, *pAttractorInfo;

/**
 * Allocate a new AttractorInfo structure for <tableSize> states
 */
extern pAttractorInfo allocAttractorInfo(unsigned long long tableSize, unsigned int numGenes);

/**
 * Free a list of attractor structures
 */
extern void freeAttractorList(pAttractor p);

/**
 * Free an AttractorInfo structure including
 * all sub-elements and the attractor list
 */
extern void freeAttractorInfo(pAttractorInfo p);

#endif
