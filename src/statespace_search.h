#ifndef STATESPACE_SEARCH_H
#define STATESPACE_SEARCH_H
#include "boolean_network.h"
#include "attractor_info.h"
#include <stdbool.h>

/**
 * Retrieves attractors only for a given set of input states supplied in <selectedStates>.
 * Here, <ceil(net->numGenes / 32)> consecutive array entries describe one state, thus
 * the array size is <ceil(net->numGenes / 32) * numberOfStates>
 * <net> describes the network structure.
 */
extern pAttractorInfo getAttractorsForStates(unsigned int * selectedStates, 
                                             unsigned int numberOfStates,
                               							 TruthTableBooleanNetwork * net);

                               							 

/**
 * Retrieves attractors from a given transition table <table> with <numberOfStates> entries.
 * 
 * Returns a list of attractors - the last element of this list is empty!
 */
extern pAttractorInfo getAttractors(unsigned long long * table, unsigned long long numberOfStates, unsigned int numberOfGenes);

/**
 * Calculate complex/loose attractors by performing <randomSteps> random transitions from
 * the <numberOfStates> states supplied in <selectedStates>.
 * If <avoidSelfLoops> is true, self loops are only considered if there are no other possible transitions.
 * If <probabilities> is not NULL, this vector holds the probabilities for each gene to be chosen
 * for a transition.
 */
pAttractorInfo getLooseAttractors(unsigned int * selectedStates, unsigned int numberOfStates,
					    TruthTableBooleanNetwork * net, unsigned int randomSteps,
					    bool avoidSelfLoops, double * probabilities);
					    
/**
 * Retrieves the result column of the state transition table.
 * <numberOfGenes> specifies the total number of genes.
 * <fixedGenes> is an array of values specifying whether gene <i> is fixed (0 or 1) or not (-1).
 * <inputGenes> provides the input genes for all transition functions and can be split up
 * for a single function according to <inputGenePositions>.
 * <transitionFunctions> provides the truth tables for all transition functions and can be split up
 * for a single function according to <transitionFunctionPositions>.
 */
unsigned long long * getTransitionTable(TruthTableBooleanNetwork * net);					    

#endif
