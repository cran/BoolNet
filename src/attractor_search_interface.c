#include <R.h>
#include <Rinternals.h>
#include "statespace_search.h"
#include "sat_search.h"
#include "common.h"

/**
 * Identification of attractors
 */

#define SYNC_MODE_STATE_SPACE 0
#define ASYNC_MODE_RANDOM 1
#define SYNC_MODE_SAT_EXHAUSTIVE 2
#define SYNC_MODE_SAT_RESTRICTED 3

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
 * simType			An integer that determines whether a synchronous, asynchronous or SAT-based search is performed
 * geneProbabilities	For asynchronous search, the probabilities of each gene to be chosen for update
 * randomSteps			For asynchronous search, the number of random transitions performed to reach a potential attractor
 * avoidSelfLoops		If set to true, self loops are only allowed if no other transitions are possible, which reduces the
 * 						number of edges in the attractors
 * returnTable			If set to true and networkType is synchronous, the transition table is included in the return value.
 * maxAttractorSize		The maximum attractor length for SAT-based search
 */
SEXP getAttractors_R(SEXP inputGenes, SEXP inputGenePositions,
		SEXP transitionFunctions, SEXP transitionFunctionPositions,
		SEXP fixedGenes, SEXP startStates, SEXP simType,
		SEXP geneProbabilities, SEXP randomSteps, SEXP avoidSelfLoops,
		SEXP returnTable, SEXP maxAttractorSize)
{
	GetRNGstate();

	// decode information in SEXP for use in C

	TruthTableBooleanNetwork network;
	network.type = TRUTHTABLE_BOOLEAN_NETWORK;
	network.numGenes = length(fixedGenes);

	network.inputGenes = INTEGER(inputGenes);
	network.inputGenePositions = INTEGER(inputGenePositions);
	network.transitionFunctions = INTEGER(transitionFunctions);
	network.transitionFunctionPositions = INTEGER(transitionFunctionPositions);
	network.fixedGenes = INTEGER(fixedGenes);
	network.nonFixedGeneBits = CALLOC(network.numGenes, sizeof(unsigned int));

	int _networkType = *INTEGER(simType);
	int _randomSteps = *INTEGER(randomSteps);

	bool _returnTable = (bool) (*INTEGER(returnTable));
	bool _avoidSelfLoops = (bool) (*INTEGER(avoidSelfLoops));

	double * _probabilities = NULL;
	if (!isNull(geneProbabilities) && length(geneProbabilities) > 0)
		_probabilities = REAL(geneProbabilities);

	// count fixed genes, and create an index array for non-fixed genes:
	// <network.nonFixedGeneBits[i]> contains the bit positions in a state
	// at which the <i>-th gene is stored - this is different from <i>
	// as fixed genes are not stored
	unsigned int numNonFixed = 0, i;

	for (i = 0; i < network.numGenes; i++)
	{
		if (network.fixedGenes[i] == -1)
		{
			network.nonFixedGeneBits[i] = numNonFixed++;
		}
	}
	pAttractorInfo res;

	if (_networkType == SYNC_MODE_STATE_SPACE
			&& (isNull(startStates) || length(startStates) == 0))
	// no start states supplied => perform exhaustive search
	{
		// calculate transition table
		unsigned long long * table = getTransitionTable(&network);

		if (table == 0)
		{
			PutRNGstate();
			return R_NilValue;
		}

		unsigned long long numStates = (unsigned long long) 1 << numNonFixed;//pow(2,numNonFixed);
		// find attractors
		res = getAttractors(table, numStates, network.numGenes);
	}
	else
	// start states supplied or SAT search
	{
		unsigned int numElts;
		if (network.numGenes % BITS_PER_BLOCK_32 == 0)
			numElts = network.numGenes / BITS_PER_BLOCK_32;
		else
			numElts = network.numGenes / BITS_PER_BLOCK_32 + 1;

		unsigned int* _startStates = (unsigned int*) INTEGER(startStates);

		if (_networkType == SYNC_MODE_STATE_SPACE)
		{
			for (unsigned int i = 0; i < length(startStates) / numElts; ++i)
			{
				removeFixedGenes(&_startStates[i * numElts], network.fixedGenes,
						network.numGenes);
			}
			res = getAttractorsForStates(_startStates,
					length(startStates) / numElts, &network);
		}
		else if (_networkType == SYNC_MODE_SAT_EXHAUSTIVE)
		{
			int _maxAttractorSize;
			if (isNull(maxAttractorSize))
			{
				_maxAttractorSize = 1;
			}
			else
				_maxAttractorSize = *INTEGER(maxAttractorSize);
			res =
					getAttractors_SAT_exhaustive(
							(BooleanNetwork *) &network,
							_maxAttractorSize,
							EXTENSION_MIXED);
		}
		else if (_networkType == SYNC_MODE_SAT_RESTRICTED)
		{
			res = getAttractors_SAT_maxLength((BooleanNetwork *) &network,
					*INTEGER(maxAttractorSize));
		}
		else
		{
			res = getLooseAttractors(_startStates,
					length(startStates) / numElts, &network, _randomSteps,
					_avoidSelfLoops, _probabilities);
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
		PROTECT(
				tableSXP = allocVector(INTSXP,res->tableSize * res->numElementsPerEntry));
		array = INTEGER(tableSXP);
		for (i = 0; i < res->tableSize; ++i)
		{
			// the transition table results do not contain fixed genes => insert them
			insertFixedGenes(&res->table[i * res->numElementsPerEntry],
					network.fixedGenes, network.numGenes);
			memcpy(&array[i * res->numElementsPerEntry],
					&res->table[i * res->numElementsPerEntry],
					res->numElementsPerEntry * sizeof(unsigned int));
		}
		SETCAR(stateInfoSXP, tableSXP);
		UNPROTECT(1);

		// write attractor assignment vector for states in transition table
		SEXP assignmentSXP;
		PROTECT(assignmentSXP = allocVector(INTSXP,res->tableSize));
		array = INTEGER(assignmentSXP);
		memcpy(array, res->attractorAssignment, res->tableSize * sizeof(int));
		SETCADR(stateInfoSXP, assignmentSXP);
		UNPROTECT(1);

		// write a vector with number of transitions from a state to an attractor
		SEXP stepSXP;
		PROTECT(stepSXP = allocVector(INTSXP,res->tableSize));
		array = INTEGER(stepSXP);
		memcpy(array, res->stepsToAttractor, res->tableSize * sizeof(int));
		SETCADDR(stateInfoSXP, stepSXP);
		UNPROTECT(1);

		// if available, write the original states
		SEXP initialStateSXP;
		if (res->initialStates == 0)
			initialStateSXP = R_NilValue;
		else
		// if start states are specified, the initial states for each transition have to be saved as well,
		// as they cannot be inferred by enumeration
		{
			PROTECT(
					initialStateSXP = allocVector(INTSXP,res->tableSize * res->numElementsPerEntry));
			array = INTEGER(initialStateSXP);
			for (i = 0; i < res->tableSize; ++i)
			{
				// the transition table results do not contain fixed genes => insert them
				insertFixedGenes(
						&res->initialStates[i * res->numElementsPerEntry],
						network.fixedGenes, network.numGenes);
				memcpy(&array[i * res->numElementsPerEntry],
						&res->initialStates[i * res->numElementsPerEntry],
						res->numElementsPerEntry * sizeof(unsigned int));
			}
			SETCADDDR(stateInfoSXP, initialStateSXP);
			UNPROTECT(1);
		}
	}
	else
	{
		stateInfoSXP = R_NilValue;
	}

	// assign to result list
	SETCAR(resSXP, stateInfoSXP);

	if (res->tableSize != 0 && _returnTable)
		UNPROTECT(1);

	// count attractors
	unsigned int numAttractors = 0;
	pAttractor el;

	for (el = res->attractorList; el->next != NULL; el = el->next)
		++numAttractors;

	// write attractors
	SEXP attractorsSXP;
	PROTECT(attractorsSXP = allocList(numAttractors));
	SEXP listPos = attractorsSXP;
	for (el = res->attractorList, i = 0; el->next != NULL; el = el->next, ++i)
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
		PROTECT(
				stateSXP = allocVector(INTSXP,el->length * el->numElementsPerEntry));
		array = INTEGER(stateSXP);
		for (i = 0; i < el->length; ++i)
		{
			if (_networkType == SYNC_MODE_STATE_SPACE)
				// insert fixed gene values, as they are missing in the calculated results
				insertFixedGenes(
						&el->involvedStates[i * el->numElementsPerEntry],
						network.fixedGenes, network.numGenes);

			memcpy(&array[i * el->numElementsPerEntry],
					&el->involvedStates[i * el->numElementsPerEntry],
					el->numElementsPerEntry * sizeof(unsigned int));

		}
		SETCAR(attractorSXP, stateSXP);

		// write basin size
		SEXP basinSXP;
		PROTECT(basinSXP = allocVector(INTSXP,1));
		array = INTEGER(basinSXP);
		array[0] = el->basinSize;
		SETCADR(attractorSXP, basinSXP);
		SETCAR(listPos, attractorSXP);
		if (el->next != NULL)
			listPos = CDR(listPos);

		if (el->transitionTableSize != 0)
		{
			SEXP attrInitialStateSXP;
			SEXP attrNextStateSXP;
			PROTECT(
					attrInitialStateSXP = allocVector(INTSXP,el->numElementsPerEntry * el->transitionTableSize));
			PROTECT(
					attrNextStateSXP = allocVector(INTSXP,el->numElementsPerEntry * el->transitionTableSize));
			unsigned int * initial = (unsigned int*) INTEGER(
					attrInitialStateSXP);
			unsigned int * next = (unsigned int*) INTEGER(attrNextStateSXP);

			TransitionTableEntry * currentState = el->table;
			for (i = 0; i < el->transitionTableSize; ++i)
			{
				memcpy(&initial[i * el->numElementsPerEntry],
						currentState->initialState,
						sizeof(unsigned int) * el->numElementsPerEntry);
				memcpy(&next[i * el->numElementsPerEntry],
						currentState->nextState,
						sizeof(unsigned int) * el->numElementsPerEntry);
				currentState = currentState->next;
			}
			SETCADDR(attractorSXP, attrInitialStateSXP);
			SETCADDDR(attractorSXP, attrNextStateSXP);
			UNPROTECT(2);
		}

		UNPROTECT(3);
	}
	UNPROTECT(1);
	SETCADR(resSXP, attractorsSXP);

	PutRNGstate();
	UNPROTECT(1);

	// free resources
	freeAttractorInfo(res);

	FREE(network.nonFixedGeneBits);

	return (resSXP);
}
