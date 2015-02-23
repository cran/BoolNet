#include "sat_search.h"

#include <stdbool.h>
#include <stddef.h>
#include <assert.h>

#include "attractor_info.h"
#include "common.h"
#include "random.h"

// Activate this switch to use the Lingeling solver instead of PicoSAT
// - requires lglib.c and lglib.h to be present in the package.
// Make sure to respect Lingeling's license terms!
//#define USE_LINGELING

#ifdef USE_LINGELING
#include "lglib.h"

#define sat_add_literal(solver, lit) lgladd(solver, lit)
#define sat_deref(solver, lit) lglderef(solver, lit)
#define sat_solve(solver) lglsat(solver)
#define sat_free(solver) lglrelease(solver)
#define SATISFIABLE LGL_SATISFIABLE

/**
 * Custom interrupt function to be used with Lingeling.
 * The function calls R's user interrupt check.
 */
void SATInterrupt(void * external_state)
{
	R_CheckUserInterrupt();
}

typedef LGL Solver;

#else
#include "picosat.h"

#define sat_add_literal(solver, lit) picosat_add(solver, lit)
#define sat_deref(solver, lit) picosat_deref(solver, lit)
#define sat_solve(solver) picosat_sat(solver, -1)
#define sat_free(solver) picosat_reset(solver)
#define SATISFIABLE PICOSAT_SATISFIABLE

typedef PicoSAT Solver;

/**
 * Custom interrupt function to be used with PicoSAT.
 * The function calls R's user interrupt check.
 */
int SATInterrupt(void * external_state)
{
	R_CheckUserInterrupt();
	return 0;
}

#endif
#include "symbolic_network.h"

/**
 * Custom memory allocator to be used with PicoSAT.
 * The allocator uses the internal memory map of BoolNet
 */
void * SATAlloc(void* mem, size_t sz)
{
	return CALLOC(sz, 1);
}

/**
 * Custom memory deallocator to be used with PicoSAT/Lingeling.
 * The deallocator uses the internal memory map of BoolNet
 */
void SATDealloc(void* mem, void * ptr, size_t sz)
{
	FREE(ptr);
}

/**
 * Custom memory resize function to be used with PicoSAT/Lingeling.
 * The function uses the internal memory map of BoolNet
 */
void * SATRealloc(void* mem, void * ptr, size_t old, size_t new)
{
	return REALLOC(ptr, new);
}

#ifdef USE_LINGELING
/**
 * Initialize the Lingeling SAT solver using appropriate memory handlers
 * and interrupt handlers.
 */
Solver * initSATSolver()
{
	Solver * sat = lglminit(NULL, &SATAlloc, &SATRealloc, &SATDealloc);
	lglonabort(sat, NULL, &SATInterrupt);
	lglsetopt(sat, "elmfull", 1);
	lglsetopt(sat, "elmschedprod", 1);
	lglsetopt(sat, "seed", intrand(~0)));
	return sat;
}
#else
/**
 * Initialize the PicoSAT SAT solver using appropriate memory handlers
 * and interrupt handlers.
 */
Solver * initSATSolver()
{
	Solver * sat = picosat_minit(NULL, &SATAlloc, &SATRealloc, &SATDealloc);
	picosat_set_interrupt(sat, NULL, &SATInterrupt);
	picosat_set_seed(sat, intrand(~0));
	return sat;
}
#endif

/**
 * Recursively encode a symbolic formula for the SAT solver. This function expects that
 * <formula> is already in CNF, being either a constant, an atom, or a conjunction
 * comprising only constants, atoms or flat disjunctions.
 * <network> is the network to which the formula belongs.
 * <geneIndex> is the index of the gene that is currently encoded.
 * <stateIndex> is the index of the state for which the formula should be encoded.
 * <sat> is the SAT solver.
 */
void encodeFormula(SymbolicBooleanNetwork * network, BooleanFormula * formula,
		int geneIndex, int stateIndex, Solver * sat)
{
	assert(formula->type != FORMULA_CONSTANT);
//	if (formula->type == FORMULA_CONSTANT)
//	{
//		Constant * constant = (Constant *) formula;
//		if ((constant->value == 1 && !constant->negated)
//				|| (constant->value == 0 && constant->negated))
//		{
//			// add tautology
//			picosat_add(sat, -(network->numGenes * stateIndex + geneIndex + 1));
//			picosat_add(sat, network->numGenes * stateIndex + geneIndex + 1);
////			picosat_assume(sat, network->numGenes * stateIndex + geneIndex + 1);
//		}
//		else
//		{
////			picosat_assume(sat, -(network->numGenes * stateIndex + geneIndex + 1));
//		}
//	}
//	else
	if (formula->type == FORMULA_ATOM)
	{
		BooleanAtom * atom = (BooleanAtom *) formula;

		// add the atom directly
		if (atom->negated)
		{
			sat_add_literal(sat,
					-(network->numGenes * (stateIndex + atom->time)
							+ atom->literal + 1));
		}
		else
		{
			sat_add_literal(sat,
					network->numGenes * (stateIndex + atom->time)
							+ atom->literal + 1);
		}
	}
	else if (formula->type == FORMULA_OPERATOR)
	{
		BooleanOperator * operator = (BooleanOperator *) formula;

		int i;
		for (i = 0; i < operator->numOperands; ++i)
		{
			// recursively encode operands
			encodeFormula(network, operator->operands[i], geneIndex, stateIndex,
					sat);
			if (operator->operator == OPERATOR_AND)
			// if this is the outer conjunction, each inner element must be a clause
			// and thus be terminated by a 0
			{
				sat_add_literal(sat, 0);
			}
		}
	}

}

/**
 * Add CNF clauses that represent the <stateIndex>th-last state of the network
 * <network> in the chain to the SAT solver <sat>.
 * <formulaIndex> specifies which of the formulae available for each gene should
 * be utilized for the state (usually the first, but there may be more than one
 * formula if the network comprises time-dependent predicates).
 */
static inline void addState_SAT(BooleanNetwork * network, Solver * sat,
		int stateIndex, unsigned int formulaIndex)
{
	unsigned int i, j, k;
	for (i = 0; i < network->numGenes; ++i)
	{
		if (network->fixedGenes[i] != -1)
		// fixed gene => add a clause with a constant 1 or 0 for that gene
		{
			if (network->fixedGenes[i] == 0)
				sat_add_literal(sat, -(network->numGenes * stateIndex + i + 1));
			else
				sat_add_literal(sat, network->numGenes * stateIndex + i + 1);
			sat_add_literal(sat, 0);
		}
		else
		// non-fixed gene
		{
			if (network->type == TRUTHTABLE_BOOLEAN_NETWORK)
			// network in truth table representation => construct CNF from truth table
			{
				TruthTableBooleanNetwork * ttNetwork =
						(TruthTableBooleanNetwork *) network;
				unsigned int numInputs = ttNetwork->inputGenePositions[i + 1]
						- ttNetwork->inputGenePositions[i];
				unsigned int funcSize = ttNetwork->transitionFunctionPositions[i
						+ 1] - ttNetwork->transitionFunctionPositions[i];
				// construct clauses for equivalence relations by "doubling" the truth table of
				// the transition function (corresponding to inactive and active output)
				for (j = 0; j < funcSize; ++j)
				{
					// in the first half, the target gene is inactive (truth table result is negated)
					// input genes must be negated for CNF
					if (ttNetwork->transitionFunctions[ttNetwork->transitionFunctionPositions[i]
							+ j])
					{
						sat_add_literal(sat,
								network->numGenes * stateIndex + i + 1);
						for (k = 0; k < numInputs; ++k)
						{
							int varIndex =
									network->numGenes * (stateIndex + 1)
											+ ttNetwork->inputGenes[ttNetwork->inputGenePositions[i]
													+ k];

							if (((j >> (numInputs - k - 1)) & 1) != 0)
								sat_add_literal(sat, -varIndex);
							else
								sat_add_literal(sat, varIndex);
						}
						sat_add_literal(sat, 0);
					}
					else
					{
						// in the second half, the target gene is active (unchanged truth table result)
						// input genes must be negated for CNF
						sat_add_literal(sat,
								-(network->numGenes * stateIndex + i + 1));
						for (k = 0; k < numInputs; ++k)
						{
							int varIndex =
									network->numGenes * (stateIndex + 1)
											+ ttNetwork->inputGenes[ttNetwork->inputGenePositions[i]
													+ k];

							if (((j >> (numInputs - k - 1)) & 1) != 0)
								sat_add_literal(sat, -varIndex);
							else
								sat_add_literal(sat, varIndex);
						}
						sat_add_literal(sat, 0);
					}
				}
			}
			else if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
			// network in symbolic representation => use CNF representation of the formulae
			// to generate the corresponding SAT solver clauses
			{
				SymbolicBooleanNetwork * symNetwork =
						(SymbolicBooleanNetwork *) network;

				// recursively encode the formula
				encodeFormula(symNetwork,
						symNetwork->cnfInteractions[i][formulaIndex], i,
						stateIndex, sat);

				if (symNetwork->cnfInteractions[i][formulaIndex]->type
						!= FORMULA_OPERATOR
						|| ((BooleanOperator *) symNetwork->cnfInteractions[i][formulaIndex])->operator
								!= OPERATOR_AND)
				// if the top-level formula is not an AND operator,
				// the result is a single clause that needs to be terminated here
				{
					sat_add_literal(sat, 0);
				}
			}
		}
	}
}

/**
 * Add clauses to the SAT solver <sat> to specify that the <state1>-th last state and
 * the <state2>-th last state of the chain should be equal.
 * <network> is the Boolean network from which the transition functions come.
 */
static inline void addEqualityCondition_SAT(Solver * sat,
		BooleanNetwork * network, unsigned int state1, unsigned int state2)
{
	unsigned int i, j;
	for (i = 0; i < network->numGenes; ++i)
	// iterate over genes
	{
		unsigned int historySize = 1;
		if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
			// determine history size
			historySize = ((SymbolicBooleanNetwork *) network)->stateSizes[i];

		// for temporal networks, the history must also be equal for the two states
		for (j = 0; j < historySize; ++j)
		{
			int gene1 = network->numGenes * (state1 + j) + i + 1;
			int gene2 = network->numGenes * (state2 + j) + i + 1;

			// equivalence consists of two clauses: (!A | B) & (A | !B)
			sat_add_literal(sat, gene1);
			sat_add_literal(sat, -gene2);
			sat_add_literal(sat, 0);

			sat_add_literal(sat, -gene1);
			sat_add_literal(sat, gene2);
			sat_add_literal(sat, 0);
		}
	}
}

/**
 * Exclude an identified attractor from the SAT search by adding its complement to the formula.
 * <sat> is a pointer to the SAT solver.
 * <network> is the employed Boolean network.
 * <attractor> is a structure holding the attractor to be excluded.
 */
static inline void excludeAttractor_SAT(Solver * sat, BooleanNetwork * network,
		pAttractor attractor)
{
	unsigned int l, i, j;
	for (l = 0; l < attractor->length; ++l)
	// iterate over attractor states
	{
		// add the complement of the current state as a clause
		// (for temporal networks, this must also include the state history)
		for (i = 0; i < network->numGenes; ++i)
		// iterate over genes
		{
			unsigned int historySize = 1;
			if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
				// determine history size of the current gene
				historySize =
						((SymbolicBooleanNetwork *) network)->stateSizes[i];
			for (j = 0; j < historySize; ++j)
			{
				// decode value of the gene
				bool val = GET_BIT_ARRAY(
						attractor->involvedStates[((attractor->length + l - j)
								% attractor->length)
								* attractor->numElementsPerEntry], i);

				// add its complement to the clause (at the last state of the chain)
				if (val)
					sat_add_literal(sat, -(j * network->numGenes + i + 1));
				else
					sat_add_literal(sat, j * network->numGenes + i + 1);
			}
		}
		// finish clause for the state
		sat_add_literal(sat, 0);
	}
}

/**
 * If a temporal network contains time-dependent predicates, this function adds a sequence of states
 * traversing the instable phase of the network to the beginning of the chain in <sat>. This is required
 * to ensure that only reachable start states are used for the attractor search.
 * <stateIndex> is index of the current beginning of the chain.
 */
static inline void addBurnInPhase(SymbolicBooleanNetwork * network,
		Solver * sat, int startStateIndex)
{
	if (network->attractorSearchStartTime > 0)
	// only do this if there is a time-dependent predicate
	{
		unsigned int j;
		for (j = 1; j <= network->attractorSearchStartTime; ++j)
			// iterate over the number of steps required for the burn-in and add states
			// utilizing the formulae that describe the behaviour at the corresponding time step
			addState_SAT((BooleanNetwork *) network, sat, startStateIndex + j,
					j);
	}
}

/**
 * Determine whether the SAT solver <sat> found an attractor. If yes, the attractor is
 * encoded in a Attractor structure and returned, and its states are excluded from future 
 * SAT searches.
 * <network> is the employed Boolean network.
 * <attractorLength> is the length of the attractor if already known, or -1 if the length
 * should be determined by the function.
 * <maxLength> is the total length of the state chain (which may be more than the attractor size).
 * Returns a pointer to an Attractor structure holding the attractor, or NULL if no attractor
 * was found in the state chain.
 */
static inline pAttractor getNextAttractor_SAT(Solver * sat,
		BooleanNetwork * network, int attractorLength, unsigned int maxLength)
{
	int i, j;

	if (attractorLength <= 0)
	// determine length if not already known
	{
		bool identical = true;
		for (attractorLength = 1; attractorLength <= maxLength;
				++attractorLength)
		// traverse the states from the next-to-first to the last to find a duplicate of the
		// first state (chain is reversed => first state is a potential attractor state)
		{
			identical = true;

			for (i = 0; i < network->numGenes; ++i)
			// compare the values of each gene in the current state and the first state
			{
				unsigned int historySize = 1;
				if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
					// determine history size of current gene
					historySize =
							((SymbolicBooleanNetwork *) network)->stateSizes[i];
				for (j = 0; j < historySize; ++j)
				// in temporal networks, do not only compare the current state,
				// but also its history to make sure this is really an attractor
				{
					if (sat_deref(sat,
							(attractorLength + j) * network->numGenes + i
									+ 1) != sat_deref(sat, j * network->numGenes + i + 1))
					{
						identical = false;
						break;
					}
				}
				if (!identical)
					break;

			}
			if (identical)
				break;
		}
		if (!identical)
			// none of the states was identical to the last one => no attractor in the chain
			return NULL;
	}

	// An attractor was found => build Attractor structure
	pAttractor res = CALLOC(1, sizeof(Attractor));
	if ((network->numGenes % BITS_PER_BLOCK_32) == 0)
		res->numElementsPerEntry = network->numGenes / BITS_PER_BLOCK_32;
	else
		res->numElementsPerEntry = network->numGenes / BITS_PER_BLOCK_32 + 1;
	res->length = attractorLength;
	res->involvedStates = (unsigned int *) CALLOC(
			res->length * res->numElementsPerEntry, sizeof(unsigned int));

	for (attractorLength = 0; attractorLength < res->length; ++attractorLength)
	// iterate over attractor states
	{
		for (i = 0; i < network->numGenes; ++i)
		// iterate over genes
		{
			// determine gene value and set the corresponding entry in the state array
			if (sat_deref(sat,
					(res->length - attractorLength - 1) * network->numGenes + i
					+ 1) == 1)
				SET_BIT_ARRAY(
						res->involvedStates[attractorLength
								* res->numElementsPerEntry], i);
		}
	}

// exclude the found attractor from future SAT searches
	excludeAttractor_SAT(sat, network, res);

	return res;
}

/**
 * Determine all attractors having at most <maxLength> states in the network <network>,
 * and return them in an AttractorInfo structure.
 */
pAttractorInfo getAttractors_SAT_maxLength(BooleanNetwork * network,
		unsigned int maxLength)
{
	if (network->type == SYMBOLIC_BOOLEAN_NETWORK
			&& ((SymbolicBooleanNetwork *) network)->attractorSearchStartTime
					> 0)
	{
		Rf_error(
				"SAT-based attractor search in networks with time-dependent predicates is only possible without attractor length restrictions!");
	}

	pAttractorInfo res = allocAttractorInfo(0, network->numGenes);
	// add dummy terminator node
	res->attractorList = (pAttractor) CALLOC(1, sizeof(Attractor));

	unsigned int i, j, length;

	unsigned int maxDelay = 1;

	if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
	// determine maximum time delay in the network
	{
		for (j = 0; j < network->numGenes; ++j)
		{
			unsigned int stateSize =
					((SymbolicBooleanNetwork *) network)->stateSizes[j];
			if (stateSize > maxDelay)
			{
				maxDelay = stateSize;
			}
		}
	}

	for (length = 1; length <= maxLength; ++length)
	// iterate over possible attractor lengths
	{

		// initialize SAT solver
		Solver * sat = initSATSolver();

		pAttractor attractor = res->attractorList;
		while (attractor->next != NULL)
		// exclude all previously identified attractors from the search
		{
			excludeAttractor_SAT(sat, network, attractor);
			attractor = attractor->next;
		}

		// create a chain of <length>+1 states
		for (i = 0; i <= length + maxDelay; ++i)
		{
			addState_SAT(network, sat, i, 0);
		}

#ifdef USE_LINGELING
		for (i = 0; i < network->numGenes; ++i)
		// freeze the first state(s) of the chain:
		// the variables of the first state are used
		// to exclude identified attractors from future searches
		// (multiple states for temporal Boolean networks with delays >1)
		{
			unsigned int stateSize = 1;
			if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
			// determine history size for current gene
				stateSize = ((SymbolicBooleanNetwork *) network)->stateSizes[i];

			for (j = 0; j < stateSize; ++j)
			{
				lglfreeze(sat, j * network->numGenes + i + 1);
			}
		}
#endif

		// specify that the first and the last state of the chain
		// should be identical
		addEqualityCondition_SAT(sat, network, 0, length);

		//lglprint(lgl, stdout);

		// call the SAT solver
		int satisfiable = sat_solve(sat);

		while (satisfiable == SATISFIABLE)
		// as long as the formula is satisfiable, there are more attractors
		// of the specified length
		{
			// encode the attractor in an Attractor structure,
			// prepend it to the result list, and exclude it from the search
			pAttractor attr = getNextAttractor_SAT(sat, network, length, 0);
			++res->numAttractors;
			attr->next = res->attractorList;
			res->attractorList = attr;

			// try to find further solutions
			satisfiable = sat_solve(sat);
		}
		// free SAT solver
		sat_free(sat);
	}
	return res;
}

/**
 * Identify all attractors using a SAT-based algorithm adapted from Dubrova et al., 2011.
 * Here, <network> is the network whose attractors are identified.
 * If <initialCycleSearchLength> > 0, the first step of the algorithm is to identify and exclude
 * all attractors of length 1 to <initialCycleSearchLength>.
 * <extensionMode> specifies the way the chain is extended:
 * EXTENSION_EXPONENTIAL corresponds to the exponential scheme by Dubrova et al.
 * EXTENSION_LINEAR corresponds to linearly increasing the length
 * EXTENSION_LINEAR_ADAPT corresponds to a linear increase whose step width is 
 * increased over the iterations
 * EXTENSION_MIXED corresponds to a mixture of linear and exponential increase,
 * where the chain length is doubled after 5 linear increases.
 * Returns an AttractorInfo structure comprising all attractors of the network.
 */
pAttractorInfo getAttractors_SAT_exhaustive(BooleanNetwork * network,
		unsigned int initialCycleSearchLength, unsigned int extensionMode)
{

	// initialize SAT solver
	Solver * sat = initSATSolver();

	int i, j, maxDelay = 1;

	if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
	// determine maximum time delay in the network
	{
		for (j = 0; j < network->numGenes; ++j)
		{
			unsigned int stateSize =
					((SymbolicBooleanNetwork *) network)->stateSizes[j];
			if (stateSize > maxDelay)
			{
				maxDelay = stateSize;
			}
		}
	}

	// check if the network comprises time-dependent predicates --
	// for such networks, the chain cannot be extended and must be rebuilt
	// each time a longer chain is needed
	bool isTimeDependent = (network->type == SYMBOLIC_BOOLEAN_NETWORK
			&& ((SymbolicBooleanNetwork *) network)->attractorSearchStartTime
					> 0);

	if (isTimeDependent)
	// for time-dependent predicates, length-based attractor search is not supported
	{
		initialCycleSearchLength = 0;
	}

	pAttractorInfo res;
	if (initialCycleSearchLength != 0)
	// first determine attractors by explicitly searching for cycles
	// if a cycle length is specified
	{
		//unsigned int maxAttractorLength = 2000/network->numGenes;
		//if (maxAttractorLength == 0)
		//  maxAttractorLength = 1;
		//else
		//if (maxAttractorLength > 10)
		//  maxAttractorLength = 10;
		res = getAttractors_SAT_maxLength(network, initialCycleSearchLength);
	}
	else
	// initialize empty attractor list
	{
		res = allocAttractorInfo(0, network->numGenes);
		// add dummy terminator node
		res->attractorList = (pAttractor) CALLOC(1, sizeof(Attractor));
	}

	// construct initial chain of maximum length 100
	int numStates = network->numGenes > 100 ? 100 : network->numGenes;

	for (i = 0; i < numStates; ++i)
	{
		addState_SAT(network, sat, i, 0);
	}
	//picosat_print(sat, stdout);

	if (initialCycleSearchLength != 0)
	// if attractors have been identified by the cycle search,
	// exclude them from the search
	{
		pAttractor attractor = res->attractorList;

		while (attractor->next != NULL)
		{
			excludeAttractor_SAT(sat, network, attractor);
			attractor = attractor->next;
		}
	}

#ifdef USE_LINGELING
	for (i = 0; i < network->numGenes; ++i)
	// freeze the first and the last state(s) of the chain:
	// the variables of the first state are used
	// to exclude identified attractors from future searches,
	// whereas the variables of the last state(s) are required for
	// an extension of the chain
	// (multiple states for temporal Boolean networks with delays >1)
	{
		unsigned int stateSize = 1;
		if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
		// determine history size of current gene
			stateSize = ((SymbolicBooleanNetwork *) network)->stateSizes[i];

		for (j = 0; j < stateSize; ++j)
		{
			lglfreeze(sat, j * network->numGenes + i + 1);

			if (!isTimeDependent)
				lglfreeze(sat, (numStates + j) * network->numGenes + i + 1);
		}
	}
#endif

	if (isTimeDependent)
		// add additional states in case of time-dependent predicates
		addBurnInPhase((SymbolicBooleanNetwork *) network, sat, numStates - 1);

	//picosat_print(sat, stdout);
	int satisfiable = sat_solve(sat);
	unsigned int numIncrements = 0;
	while (satisfiable == SATISFIABLE)
	// as long as the formula is satisfiable, there may be more attractors
	{
		// check whether the identified chain comprises an attractor,
		// and decode it if yes
		pAttractor attr = getNextAttractor_SAT(sat, network, -1,
				numStates - maxDelay + 1);

		if (attr == NULL)
		// no attractor found => chain is too short and must be prolonged
		{

			if (isTimeDependent)
			// if additional burn-in states were added, the chain must be rebuilt "from scratch",
			// i.e. the previous result cannot be recycled.
			{
				// free old solver
				sat_free(sat);

				// create new solver
				sat = initSATSolver();

				// exclude already identified attractors
				pAttractor attractor = res->attractorList;

				while (attractor->next != NULL)
				{
					excludeAttractor_SAT(sat, network, attractor);
					attractor = attractor->next;
				}

				// start chain at time point 0
				i = 0;
			}
			else
			{
#ifdef USE_LINGELING
				if (!isTimeDependent)
				{
					// melt the variables of the last state(s), as
					// the chain will be extended
					for (i = 0; i < network->numGenes; ++i)
					{
						unsigned int stateSize = 1;
						if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
							stateSize =
									((SymbolicBooleanNetwork *) network)->stateSizes[i];
						for (j = 0; j < stateSize; ++j)
							lglmelt(sat,
									(numStates + j) * network->numGenes + i
									 + 1);
					}

					// extend chain at previous time point
					i = numStates;
				}
#endif

			}

			switch (extensionMode)
			// determine the new chain length depending on the extension mode
			{
			case EXTENSION_EXPONENTIAL:
				// double the chain length
				numStates *= 2;
				break;
			case EXTENSION_LINEAR:
				// add the initial chain length
				numStates += network->numGenes > 100 ? 100 : network->numGenes;
				break;
			case EXTENSION_LINEAR_ADAPT:
				// add a multiple of the initial chain length depending
				// on the number of increments already performed
				numStates += (numIncrements++ / 5 + 1)
						* (network->numGenes > 100 ? 100 : network->numGenes);
				break;
			case EXTENSION_MIXED:
				// perform an exponential increase every five increments and a
				// linear increase otherwise
				if (++numIncrements % 5 == 0)
					numStates *= 2;
				else
					numStates +=
							network->numGenes > 100 ? 100 : network->numGenes;
				break;
			}

			// add new states to the chain
			for (; i < numStates; ++i)
			{
				addState_SAT(network, sat, i, 0);
			}

			if (!isTimeDependent)
			{
#ifdef USE_LINGELING
				// freeze the new end of the chain
				// (possibly multiple states in temporal Boolean networks)
				for (i = 0; i < network->numGenes; ++i)
				{
					unsigned int stateSize = 1;
					if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
						// determine history size of current gene
						stateSize =
								((SymbolicBooleanNetwork *) network)->stateSizes[i];
					for (j = 0; j < stateSize; ++j)
						lglfreeze(sat,
							(numStates + j) * network->numGenes + i + 1);

				}
#endif
			}
			else
			// network has time-dependent predicates
			{
				// add additional states in case of time-dependent predicates
				addBurnInPhase((SymbolicBooleanNetwork *) network, sat,
						numStates - 1);
#ifdef USE_LINGELING
				// as a new chain has been constructed from scratch, we must
				// again freeze the start state(s)
				for (i = 0; i < network->numGenes; ++i)
				{
					unsigned int stateSize = 1;
					if (network->type == SYMBOLIC_BOOLEAN_NETWORK)
					stateSize =
					((SymbolicBooleanNetwork *) network)->stateSizes[i];
					for (j = 0; j < stateSize; ++j)
					{
						lglfreeze(sat, j * network->numGenes + i + 1);
					}
				}
#endif
			}
		}
		else
		// an attractor has been found => add it
		{
			++res->numAttractors;
			attr->next = res->attractorList;
			res->attractorList = attr;
		}

		// look for more potential solutions
		satisfiable = sat_solve(sat);
	}
	// free the SAT solver
	sat_free(sat);
	return res;
}
