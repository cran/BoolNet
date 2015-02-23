#ifndef SYMBOLIC_SIMULATOR_H
#define SYMBOLIC_SIMULATOR_H

#include <stdlib.h>
#include <stdbool.h>
#include "boolean_network.h"

#define FORMULA_ATOM 0
#define FORMULA_OPERATOR 1
#define FORMULA_CONSTANT 2

#define OPERATOR_AND 0
#define OPERATOR_OR 1
#define OPERATOR_MAJ 2
#define OPERATOR_SUMIS 3
#define OPERATOR_SUMGT 4
#define OPERATOR_SUMLT 5
#define OPERATOR_TIMEIS 6
#define OPERATOR_TIMEGT 7
#define OPERATOR_TIMELT 8

#define MODE_EXHAUSTIVE 0
#define MODE_SUPPLIED 1
#define MODE_RANDOM 2

#define HASH_UNSET ~0

/**
 * "Base class" for Boolean Formulae
 */
typedef struct
{
	// The formula type which allows
	// for casting to the corresponding "derived" type
	unsigned char type;

	// Is the formula negated or not?
	bool negated;
} BooleanFormula;

/**
 * "Derived class" for operators
 */
typedef struct
{
	// here: type = FORMULA_OPERATOR
	unsigned char type;

	bool negated;

	// The operator type
	unsigned char operator;

	// The number of operands
	unsigned int numOperands;

	// The operands
	BooleanFormula ** operands;
} BooleanOperator;

/**
 * "Derived class" for literals
 */
typedef struct
{
	// here: type = FORMULA_ATOM
	unsigned char type;

	bool negated;

	// The index of the variable that this object represents
	int literal;

	// The temporal difference
	unsigned int time;

} BooleanAtom;

/**
 * "Derived class" for constants
 */
typedef struct
{
	// here: type = FORMULA_CONSTANT
	unsigned char type;

	bool negated;

	// The value of the constant
	int value;
} Constant;

/**
 * A structure holding a symbolic network
 * with all its network trees
 */
typedef struct
{
	// here: type = SYMBOLIC_BOOLEAN_NETWORK
	unsigned char type;

	// The number of genes in the network
	unsigned int numGenes;

	// A vector specifying whether genes are fixed (0/1) or not (-1)
	int * fixedGenes;

	// The formulae representing the transition functions
	BooleanFormula ** interactions;

	// The formulae representing the transition functions in CNF
	// (for SAT-based attractor search)
	// Usually, this is a single formula per gene.
	// In the presence of time-dependent predicates, there are
	// <attractorSearchStartTime> + 1 formulae per gene.
	BooleanFormula *** cnfInteractions;

	// A vector specifying the maximum delay for each gene
	unsigned int * stateSizes;

	// The sum of all elements in stateSizes
	unsigned int totalStateSize;

	// The point of time at which attractor search can start,
	// as temporal predicates do not change any more
	unsigned int attractorSearchStartTime;

	// A vector specifying the indices in the state vector
	// at which the histories for each of the genes start
	unsigned int * stateOffsets;

	// A vector marking each entry in the state vector as fixed or not,
	// corresponding to the genes it belongs to
	int * stateFixed;

} SymbolicBooleanNetwork;

/**
 * A structure holding a state for temporal network simulation
 */
typedef struct
{
	// The index of the start state to which this state belongs
	unsigned long long startState;

	// The current time step
	unsigned int timeStep;

	// The state vector (also containing previous states for certain genes if required)
	unsigned char state[];
} TemporalState;

extern void freeSymbolicNetwork(SymbolicBooleanNetwork * network);

static inline BooleanOperator * allocOperator(unsigned int operator,
bool negated, unsigned int numOperands, BooleanFormula ** operands)
{
	BooleanOperator * res = calloc(1, sizeof(BooleanOperator));
	res->negated = negated;
	res->type = FORMULA_OPERATOR;
	res->operator = operator;
	res->numOperands = numOperands;
	if (operands == NULL)
		res->operands = calloc(numOperands, sizeof(BooleanFormula *));
	else
		res->operands = operands;
	return res;
}

static inline BooleanAtom * allocAtom(int literal, unsigned int time,
		bool negated)
{
	BooleanAtom * res = calloc(1, sizeof(BooleanAtom));
	res->type = FORMULA_ATOM;
	res->negated = negated;
	res->literal = literal;
	res->time = time;
	return res;
}

static inline Constant * allocConstant(int value, bool negated)
{
	Constant * res = calloc(1, sizeof(Constant));
	res->type = FORMULA_CONSTANT;
	res->negated = negated;
	res->value = value;
	return res;
}

/**
 * Free the expression tree <formula>
 */
extern void freeFormula(BooleanFormula * formula);

/**
 * Free the internal network structure <network>
 */
extern void freeSymbolicNetwork(SymbolicBooleanNetwork * network);

extern BooleanFormula * copyFormula(BooleanFormula * formula, bool negate, unsigned int timeOffset);

/**
 * Evaluate an expression tree <formula> on the previous state <state>.
 * <stateOffsets> is an array describing where the history of the i-th gene starts.
 * <numGenes> is the total number of genes in the network.
 * Returns a logical value resulting from the evaluation.
 */
extern unsigned char evaluate(BooleanFormula * formula, TemporalState * state,
		unsigned int * stateOffsets, const unsigned int numGenes);

extern BooleanFormula * convertToCNF(BooleanFormula * formula, bool negate, unsigned int time);

#endif
