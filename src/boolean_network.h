#ifndef BOOLEAN_NETWORK_H
#define BOOLEAN_NETWORK_H

#define TRUTHTABLE_BOOLEAN_NETWORK 0
#define PROBABILISTIC_BOOLEAN_NETWORK 1
#define SYMBOLIC_BOOLEAN_NETWORK 2

/**
 * Basic structure for different types of
 * Boolean network (extended by derived structures
 * TruthTableBooleanNetwork, ProbabilisticBooleanNetwork and SymbolicBooleanNetwork)
 */
typedef struct {
	// the network type
	unsigned char type;

	// the number of genes in the network
	unsigned int numGenes;

	// a vector specifying whether the genes are fixed:
	// -1 means the gene is not fixed, 1 and 0 means the
	// genes are fixed to the corresponding values
	int * fixedGenes;

} BooleanNetwork;

/**
 * Internal structure describing a Boolean network
 * with a truth table representation
 */
typedef struct
{
	// here: type = TRUTHTABLE_BOOLEAN_NETWORK
	unsigned char type;

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

} TruthTableBooleanNetwork;


typedef struct
{
	int * inputGenes;
	int * transitionFunction;

	unsigned int numGenes;

	double probability;

	unsigned int functionIndex;
} PBNFunction;

/**
 * Internal structure describing a Probabilistic Boolean network
 */
typedef struct
{
	// here: type = PROBABILISTIC_BOOLEAN_NETWORK
	unsigned char type;

	// the number of genes in the network
	unsigned int numGenes;

	// a vector specifying whether the genes are fixed:
	// -1 means the gene is not fixed, 1 and 0 means the
	// genes are fixed to the corresponding values
	int * fixedGenes;

	// the number of non-fixed genes in the network
	unsigned int numNonFixedGenes;

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

// see also symbolic_boolean_network.h for the definition of symbolic networks

#endif
