/**
 * C code for a generic symbolic simulator for Boolean networks
 *
 * This is part of the BoolNet R package.
 *
 * Copyright 2013 by Christoph Muessel
 *
 * Contact christoph.muessel@uni-ulm.de
 *
 */

#include "random.h"
#include "common.h"

#include "symbolic_network.h"

/**
 * Free the expression tree <formula>
 */
void freeFormula(BooleanFormula * formula)
{
	if (formula->type == FORMULA_ATOM || formula->type == FORMULA_CONSTANT)
	{
		free(formula);
	}
	else
	{
		BooleanOperator* operator = (BooleanOperator*) formula;

		unsigned int i;
		for (i = 0; i < operator->numOperands; ++i)
		{
			freeFormula(operator->operands[i]);
		}
		free(operator->operands);
		free(operator);
	}
}

/**
 * Free the internal network structure <network>
 */
void freeSymbolicNetwork(SymbolicBooleanNetwork * network)
{
	free(network->stateSizes);
	free(network->stateOffsets);
	free(network->fixedGenes);
	free(network->stateFixed);

	unsigned int i, j;
	for (i = 0; i < network->numGenes; ++i)
	{
		freeFormula(network->interactions[i]);
		if (network->cnfInteractions != NULL)
		{
			for (j = 0; j <= network->attractorSearchStartTime; ++j)
				freeFormula(network->cnfInteractions[i][j]);
			free(network->cnfInteractions[i]);
		}
	}
	free(network->interactions);
	if (network->cnfInteractions != NULL)
		free(network->cnfInteractions);
	free(network);
}

/**
 * Evaluate an expression tree <formula> on the previous state <state>.
 * <stateOffsets> is an array describing where the history of the i-th gene starts.
 * <numGenes> is the total number of genes in the network.
 * Returns a logical value resulting from the evaluation.
 */
unsigned char evaluate(BooleanFormula * formula, TemporalState * state,
		unsigned int * stateOffsets, const unsigned int numGenes)
{
	if (formula->type == FORMULA_CONSTANT)
	// return the constant value
	{
		Constant* constant = (Constant*) formula;
		char value = constant->value;

		if (constant->negated)
			return !value;
		else
			return value;
	}
	else if (formula->type == FORMULA_ATOM)
	// return the value of the atom in the current state
	{
		unsigned char res;
		BooleanAtom* atom = (BooleanAtom*) formula;

		res = state->state[stateOffsets[atom->literal] + atom->time];

		if (atom->negated)
			return !res;
		else
			return res;
	}
	else
	// recursively evaluate operands and summarize them using the operators
	{
		BooleanOperator* operator = (BooleanOperator*) formula;
		unsigned int i, threshold, sum, comparisonTimeStep;

		switch (operator->operator)
		{
		case OPERATOR_OR:
			for (i = 0; i < operator->numOperands; ++i)
			{
				if (evaluate(operator->operands[i], state, stateOffsets,
						numGenes))
				{
					return !operator->negated;
				}
			}
			return operator->negated;
		case OPERATOR_AND:
			for (i = 0; i < operator->numOperands; ++i)
			{
				if (!evaluate(operator->operands[i], state, stateOffsets,
						numGenes))
				{
					return operator->negated;
				}
			}
			return !operator->negated;
		case OPERATOR_MAJ:
			threshold = floor(operator->numOperands / 2.0);
			sum = 0;

			for (i = 0; i < operator->numOperands; ++i)
			{
				if (evaluate(operator->operands[i], state, stateOffsets,
						numGenes))
				{
					++sum;
					if (sum > threshold)
						return !operator->negated;
				}
			}
			return operator->negated;
		case OPERATOR_SUMIS:
			threshold = ((Constant *) operator->operands[operator->numOperands
					- 1])->value;
			sum = 0;

			for (i = 0; i < operator->numOperands - 1; ++i)
			{
				if (evaluate(operator->operands[i], state, stateOffsets,
						numGenes))
				{
					++sum;
				}
			}
			if (sum == threshold)
				return !operator->negated;
			else
				return operator->negated;
		case OPERATOR_SUMGT:
			threshold = ((Constant *) operator->operands[operator->numOperands
					- 1])->value;
			sum = 0;

			for (i = 0; i < operator->numOperands - 1; ++i)
			{
				if (evaluate(operator->operands[i], state, stateOffsets,
						numGenes))
				{
					++sum;
					if (sum > threshold)
						return !operator->negated;
				}
			}
			return operator->negated;
		case OPERATOR_SUMLT:
			threshold = ((Constant *) operator->operands[operator->numOperands
					- 1])->value;
			sum = 0;

			for (i = 0; i < operator->numOperands - 1; ++i)
			{
				if (evaluate(operator->operands[i], state, stateOffsets,
						numGenes))
				{
					++sum;
					if (sum >= threshold)
						return operator->negated;
				}
			}
			return !operator->negated;
		case OPERATOR_TIMEIS:
			comparisonTimeStep = ((Constant *) operator->operands[0])->value;
			if (operator->negated)
				return (state->timeStep != comparisonTimeStep - 1);
			else
				return (state->timeStep == comparisonTimeStep - 1);
		case OPERATOR_TIMEGT:
			comparisonTimeStep = ((Constant *) operator->operands[0])->value;
			if (operator->negated)
				return (state->timeStep <= comparisonTimeStep - 1);
			else
				return (state->timeStep > comparisonTimeStep - 1);
		case OPERATOR_TIMELT:
			comparisonTimeStep = ((Constant *) operator->operands[0])->value;
			if (operator->negated)
				return (state->timeStep >= comparisonTimeStep - 1);
			else
				return (state->timeStep < comparisonTimeStep - 1);
		default:
			Rf_error("Unknown operator!");
		}
	}
}

/**
 * Print a Boolean expression tree (for debugging purposes only)
 */
void printFormula(BooleanFormula * formula)
{
	if (formula->type == FORMULA_ATOM)
	{
		BooleanAtom* atom = (BooleanAtom*) formula;

		if (atom->negated)
			Rprintf("!");

		Rprintf("var%d", atom->literal);

		if (atom->time != 0)
			Rprintf("[%d]", -atom->time - 1);
	}
	else if (formula->type == FORMULA_CONSTANT)
	{
		Constant* constant = (Constant*) formula;

		if (constant->negated)
			Rprintf("!");

		Rprintf("%d", constant->value);
	}
	else
	{
		BooleanOperator* operator = (BooleanOperator*) formula;

		if (operator->negated)
			Rprintf("!");

		if (operator->operator == OPERATOR_MAJ)
			Rprintf("maj");
		else if (operator->operator == OPERATOR_SUMGT)
			Rprintf("sumgt");
		else if (operator->operator == OPERATOR_TIMEIS)
			Rprintf("timeis");
		else if (operator->operator == OPERATOR_TIMEGT)
			Rprintf("timegt");
		else if (operator->operator == OPERATOR_TIMELT)
			Rprintf("timelt");

		Rprintf("(");

		unsigned int i;
		for (i = 0; i < operator->numOperands; ++i)
		{
			printFormula(operator->operands[i]);

			if (i < operator->numOperands - 1)
			{
				if (operator->operator == OPERATOR_OR)
					Rprintf(" | ");
				else if (operator->operator == OPERATOR_AND)
					Rprintf(" & ");
				else
					Rprintf(", ");
			}
		}

		Rprintf(")");
	}
}

BooleanFormula * copyFormula(BooleanFormula * formula, bool negate, unsigned int timeOffset)
{
	if (formula->type == FORMULA_ATOM)
	{
		BooleanAtom * res = calloc(1, sizeof(BooleanAtom));
		memcpy(res, formula, sizeof(BooleanAtom));
		if (negate)
			res->negated = !res->negated;
		res->time += timeOffset;
		return (BooleanFormula *) res;
	}

	if (formula->type == FORMULA_CONSTANT)
	{
		Constant * res = calloc(1, sizeof(Constant));
		memcpy(res, formula, sizeof(Constant));
		if (negate)
			res->negated = !res->negated;
		return (BooleanFormula *) res;
	}

	unsigned int i;
	BooleanOperator * operator = (BooleanOperator*) formula;

	bool negated = (negate && !operator->negated)
			|| (!negate && operator->negated);
	BooleanOperator * res = allocOperator(operator->operator, negated,
			operator->numOperands, NULL);

	for (i = 0; i < res->numOperands; ++i)
	{
		res->operands[i] = copyFormula(operator->operands[i], false, timeOffset);
	}
	return (BooleanFormula *) res;
}

BooleanFormula * convertToCNF(BooleanFormula * formula, bool negate, unsigned int time)
{
	bool innerNegate = (!negate && formula->negated)
			|| (negate && !formula->negated);

	if (formula->type == FORMULA_ATOM || formula->type == FORMULA_CONSTANT)
	{
		return copyFormula(formula, negate, 0);
	}

	BooleanOperator* operator = (BooleanOperator*) formula;
	if (operator->operator == OPERATOR_TIMEIS
			|| operator->operator == OPERATOR_TIMELT
			|| operator->operator == OPERATOR_TIMEGT)
	// temporal operators => compare to time
	{
		int value;
		int comparisonTimeStep = ((Constant *) operator->operands[0])->value - 1;
		switch (operator->operator)
		{
			case OPERATOR_TIMEGT:
				value = (time > comparisonTimeStep);
				break;
			case OPERATOR_TIMELT:
				value = (time < comparisonTimeStep);
				break;
			default:
				value = (time == comparisonTimeStep);
				break;
		}
		if (innerNegate)
			value = !value;
		Constant * res = allocConstant(value,
				false);

		return (BooleanFormula *) res;
	}

	unsigned int operandBufferSize, numOperands, i, j;
	BooleanFormula ** operands;
	ALLOC_BUFFER(operands, BooleanFormula *, operandBufferSize, numOperands);

	if (operator->operator == OPERATOR_SUMGT
			|| operator->operator == OPERATOR_SUMLT
			|| operator->operator == OPERATOR_SUMIS
			|| operator->operator == OPERATOR_MAJ)
	// counting operators => transform to canonical CNF
	{
		unsigned int stateSize, threshold;
		if (operator->operator != OPERATOR_MAJ)
		{
			stateSize = operator->numOperands - 1;
			threshold = ((Constant*) operator->operands[stateSize])->value;
		}
		else
		{
			stateSize = operator->numOperands;
			threshold = floor(stateSize / 2.0);
		}
		unsigned char * state = calloc(stateSize, sizeof(unsigned char));

		do
		{
			bool includeTerm = false;
			unsigned int count = 0;
			for (i = 0; i < stateSize; ++i)
			{
				if (state[i])
					++count;
			}
			switch (operator->operator)
			{
			case OPERATOR_SUMGT:
			case OPERATOR_MAJ:
				includeTerm = (count <= threshold);
				break;
			case OPERATOR_SUMLT:
				includeTerm = (count >= threshold);
				break;
			case OPERATOR_SUMIS:
				includeTerm = (count != threshold);
				break;
			}
			if (innerNegate)
				includeTerm = !includeTerm;

			if (includeTerm)
			{
				BooleanFormula ** subOperands = calloc(stateSize,
						sizeof(BooleanFormula *));
				for (i = 0; i < stateSize; ++i)
				{
					subOperands[i] = copyFormula(operator->operands[i],
							state[i], 0);
				}
				BooleanOperator * subOp = allocOperator(OPERATOR_OR, false,
						stateSize, subOperands);
				PUT_BUFFER(operands, BooleanFormula *, operandBufferSize,
						numOperands, subOp);
			}
		} while (getNextState(state, NULL, stateSize));

		free(state);

		BooleanOperator * res = allocOperator(OPERATOR_AND, false, numOperands,
				operands);
		BooleanFormula * res2 = convertToCNF((BooleanFormula *) res, false, time);

		freeFormula((BooleanFormula *) res);
		return res2;

	}

	// AND or OR
	unsigned int innerAnds = 0, op;

	if (innerNegate)
	{
		if (operator->operator == OPERATOR_AND)
			op = OPERATOR_OR;
		else
			op = OPERATOR_AND;
	}
	else
		op = operator->operator;

	Constant * overrideResult = NULL;
	for (i = 0; i < operator->numOperands; ++i)
	{
		BooleanFormula * operand = convertToCNF(operator->operands[i],
				innerNegate, time);
		if (operand->type == FORMULA_OPERATOR)
		{
			BooleanOperator * subOp = (BooleanOperator *) operand;
			if (subOp->operator == op || subOp->numOperands == 1)
			{
				for (j = 0; j < subOp->numOperands; ++j)
				{
					PUT_BUFFER(operands, BooleanFormula *, operandBufferSize,
							numOperands, subOp->operands[j]);
				}

				free(subOp->operands);
				free(subOp);
			}
			else
			{
				if (subOp->operator == OPERATOR_AND)
					++innerAnds;
				PUT_BUFFER(operands, BooleanFormula *, operandBufferSize,
						numOperands, subOp);
			}
		}
		else
		{
			if (operand->type == FORMULA_CONSTANT)
			{
				Constant * c = (Constant * )operand;
				if ((op == OPERATOR_OR && ((c->value == 1) ^ c->negated)) ||
					(op == OPERATOR_AND && ((c->value == 0) ^ c->negated)))
					overrideResult = (Constant *)copyFormula((BooleanFormula *)c, false, 0);
				else
					freeFormula(operand);
			}
			else
				PUT_BUFFER(operands, BooleanFormula *, operandBufferSize,
						numOperands, operand);
		}
	}

	if (overrideResult != NULL)
	{
		for (i = 0; i < numOperands; ++i)
		{
			freeFormula(operands[i]);
		}
		free(operands);
		return (BooleanFormula *)overrideResult;
	}

	if (op == OPERATOR_OR && innerAnds > 0)
	{
		unsigned int transformedOperandBufferSize, numTransformedOperands;
		BooleanFormula ** transformedOperands;
		ALLOC_BUFFER(transformedOperands, BooleanFormula *,
				transformedOperandBufferSize, numTransformedOperands);
		unsigned int * comb = calloc(numOperands, sizeof(unsigned int));
		unsigned int * counts = calloc(numOperands, sizeof(unsigned int));

		for (i = 0; i < numOperands; ++i)
		{
			if (operands[i]->type == FORMULA_OPERATOR)
				counts[i] = ((BooleanOperator *) operands[i])->numOperands;
			else
				counts[i] = 1;
		}

		unsigned int pos = 0;
		while (true)
		{
			BooleanOperator * newOp = allocOperator(OPERATOR_OR, false, numOperands, NULL);
			for (i = 0; i < numOperands; ++i)
			{
				if (operands[i]->type == FORMULA_OPERATOR)

					newOp->operands[i] =
							copyFormula(
									((BooleanOperator *) operands[i])->operands[comb[i]],
									false, 0);
				else
					newOp->operands[i] = copyFormula(operands[i], false, 0);
			}
			PUT_BUFFER(transformedOperands, BooleanFormula *,
					transformedOperandBufferSize, numTransformedOperands, newOp);

			if (comb[pos] < counts[pos] - 1)
				++comb[pos];
			else
			{
				if (pos == numOperands - 1)
					break;

				while (pos < numOperands && comb[pos] == counts[pos] - 1)
					++pos;

				if (pos == numOperands)
					break;

				for (i = 0; i < pos; ++i)
				{
					comb[i] = 0;
				}
				++comb[pos];
				pos = 0;
			}
		}
		free(comb);
		free(counts);

		for (i = 0; i < numOperands; ++i)
			freeFormula(operands[i]);
		free(operands);

		BooleanOperator * res = allocOperator(OPERATOR_AND, false,
				numTransformedOperands, transformedOperands);
		BooleanFormula * res2 = convertToCNF((BooleanFormula *) res, false, time);

		freeFormula((BooleanFormula *) res);

		return res2;
	}
	else
	{
		BooleanOperator * res = allocOperator(op, false, numOperands, operands);
		return (BooleanFormula *) res;
	}
}

SEXP convertToCNF_R(SEXP network)
{
	SymbolicBooleanNetwork * _network = R_ExternalPtrAddr(network);

	if (_network == NULL)
		error("Internal network structures not supplied to C handler!");

	unsigned int k;
	for (k = 0; k < _network->numGenes; ++k)
	{
		Rprintf("var%d = ", k);
		BooleanFormula * converted = convertToCNF(_network->interactions[k],
		false, 0);
		printFormula(converted);
		freeFormula(converted);
		Rprintf("\n");
	}
	return R_NilValue;
}
