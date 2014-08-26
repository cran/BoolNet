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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "random.h"
#include "common.h"
#include <time.h>

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

typedef struct
{
  unsigned char type;
  bool negated;  
} BooleanFormula;

typedef struct
{
  unsigned char type;
  bool negated;  
  
  unsigned char operator;
  unsigned int numOperands;
  BooleanFormula ** operands;
} BooleanOperator;

typedef struct
{
  unsigned char type;
  bool negated;
  
  int literal;
  unsigned int time;
      
} BooleanAtom;

typedef struct
{
  unsigned char type;
  bool negated;  
  int value;      
} Constant;

typedef struct 
{
  unsigned int numGenes;  
  //char ** geneNames;  
  BooleanFormula ** interactions;
  
  unsigned int * stateSizes;
  unsigned int totalStateSize;
  unsigned int attractorSearchStartTime;
  
  unsigned int * stateOffsets;
  int * fixed;
  int * stateFixed;  
    
} SymbolicBooleanNetwork;

typedef struct TS
{
  unsigned long long startState;
  unsigned int timeStep;
  unsigned char state[];
} TemporalState;

typedef struct ASLE
{
  TemporalState **  states;
  unsigned int numStates;
  struct ASLE * next;
  unsigned int index;  
} AttractorStateListElement;

typedef struct
{
  AttractorStateListElement * head;
  AttractorStateListElement * tail;  
  unsigned int size;
} AttractorStateList;

typedef struct SH
{
  UT_hash_handle hh;
   //unsigned char * nextState;
   // unsigned char initialState[];
  struct SH * nextState;
  unsigned int sequenceIndex;
  
  TemporalState * initialState;
} StateHash;

typedef struct 
{
  UT_hash_handle hh;
  AttractorStateListElement * attractor;
  TemporalState * state;
} AttractorHash;


typedef struct 
{
  ArrayListElement * hashStructs;
  ArrayListElement * stateStructs;
  
  StateHash * table;
  unsigned int stateSize;
  unsigned int internalStateSize;
  unsigned int hashSize;  
  
  unsigned int currentEntryHash;
  unsigned int poolArraySize;
}
StateHashTable;

typedef struct 
{
  ArrayListElement * hashStructs;
  ArrayListElement * stateStructs;
  
  AttractorHash * table;
  unsigned int stateSize;
  unsigned int internalStateSize;
  
  unsigned int currentEntryHash;
  unsigned int poolArraySize;
}
AttractorHashTable;


void freeSymbolicNetwork(SymbolicBooleanNetwork * network);
static void finalizeSymbolicNetwork(SEXP network);
void printFormula(BooleanFormula * formula);

#define STATEHASH_ARRAY_SIZE 1024

/**
 * Functions to handle hash tables for states and attractors
 */

StateHashTable * allocStateHashTable(unsigned int stateSize)
{
  StateHashTable * res = CALLOC(1,sizeof(StateHashTable));
  
  res->hashStructs = NULL;
  res->stateStructs = NULL;
  
  res->stateSize = stateSize;
  
  // ensure proper alignment of states in array
  if (stateSize % sizeof(unsigned long long) == 0)
    res->internalStateSize = stateSize;
  else
    res->internalStateSize = ((stateSize / sizeof(unsigned long long)) * 
                              sizeof(unsigned long long) + sizeof(unsigned long long));
  
  res->hashSize = offsetof(TemporalState, state) +
                  stateSize * sizeof(unsigned char);
  
  res->table = NULL;
  res->currentEntryHash = 0;
  res->poolArraySize = STATEHASH_ARRAY_SIZE;
  return res;
}

AttractorHashTable * allocAttractorHashTable(unsigned int stateSize)
{
  AttractorHashTable * res = CALLOC(1,sizeof(AttractorHashTable));
  
  res->hashStructs = NULL;
  res->stateStructs = NULL;
  
  res->stateSize = stateSize;
  
  // ensure proper alignment of states in array
  if (stateSize % sizeof(unsigned long long) == 0)
    res->internalStateSize = stateSize;
  else
    res->internalStateSize = ((stateSize / sizeof(unsigned long long)) * 
                              sizeof(unsigned long long) + sizeof(unsigned long long));
  
  res->table = NULL;
  res->currentEntryHash = 0;
  res->poolArraySize = STATEHASH_ARRAY_SIZE;
  return res;
}

void freeStateHashTable(StateHashTable * hash)
{
  StateHash * entry, * tmp;
  HASH_ITER(hh, hash->table, entry, tmp) 
  {
    HASH_DEL(hash->table, entry);
  }
  freeArrayList(hash->hashStructs);
  freeArrayList(hash->stateStructs);    
  FREE(hash);
}

void freeAttractorHashTable(AttractorHashTable * hash)
{
  AttractorHash * entry, * tmp;
  HASH_ITER(hh, hash->table, entry, tmp) 
  {
    HASH_DEL(hash->table, entry);
  }
  freeArrayList(hash->hashStructs);  
  freeArrayList(hash->stateStructs);    
  FREE(hash);
}

/**
 * Finalizer that ensures that the internal C structures
 * of the R network object <network> are freed when the
 * garbage collector frees <network>.
 */
static void finalizeSymbolicNetwork(SEXP network)
{
    if (NULL == R_ExternalPtrAddr(network))
        return;
    //Rprintf("finalizing\n");
    SymbolicBooleanNetwork * _network = R_ExternalPtrAddr(network);
    
    freeSymbolicNetwork(_network);
    
    R_ClearExternalPtr(network);
}

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
  free(network->fixed);
  free(network->stateFixed);  
  
  unsigned int i;
  for (i = 0; i < network->numGenes; ++i)
  {
    freeFormula(network->interactions[i]);
  }
  free(network->interactions);
  free(network);
}

/**
 * Decode a sequence of states <state> as a state-memory representation and
 * returns it in <result>.
 * <network> is the network to which the state belongs.
 * <stateLength> is the number of entries in the state vector.
 */
inline static void decodeState(SymbolicBooleanNetwork * network, 
                        int * state, 
                        unsigned int stateLength, 
                        TemporalState * result)
{
  unsigned int numStates =  stateLength/network->numGenes;
  unsigned int i, j, k;
  
  memset(result,0,sizeof(TemporalState) + network->totalStateSize * sizeof(unsigned char));

  for (i = 0; i < network->numGenes; ++i)
  {
    for (j = 0; j < network->stateSizes[i]; ++j)
    {
      
      if (j < numStates)
        result->state[network->stateOffsets[i] + j] = state[i * numStates + j];
      else
        result->state[network->stateOffsets[i] + j] = state[(i + 1) * numStates - 1];

      //Rprintf("Gene %d state %d Idx %d Val %d\n", i, j, i * numStates + j, result->state[network->stateOffsets[i] + j]);        
    }
  }
    
}

/** 
 * Check whether a stat transition pair has already been traversed, and store it in the hash
 * table otherwise.
 * <hash> is the hash table for the lookup.
 * <nextState> is the next state after the state transition.
 * <initialState> receives the initial state of the transition if it was found.
 * If <dropTimeInfo> is true, the absolute time point is not stored in the table.
 * <sequenceIndex> is the index of the sequence to which the transition belongs.
 * Returns true if the state was found, and false otherwise.
 */
bool findOrStore(StateHashTable * hash, StateHash ** initialState, TemporalState * nextState, 
                 bool dropTimeInfo, unsigned int sequenceIndex)
{
  //TemporalState * search = CALLOC(1, sizeof(StateHash) + sizeof(unsigned char) * hash->stateSize);
  //memcpy(search, nextState, sizeof(TemporalState) +  sizeof(unsigned char) * hash->stateSize);
  
  //if (dropTimeInfo)
  //{
  //  search->timeStep = HASH_UNSET;
  //  search->startState = HASH_UNSET;    
  //}  
  
  StateHash * res = NULL;
  
 
  if (dropTimeInfo)
  {
    HASH_FIND(hh, hash->table, nextState->state, hash->stateSize, res);
  }
  else
  {
    HASH_FIND(hh, hash->table, nextState, hash->hashSize, res);
  }

  if (res != NULL)
  {
    if (*initialState != NULL)
      (*initialState)->nextState = res;
    *initialState = res;
    return true;
  }
  
  if (hash->currentEntryHash % hash->poolArraySize == 0)
  {
	    allocNewArray(&hash->hashStructs, hash->poolArraySize, 
	                  sizeof(StateHash));
      allocNewArray(&hash->stateStructs, hash->poolArraySize, 
	                  sizeof(TemporalState) + sizeof(unsigned char) * hash->internalStateSize);  
  }
  res = (StateHash *) (((unsigned char * )hash->hashStructs->array) + sizeof(StateHash)
                                                  * (hash->currentEntryHash % hash->poolArraySize));
  res->initialState = (TemporalState *) (((unsigned char * )hash->stateStructs->array) + (sizeof(TemporalState) +
                                   sizeof(unsigned char) * hash->internalStateSize) 
                                                  * (hash->currentEntryHash % hash->poolArraySize));
  
  memcpy(res->initialState, nextState, sizeof(TemporalState) +  sizeof(unsigned char) * hash->stateSize);
  res->sequenceIndex = sequenceIndex;
  ++hash->currentEntryHash;  

  if (*initialState != NULL)
  {
    (*initialState)->nextState = res;
  }
  
  //unsigned int i;
  //for (i = 0; i < hash->hashSize; ++i)
  //  Rprintf("%x",((unsigned char *)&res->initialState)[i]);

  if (dropTimeInfo)
  {  
    HASH_ADD(hh, hash->table, initialState->state, hash->stateSize, res);
  }
  else
  {
    HASH_ADD_KEYPTR(hh, hash->table, res->initialState, hash->hashSize, res);
  }
  *initialState = res;
  return false;
}

AttractorHash * addAttractorHashEntry(AttractorHashTable * hash, TemporalState * state,
                                      AttractorStateListElement * listElement)
{
  if (hash->currentEntryHash % hash->poolArraySize == 0)
  {
	    allocNewArray(&hash->hashStructs, hash->poolArraySize, 
	                  sizeof(AttractorHash));
      allocNewArray(&hash->stateStructs, hash->poolArraySize, 
	                  sizeof(TemporalState) + sizeof(unsigned char) * hash->internalStateSize);  	                  
  }
  AttractorHash * res = (AttractorHash *) (((unsigned char * )hash->hashStructs->array) + 
                                            sizeof(AttractorHash) 
                                                  * (hash->currentEntryHash % hash->poolArraySize));
  
  res->state = (TemporalState *) (((unsigned char * )hash->stateStructs->array) + (sizeof(TemporalState) +
                                   sizeof(unsigned char) * hash->internalStateSize) 
                                                  * (hash->currentEntryHash % hash->poolArraySize));
  
  memcpy(res->state, state, sizeof(TemporalState) +  sizeof(unsigned char) * hash->stateSize);
  res->attractor = listElement;
  ++hash->currentEntryHash;
  
  HASH_ADD(hh, hash->table, state->state, hash->stateSize, res);
  return res;
}

/**
 * Find the attractor belonging to state <state> in the hash table <hash>
 * and return it if it exists.
 */
AttractorStateListElement * getAttractorForState(AttractorHashTable * hash, TemporalState * state)
{
  AttractorHash * res = NULL;
 
  HASH_FIND(hh, hash->table, state->state, hash->stateSize, res);
  
  if (res == NULL)
  {
    return NULL;
  }
  else
  {
    return res->attractor;
  }
}

/**
 * Allocate a list structure that stores attractors
 */ 
AttractorStateList * allocAttractorStateList()
{
  AttractorStateList * res = CALLOC(1,sizeof(AttractorStateList));
  res->size = 0;
  res->head = res->tail = NULL;
  return res;
}

/**
 * Add an empty attractor entry to the attractor list <list>.
 * <numStates> is the number of states of the attractor.
 * Returns the empty entry.
 */
AttractorStateListElement * addAttractor(AttractorStateList * list, unsigned int numStates)
{
  AttractorStateListElement * res = CALLOC(1,sizeof(AttractorStateListElement));
  res->numStates = numStates;
  res->states = CALLOC(numStates, sizeof(TemporalState *));
  
  if (list->head == NULL)
  {
    list->head = res;
    res->index = 0;
  }
  else
  {
    list->tail->next = res;
    res->index = list->tail->index + 1;
  }
  ++list->size;
  list->tail = res;
  return res;
}

/**
 * Free all entries of the attractor list <list>
 */
void freeAttractorStateList(AttractorStateList * list)
{
  AttractorStateListElement * current = list->head, * tmp;
  
  while (true)
  {
    if (current == NULL)
      break;
      
    tmp = current;
    current = current->next;
    
    FREE(tmp->states);
    FREE(tmp);
  }
  
  FREE(list);
}

/**
 * Parse an R list <formula> specifying an expression tree and convert it to the
 * internal representation.
 * <memorySizes> is a vector specifying the size of the history for each gene.
 * <attractorSearchStartTime> is the time step at which the last predicate that
 * depends on the absolute time point has obtained its final value.
 * <geneOccurs> is a Boolean vector that is set to true if the gene occurs in the formula
 * and to false otherwise.
 * Returns an internal tree representation of <formula>.
 */
BooleanFormula * parseRTree(SEXP formula, unsigned int * memorySizes, unsigned int * attractorSearchStartTime, bool * geneOccurs)
{
  if (strcmp(CHAR(STRING_ELT(getListElement(formula, "type"), 0)), "atom") == 0)
  {
    BooleanAtom * res = calloc(1, sizeof(BooleanAtom));
    res->type = FORMULA_ATOM;
    res->negated = *LOGICAL(getListElement(formula, "negated"));
    res->literal = *INTEGER(getListElement(formula, "index")) - 1;
    res->time = - *INTEGER(getListElement(formula, "time")) - 1;
    
    if (res->literal >= 0) 
    {
      if (geneOccurs != NULL)
        geneOccurs[res->literal] = true;
        
      if (memorySizes[res->literal] < res->time + 1)
      {
        memorySizes[res->literal] = res->time + 1;
      }
    }
    
    return (BooleanFormula *)res;
  }
  else
  if (strcmp(CHAR(STRING_ELT(getListElement(formula, "type"), 0)), "constant") == 0)
  {
    Constant * res = calloc(1, sizeof(Constant));
    res->type = FORMULA_CONSTANT;
    res->negated = *LOGICAL(getListElement(formula, "negated"));
    res->value = *INTEGER(getListElement(formula, "value"));
    
    return (BooleanFormula *)res;
  }
  else
  {
    BooleanOperator * res = calloc(1, sizeof(BooleanOperator));
    res->type = FORMULA_OPERATOR;
    res->negated = *LOGICAL(getListElement(formula, "negated"));

    SEXP operands = getListElement(formula, "operands");
    res->operands = calloc(length(operands),sizeof(BooleanFormula * ));
    res->numOperands = length(operands);
    
    const char * operator = CHAR(STRING_ELT(getListElement(formula, "operator"),0));
    
    if (strcmp(operator, "|") == 0)
      res->operator = OPERATOR_OR;
    else
    if (strcmp(operator, "&") == 0)
      res->operator = OPERATOR_AND;
    else
    if (strcmp(operator, "maj") == 0)
      res->operator = OPERATOR_MAJ;
    else
    if (strcmp(operator, "sumis") == 0)
      res->operator = OPERATOR_SUMIS;      
    else
    if (strcmp(operator, "sumgt") == 0)
      res->operator = OPERATOR_SUMGT;
    else
    if (strcmp(operator, "sumlt") == 0)
      res->operator = OPERATOR_SUMLT;    
    else
    if (strcmp(operator, "timeis") == 0)
      res->operator = OPERATOR_TIMEIS;      
    else
    if (strcmp(operator, "timegt") == 0)
      res->operator = OPERATOR_TIMEGT;
    else
    if (strcmp(operator, "timelt") == 0)
      res->operator = OPERATOR_TIMELT;      
      
    unsigned int i;
    for (i = 0; i < length(operands); ++i)
    {
      res->operands[i] = parseRTree(VECTOR_ELT(operands, i), 
                                    memorySizes, 
                                    attractorSearchStartTime, 
                                    geneOccurs);
    }
    
    if (res->operator == OPERATOR_TIMEIS || res->operator == OPERATOR_TIMEGT || 
        res->operator == OPERATOR_TIMELT)
    {
      if (res->numOperands < 1 || res->operands[0]->type != FORMULA_CONSTANT)
        error("Time operator has an invalid specification!");
      
      Constant * constant = ((Constant*)res->operands[0]);
      unsigned int value = constant->value;
      
      if (res->operator == OPERATOR_TIMELT)
        --value;
      
      if (value >= *attractorSearchStartTime)
        *attractorSearchStartTime = value;
    }
    return (BooleanFormula *)res;
  }
}

/**
 * Evaluate an expression tree <formula> on the previous state <state>.
 * <stateOffsets> is an array describing where the history of the i-th gene starts.
 * <numGenes> is the total number of genes in the network.
 * Returns a logical value resulting from the evaluation.
 */ 
unsigned char evaluate(BooleanFormula * formula, TemporalState * state, 
                       unsigned int * stateOffsets,
                       const unsigned int numGenes)
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
  else
  if (formula->type == FORMULA_ATOM)
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
    
    switch(operator->operator)
    {
      case OPERATOR_OR:
        for (i = 0; i < operator->numOperands; ++i)
        {
          if (evaluate(operator->operands[i], state, stateOffsets, numGenes))
          {
              return !operator->negated;  
          }
        }
        return operator->negated;
      case OPERATOR_AND:
        for (i = 0; i < operator->numOperands; ++i)
        {
          if (!evaluate(operator->operands[i], state, stateOffsets, numGenes))
          {
              return operator->negated;  
          }
        }
        return !operator->negated;
      case OPERATOR_MAJ:
        threshold = floor(operator->numOperands/2.0);
        sum = 0;
        
        for (i = 0; i < operator->numOperands; ++i)
        {
          if (evaluate(operator->operands[i], state, stateOffsets, numGenes))
          {
            ++sum;
            if (sum > threshold)
              return !operator->negated;
          }
        }
        return operator->negated;
      case OPERATOR_SUMIS:
        threshold = ((Constant *)operator->operands[operator->numOperands-1])->value;
        sum = 0;
        
        for (i = 0; i < operator->numOperands - 1; ++i)
        {
          if (evaluate(operator->operands[i], state, stateOffsets, numGenes))
          {
            ++sum;
          }
        }        
        if (sum == threshold)
          return !operator->negated;
        else
          return operator->negated;          
      case OPERATOR_SUMGT:
        threshold = ((Constant *)operator->operands[operator->numOperands-1])->value;
        sum = 0;
        
        for (i = 0; i < operator->numOperands - 1; ++i)
        {
          if (evaluate(operator->operands[i], state, stateOffsets, numGenes))
          {
            ++sum;
            if (sum > threshold)
              return !operator->negated;
          }
        }
        return operator->negated;
      case OPERATOR_SUMLT:
        threshold = ((Constant *)operator->operands[operator->numOperands-1])->value;
        sum = 0;
        
        for (i = 0; i < operator->numOperands - 1; ++i)
        {
          if (evaluate(operator->operands[i], state, stateOffsets, numGenes))
          {
            ++sum;
            if (sum >= threshold)
              return operator->negated;
          }
        }
        return !operator->negated;
      case OPERATOR_TIMEIS:
        comparisonTimeStep = ((Constant * ) operator->operands[0])->value;
        return (state->timeStep == comparisonTimeStep - 1);
      case OPERATOR_TIMEGT:
        comparisonTimeStep = ((Constant * ) operator->operands[0])->value;
        return (state->timeStep > comparisonTimeStep - 1);
      case OPERATOR_TIMELT:
        comparisonTimeStep = ((Constant * ) operator->operands[0])->value;
        return (state->timeStep < comparisonTimeStep - 1);
      default:
        Rf_error("Unknown operator!");
    }
  }
}

/**
 * Perform a state transition of the network <network> from state <currentState>
 * at the absolute time point <timeStep>.
 * Assigns the successor state to <nextState>.
 */
void symbolicStateTransition(SymbolicBooleanNetwork * network, 
                     TemporalState * currentState, 
                     TemporalState * nextState,
                     unsigned int  * timeStep)
{
  unsigned int i;
  for (i = 0; i < network->numGenes; ++i)
  {
    if (network->stateSizes[i] > 1)
      memcpy(&nextState->state[network->stateOffsets[i] + 1], 
             &currentState->state[network->stateOffsets[i]], 
             (network->stateSizes[i] - 1) * sizeof(unsigned char));

    if (network->fixed[i] != -1)
      nextState->state[network->stateOffsets[i]] = network->fixed[i];
    else
      nextState->state[network->stateOffsets[i]] = evaluate(network->interactions[i], 
                                                     currentState, 
                                                     network->stateOffsets,
                                                     network->numGenes);
  }

  nextState->timeStep = *timeStep;
  nextState->startState = currentState->startState;    

  ++ *timeStep;
}

/**
 * R interface function to perform a state transition of 
 * the network <network> from state <previousState>
 * at the absolute time point <timeStep>.
 * Returns the successor state.
 */
SEXP symbolicStateTransition_R(SEXP network, SEXP previousState, SEXP timeStep)
{
  SymbolicBooleanNetwork * _network = R_ExternalPtrAddr(network);
  
  if (_network == NULL)
    error("Internal network structures not supplied to C handler!");
  
  unsigned int _timeStep = *INTEGER(timeStep), k;
  TemporalState * current = CALLOC(1,sizeof(TemporalState) + 
                                     _network->totalStateSize * sizeof(unsigned char));
  TemporalState * next = CALLOC(1,sizeof(TemporalState) + 
                                     _network->totalStateSize * sizeof(unsigned char));
  decodeState(_network, 
              INTEGER(previousState), 
              length(previousState), 
              current);

  current->timeStep = _timeStep;
  symbolicStateTransition(_network, current, next, &_timeStep);
 
  SEXP retSXP;
  
  PROTECT(retSXP = allocVector(INTSXP, _network->numGenes));
    
  int * resultState = INTEGER(retSXP);
      
  for (k = 0; k < _network->numGenes; ++k)
  {
     resultState[k] = 
        next->state[_network->stateOffsets[k]];
  }
  UNPROTECT(1);
  FREE(next);
  FREE(current);  
  
  return retSXP;
}

/**
 * Construct the internal expression tree representation
 * of a symbolic network specified by the R object <object>
 */
SEXP constructNetworkTrees(SEXP object)
{
  SymbolicBooleanNetwork * network = calloc(1,sizeof(SymbolicBooleanNetwork));
  SEXP interactions = getListElement(object, "interactions");
  SEXP fixed = getListElement(object, "fixed");  
    
  network->numGenes = length(interactions);
  network->attractorSearchStartTime = 0;
  network->stateSizes = calloc(network->numGenes, sizeof(unsigned int));
  
  network->stateOffsets = calloc(network->numGenes + 1, sizeof(unsigned int));
  network->fixed = calloc(network->numGenes, sizeof(int));
  network->interactions = calloc(network->numGenes, sizeof(BooleanFormula *));
  
  unsigned int i, j;
  for (i = 0; i < network->numGenes; ++i)
  {
    network->fixed[i] = INTEGER(fixed)[i];
    network->stateSizes[i] = 1;
  }
    
  for (i = 0; i < network->numGenes; ++i)
  {  
    network->interactions[i] = parseRTree(VECTOR_ELT(interactions, i), 
                                          network->stateSizes, &network->attractorSearchStartTime, NULL);
    //printFormula(network->interactions[i]);
    //Rprintf("\n");
  }

  network->totalStateSize = 0;
  for (i = 0; i < network->numGenes; ++i)
  {
    network->stateOffsets[i] = network->totalStateSize;
    network->totalStateSize += network->stateSizes[i];
  }
  network->stateOffsets[network->numGenes] = network->totalStateSize;
  network->stateFixed = calloc(network->totalStateSize, sizeof(int));    
  for (i = 0; i < network->numGenes; ++i)
  {
    for (j = 0; j < network->stateSizes[i]; ++j)
    {
      network->stateFixed[network->stateOffsets[i] + j] = network->fixed[i];
    }
  }
  
  SEXP res = PROTECT(R_MakeExternalPtr(network, install("CStructures"), R_NilValue));
  R_RegisterCFinalizerEx(res, finalizeSymbolicNetwork, true);
  UNPROTECT(1);
  return res;
}

/**
 * Convert a Boolean expression tree <tree> to a
 * truth table representation.
 * <numGenes> is the total number of genes in the network.
 * Returns a binary vector representing the truth table result column.
 */
SEXP getTruthTable(SEXP tree, SEXP numGenes)
{
  int _numGenes = *INTEGER(numGenes);
  unsigned int attractorSearchStartTime = 0;
  unsigned int * stateSizes = CALLOC(_numGenes, sizeof(unsigned int));
  unsigned int * stateOffsets = CALLOC(_numGenes, sizeof(unsigned int));  
  bool * geneOccurs = CALLOC(_numGenes, sizeof(bool));  
  
  unsigned int i, j;
  for (i = 0; i < _numGenes; ++i)
  {
    stateSizes[i] = 1;
  }

  BooleanFormula * _tree = parseRTree(tree, 
                                      stateSizes, 
                                      &attractorSearchStartTime,
                                      geneOccurs);

  if (attractorSearchStartTime > 0)
  {
    freeFormula(_tree);
    FREE(stateSizes);
    FREE(stateOffsets);
    FREE(geneOccurs);
    Rf_error("Temporal operators are not allowed in the truth table representation!");
  }
  
  unsigned int geneCount = 0;
  for (i = 0; i < _numGenes; ++i)
  {
    if (stateSizes[i] > 1)
    {
      freeFormula(_tree);
      FREE(stateSizes);
      FREE(stateOffsets);
      FREE(geneOccurs);
      Rf_error("Temporal operators are not allowed in the truth table representation!");
    }
  
    stateOffsets[i] = geneCount;
    if (geneOccurs[i])
      ++geneCount;    
  }

  SEXP res = PROTECT(allocList(2));  
  SEXP table = PROTECT(allocVector(INTSXP, 1 << geneCount));
  SEXP genes = PROTECT(allocVector(INTSXP, (geneCount == 0? 1 : geneCount)));  
  
  if (geneCount == 0)
    *INTEGER(genes) = 0;
  else
  {
    j = 0;
    for (i = 0; i < _numGenes; ++i)
    {
      if (geneOccurs[i])
      {
        INTEGER(genes)[j++] = i + 1;
      }
    }
  }
  
  TemporalState * current = CALLOC(1,sizeof(TemporalState) + 
                                    geneCount * sizeof(unsigned char));
                                      
  memset(current, 0, sizeof(TemporalState) + geneCount * sizeof(unsigned char));
  
  int * tableEntry = INTEGER(table);
  do
  {
    current->timeStep = 0;
    *tableEntry = evaluate(_tree, current, stateOffsets,
                           geneCount);
    ++tableEntry;
  }
  while (getNextState(current->state, NULL, geneCount));
  SETCAR(res, genes);
  SETCADR(res, table);  

  freeFormula(_tree);
  FREE(stateSizes);
  FREE(stateOffsets);
  FREE(geneOccurs);

  UNPROTECT(3);
  return res;
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
    
    Rprintf("var%d",atom->literal);
    
    if (atom->time != 0)
      Rprintf("[%d]", -atom->time - 1);
  }
  else
  if (formula->type == FORMULA_CONSTANT)
  {
    Constant* constant = (Constant*) formula;
    
    if (constant->negated)
      Rprintf("!");
    
    Rprintf("%d",constant->value);    
  }  
  else
  {
    BooleanOperator* operator = (BooleanOperator*) formula;
    
    if (operator->negated)
      Rprintf("!");
    
    if (operator->operator == OPERATOR_MAJ)
      Rprintf("maj");
    else
    if (operator->operator == OPERATOR_SUMGT)
      Rprintf("sumgt");
    else
    if (operator->operator == OPERATOR_TIMEIS)          
      Rprintf("timeis");
    else
    if (operator->operator == OPERATOR_TIMEGT)          
      Rprintf("timegt");
    else
    if (operator->operator == OPERATOR_TIMELT)          
      Rprintf("timelt");
      
    Rprintf("(");
    
    unsigned int i;
    for (i = 0; i < operator->numOperands; ++i)
    {
      printFormula(operator->operands[i]);
      
      if (i < operator->numOperands-1)
      {
        if (operator->operator == OPERATOR_OR)
          Rprintf(" | ");
        else
        if (operator->operator == OPERATOR_AND)        
          Rprintf(" & ");
        else
          Rprintf(", ");
      }
    }
    
    Rprintf(")");    
  }
}

/**
 * Symbolic simulator that identifies attractors and generates sequences of states
 * and attractor graphs. 
 * <network> is a pointer to the C structures of a symbolic network.
 * <states> is either a list of start states, an integer specifying the number of
 * random start states, or NULL for an exhaustive search.
 * <maxTransitions> is the maximum number of state transitions performed to find an attractor.
 * <returnSequences>, <returnGraph> and <returnAttractors> specify whether the return
 * value comprises the state sequences, the transition graph, and the attractors respectively.
 */
SEXP simulateStates(SEXP network, SEXP states, SEXP maxTransitions, 
                    SEXP returnSequences, SEXP returnGraph, SEXP returnAttractors)
{
  SymbolicBooleanNetwork * _network = R_ExternalPtrAddr(network);
  
  if (_network == NULL)
    error("Internal network structures not supplied to C handler!");
    
  unsigned int _maxTransitions = *INTEGER(maxTransitions);  

  unsigned long long _numStates = length(states);
    
  unsigned long long i, j, k;
  unsigned int mode;
  unsigned char * randomStartStates = NULL;
  if (!Rf_isNull(states))
  {
    if (IS_INTEGER(states))
    // a number of start states to generate randomly is supplied
    {
      _numStates = *INTEGER(states);
      mode = MODE_RANDOM;
      GetRNGstate();
      randomStartStates = CALLOC(_numStates, sizeof(unsigned char) * _network->totalStateSize);
    }
    else
    // the start states themselves are supplied
    {  
       _numStates = length(states);
       mode = MODE_SUPPLIED;
    }
  }
  else
  // exhaustive search
  {
    unsigned int numNonFixed = 0;
    
    for (j = 0; j < _network->numGenes; ++j)
    {
      if (_network->fixed[j] == -1)
        numNonFixed += _network->stateSizes[j];
    }
    
    _numStates = (unsigned int)1 << numNonFixed;
    mode = MODE_EXHAUSTIVE;
  }
  
  // allocate required structures
  bool _returnSequences = *LOGICAL(returnSequences);
  bool _returnGraph = *LOGICAL(returnGraph);
  bool _returnAttractors = *LOGICAL(returnAttractors);  
  
  StateHashTable * hash = allocStateHashTable(_network->totalStateSize);
  StateHashTable * usedStates = NULL;
  
   
  
  AttractorHashTable * attractorHash = allocAttractorHashTable(_network->totalStateSize);
  
  TemporalState * current = CALLOC(1,sizeof(TemporalState) + 
                                     _network->totalStateSize * sizeof(unsigned char));
  TemporalState * currentStart = CALLOC(1,sizeof(TemporalState) + 
                                     _network->totalStateSize * sizeof(unsigned char));
  TemporalState * next = CALLOC(1,sizeof(TemporalState) + 
                                     _network->totalStateSize * sizeof(unsigned char));

  AttractorStateList * attractors = allocAttractorStateList();
  unsigned long long stateCount = 0;
  unsigned long long * sequenceSizes = CALLOC(_numStates, sizeof(unsigned long long));
  int * attractorIndices = CALLOC(_numStates, sizeof(int));  
  
  if (mode == MODE_EXHAUSTIVE)
  // generate first initial state for exhaustive search (all non-fixed genes set to 0)
  {
    for (j = 0; j < _network->totalStateSize; ++j)
    {
      if (_network->stateFixed[j] == -1)
        currentStart->state[j] = 0;
      else
        currentStart->state[j] = _network->stateFixed[j];
    }
  }
  else
  if (mode == MODE_RANDOM)
  {
    usedStates = allocStateHashTable(_network->totalStateSize);
  }
  
  for (i = 0; i < _numStates; ++i)
  // iterate over start states
  { 
    R_CheckUserInterrupt();
    sequenceSizes[i] = 1;
    attractorIndices[i] = ~0;

    SEXP state;
    switch (mode)
    {
      case MODE_SUPPLIED:
        // take next start state from the input list
        state = VECTOR_ELT(states, i);
        decodeState(_network, 
                     INTEGER(state), 
                     length(state), 
                     currentStart);
        break;
      case MODE_RANDOM:
          // generate a random start state
          while (true)
          {
            StateHash * dummy = NULL;
            for (j = 0; j < _network->totalStateSize; ++j)
            {
                if (_network->stateFixed[j] == -1)
                  currentStart->state[j] = 
                    randomStartStates[i * _network->totalStateSize +  j] = intrand(2);
                else
                  currentStart->state[j] = 
                    randomStartStates[i * _network->totalStateSize +  j] = _network->stateFixed[j];
            }
            
            // check whether the state has already been choosen as a start state
            if (!findOrStore(usedStates, &dummy, currentStart, true, 0))
              break;
            
          }
          break;
    }
    
    currentStart->startState = i;
    currentStart->timeStep = 0;
    
    memcpy(current, currentStart, sizeof(TemporalState) + 
           _network->totalStateSize * sizeof(unsigned char));
    
    unsigned int timeStep = 1;
    
    StateHash * previous = NULL;
    
    if (!findOrStore(hash, &previous, current, _network->attractorSearchStartTime == 0, i))
      ++stateCount;
    
    for (j = 0; j < _maxTransitions || _maxTransitions == 0; ++j)
    // calculate successor states
    {
    
      symbolicStateTransition(_network, current, next, &timeStep);

      if (findOrStore(hash, &previous, next, _network->attractorSearchStartTime < timeStep, i))
      // a previously found state has been reached
      {
        if (next->startState == previous->initialState->startState)
        // a new attractor has been found
        {
          unsigned int attractorSize = timeStep - previous->initialState->timeStep - 1;
          AttractorStateListElement * el = addAttractor(attractors, attractorSize);
          attractorIndices[i] = el->index;
          
          for (k = 0; k < attractorSize; ++k)
          // enumerate attractor states
          {
            addAttractorHashEntry(attractorHash, previous->initialState, el);
            el->states[k] = previous->initialState;
            previous = previous->nextState;
          }

          if (attractorSize + 1 == timeStep)
          // If the whole sequence is an attractor, add the initial state at the end
            ++sequenceSizes[i];
        }
        else
        // we have reached a sequence we traversed before        
        // => determine sequence size etc. from previous run
        if (_returnSequences)
        {
          while (sequenceSizes[i] < _maxTransitions || _maxTransitions == 0)
          {
             if (previous == NULL)
             // not enough transitions)
                break;

             AttractorStateListElement * entry = getAttractorForState(attractorHash, 
                                                                      previous->initialState);
             if (entry != NULL)
             {
               sequenceSizes[i] += entry->numStates;
               attractorIndices[i] = entry->index;
               break;
             }
             ++sequenceSizes[i];
             previous = previous->nextState;             
          }
        }
        break;
      }
      
      ++sequenceSizes[i];
      ++stateCount;  
        
      memcpy(current, next, sizeof(TemporalState) + _network->totalStateSize * sizeof(unsigned char));
    }

    if (i < _numStates - 1 && mode == MODE_EXHAUSTIVE)
    // generate next state for exhaustive search
    {
          getNextState(currentStart->state, _network->stateFixed, _network->totalStateSize);
    }    
  }
 
  SEXP retSXP, attractorListSXP, initialStatesSXP, nextStatesSXP, attractorAssignmentSXP, 
       transitionSXP, iterSXP, graphSXP;
  
  PROTECT(retSXP = allocList(4));
  
  if (_returnAttractors)
  // store attractor sequences in results
  {
    PROTECT(iterSXP = attractorListSXP = allocList(attractors->size));
    AttractorStateListElement * el = attractors->head;
    for (i = 0; i < attractors->size; ++i)
    {
      SEXP attractorSXP;
      PROTECT(attractorSXP = allocVector(INTSXP, el->numStates * _network->numGenes));
      int * attractorStates = INTEGER(attractorSXP);
      for (j = 0; j < el->numStates; ++j)
      {
        for (k = 0; k < _network->numGenes; ++k)
        {
          attractorStates[j * _network->numGenes + k] = el->states[j]->state[_network->stateOffsets[k]];
        }
      }
      SETCAR(iterSXP, attractorSXP);
      UNPROTECT(1);

      iterSXP = CDR(iterSXP);
      el = el->next;
    }
   	SETCADDR(retSXP,attractorListSXP);

    if (_returnSequences)
    // also store indices of attractors for the sequences
    { 
      SEXP attractorIndexSXP;  	
      PROTECT(attractorIndexSXP = allocVector(INTSXP, _numStates));
      memcpy(INTEGER(attractorIndexSXP),attractorIndices,sizeof(int) * _numStates);
      SETCADDDR(retSXP,attractorIndexSXP);
      UNPROTECT(1);
    }
   	
   	UNPROTECT(1);
  }
  
  if (_returnSequences)
  {
    // calculate maximum delay in the network
    unsigned int maxTimeDiff = 1;
    
    for (k = 0; k < _network->numGenes; ++k)
    {
      if (_network->stateSizes[k] > maxTimeDiff)
        maxTimeDiff = _network->stateSizes[k];
    }
    
    PROTECT(iterSXP = transitionSXP = allocList(_numStates));
        
    if (mode == MODE_EXHAUSTIVE)
    // again, regenerate first start state for exhaustive search
    {
      for (j = 0; j < _network->totalStateSize; ++j)
      {
        if (_network->stateFixed[j] == -1)
          currentStart->state[j] = 0;
        else
          currentStart->state[j] = _network->stateFixed[j];
      }
    }
    
    for (i = 0; i < _numStates; ++i)
    // iterate over states
    {
      R_CheckUserInterrupt();
      SEXP transitionStatesSXP;
      PROTECT(transitionStatesSXP = allocVector(INTSXP,
                                                ((maxTimeDiff - 1) + sequenceSizes[i]) * _network->numGenes));
      int * transitionStates = INTEGER(transitionStatesSXP);
      
      SEXP state;
      switch (mode)
      {
        case MODE_SUPPLIED:
          state = VECTOR_ELT(states, i);
          decodeState(_network, 
                       INTEGER(state), 
                       length(state), 
                       currentStart);
          break;
        case MODE_RANDOM:
          memcpy(currentStart->state, &randomStartStates[i * _network->totalStateSize], 
                 sizeof(unsigned char) * _network->totalStateSize);
          break;
      }
      
      currentStart->startState = i;
      currentStart->timeStep = 0;
      
      unsigned int currentState = 0;
      for (j = 0; j < maxTimeDiff; ++j)
      {
        bool nextState = false;
        for (k = 0; k < _network->numGenes; ++k)
        {
          if (mode == MODE_SUPPLIED)
          // store supplied initial states in results
          {
              SEXP state = VECTOR_ELT(states, i);
              int * _state = INTEGER(state);
              unsigned int numStates =  length(state)/_network->numGenes;
              if ((maxTimeDiff - j) <= numStates)
              {
                transitionStates[j * _network->numGenes + k] = 
                    _state[k * numStates + (numStates - currentState - 1)];
                nextState = true;
              }
              else
                transitionStates[j * _network->numGenes + k] = _state[(k + 1) * numStates - 1];
          }
          else
          // construct list of initial states from internal vector representation          
          if (_network->stateSizes[k] >= maxTimeDiff - j)
            transitionStates[j * _network->numGenes + k] = 
                currentStart->state[_network->stateOffsets[k] + (maxTimeDiff - j - 1)];
          else
            transitionStates[j * _network->numGenes + k] = 
                currentStart->state[_network->stateOffsets[k] + _network->stateSizes[k] - 1];
                             
        }
        
        if (nextState)
        {
          ++currentState;
          nextState = false;
        }
      }
    
      memcpy(current, currentStart, sizeof(TemporalState) + 
             _network->totalStateSize * sizeof(unsigned char));
             
      // get successor states from hash table
      StateHash * res;
      if (_network->attractorSearchStartTime == 0)    
        HASH_FIND(hh, hash->table, current->state, hash->stateSize, res);
      else
        HASH_FIND(hh, hash->table, current, hash->hashSize, res);
      
      for (j = 0; j < sequenceSizes[i]; ++j)
      {
        if (res != NULL)
        // encode successor state for results
        {
          for (k = 0; k < _network->numGenes; ++k)
          {
            transitionStates[(maxTimeDiff + j - 1) * _network->numGenes + k] = res->initialState->state[_network->stateOffsets[k]];
          }
          res = res->nextState;
        }
        else
          error("Did not find a required successor state!");
      }
      SETCAR(iterSXP, transitionStatesSXP);
      iterSXP = CDR(iterSXP);
      UNPROTECT(1);
      
      if (mode == MODE_EXHAUSTIVE)
      // generate next start state for exhaustive search
      {
        getNextState(currentStart->state, _network->stateFixed, _network->totalStateSize);
      }    
    }

    SETCAR(retSXP,transitionSXP);
    UNPROTECT(1);
  }
  
  if (mode == MODE_RANDOM)
  {
    FREE(randomStartStates);  
    PutRNGstate();
  }

  StateHash * entry, * tmp;  
  if (_returnGraph)
  // generate the graph representation
  {
    PROTECT(graphSXP = allocList(3));  
    PROTECT(initialStatesSXP = allocVector(INTSXP,stateCount * _network->numGenes));
    PROTECT(nextStatesSXP = allocVector(INTSXP,stateCount * _network->numGenes));
    PROTECT(attractorAssignmentSXP = allocVector(INTSXP,stateCount));    
  
    int * initialStates = INTEGER(initialStatesSXP);
    int * nextStates = INTEGER(nextStatesSXP);  
    int * attractorAssignment = INTEGER(attractorAssignmentSXP);    
    
    // add all state pairs in the hash table to the graph
    i = 0;
    HASH_ITER(hh, hash->table, entry, tmp) 
    {
      R_CheckUserInterrupt();
      for (k = 0; k < _network->numGenes; ++k)
      {
        initialStates[i * _network->numGenes + k] = entry->initialState->state[_network->stateOffsets[k]];
        
        if (entry->nextState != NULL)
          nextStates[i * _network->numGenes + k] = entry->nextState->initialState->state[_network->stateOffsets[k]];
        else
          nextStates[i * _network->numGenes + k] = entry->initialState->state[_network->stateOffsets[k]];
        
      }
      attractorAssignment[i] = attractorIndices[entry->sequenceIndex];
      ++i;
    }
   	SETCAR(graphSXP,initialStatesSXP);
  	SETCADR(graphSXP,nextStatesSXP);
  	SETCADDR(graphSXP,attractorAssignmentSXP);   	  	   	
   	SETCADR(retSXP,graphSXP);
   	UNPROTECT(4);
  }
  
  UNPROTECT(1);  
  
  
  if (mode == MODE_RANDOM)
    freeStateHashTable(usedStates);
  
  freeStateHashTable(hash);
  freeAttractorHashTable(attractorHash);
    
  FREE(sequenceSizes);
  FREE(next);
  FREE(current);  
  FREE(currentStart);
  FREE(attractorIndices); 
 
  freeAttractorStateList(attractors);
  
  return retSXP;
}

/**
 * Check whether <ptr> is a null pointer.
 */
SEXP checkNullPointer(SEXP ptr)
{
  void * p = R_ExternalPtrAddr(ptr);
  
  SEXP ret;
  PROTECT(ret = allocVector(LGLSXP,1));
   *LOGICAL(ret) = (p == NULL);
  UNPROTECT(1);
  
  return(ret);
}

