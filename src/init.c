#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void bin2decC(void *, void *, void *);
extern void dec2binC(void *, void *, void *);
extern void freeAllMemory(void);

/* .Call calls */
extern SEXP checkNullPointerC(SEXP);
extern SEXP constructNetworkTrees_R(SEXP);
extern SEXP getAttractors_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getTruthTable_R(SEXP, SEXP);
extern SEXP markovSimulation_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP reconstructNetwork_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simulateStates_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP symbolicSATSearch_R(SEXP, SEXP, SEXP);
extern SEXP symbolicStateTransition_R(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"bin2decC",       (DL_FUNC) &bin2decC,       3},
    {"dec2binC",       (DL_FUNC) &dec2binC,       3},
    {"freeAllMemory", (DL_FUNC) &freeAllMemory, 0},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"checkNullPointerC",          (DL_FUNC) &checkNullPointerC,           1},
    {"constructNetworkTrees_R",   (DL_FUNC) &constructNetworkTrees_R,    1},
    {"getAttractors_R",           (DL_FUNC) &getAttractors_R,           12},
    {"getTruthTable_R",           (DL_FUNC) &getTruthTable_R,            2},
    {"markovSimulation_R",        (DL_FUNC) &markovSimulation_R,        11},
    {"reconstructNetwork_R",      (DL_FUNC) &reconstructNetwork_R,      10},
    {"simulateStates_R",          (DL_FUNC) &simulateStates_R,           6},
    {"symbolicSATSearch_R",       (DL_FUNC) &symbolicSATSearch_R,        3},
    {"symbolicStateTransition_R", (DL_FUNC) &symbolicStateTransition_R,  3},
    {NULL, NULL, 0}
};

void R_init_BoolNet(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
