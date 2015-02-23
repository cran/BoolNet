#include "attractor_info.h"
#include "common.h"

/**
 * Allocate a new AttractorInfo structure for <tableSize> states
 */
pAttractorInfo allocAttractorInfo(unsigned long long tableSize, unsigned int numGenes)
{
	pAttractorInfo res = (pAttractorInfo)CALLOC(1,sizeof(AttractorInfo));
	if ((numGenes % BITS_PER_BLOCK_32) == 0)
		res->numElementsPerEntry = numGenes/BITS_PER_BLOCK_32;
	else
		res->numElementsPerEntry = numGenes/BITS_PER_BLOCK_32 + 1;
	res->numAttractors = 0;
	res->table = NULL;
	res->tableSize = tableSize;
	res->initialStates = NULL;
	res->table = (unsigned int*) CALLOC(tableSize * res->numElementsPerEntry,sizeof(unsigned int));
	res->attractorAssignment = (unsigned int*) CALLOC(tableSize,sizeof(unsigned int));
	res->stepsToAttractor = (unsigned int*) CALLOC(tableSize,sizeof(unsigned int));
	return res;
}

/**
 * Free a list of attractor structures
 */
void freeAttractorList(pAttractor p)
{
	do
	{
		pAttractor next = p->next;
		if (p->involvedStates != NULL)
  		FREE(p->involvedStates);
		FREE(p);
		p = next;
	}
	while(p != NULL);
}

/**
 * Free an AttractorInfo structure including
 * all sub-elements and the attractor list
 */
void freeAttractorInfo(pAttractorInfo p)
{
	if (p->initialStates != 0)
		FREE(p->initialStates);
	FREE(p->table);
	FREE(p->attractorAssignment);
	FREE(p->stepsToAttractor);
	freeAttractorList(p->attractorList);
	FREE(p);
}
