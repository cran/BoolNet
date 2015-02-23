#ifndef SAT_SEARCH_H
#define SAT_SEARCH_H

#include "boolean_network.h"
#include "attractor_info.h"

#define EXTENSION_EXPONENTIAL 0
#define EXTENSION_LINEAR 1
#define EXTENSION_LINEAR_ADAPT 2
#define EXTENSION_MIXED 3

/**
 * Determine all attractors having at most <maxLength> states in the network <network>,
 * and return them in an AttractorInfo structure.
 */
extern pAttractorInfo getAttractors_SAT_exhaustive(BooleanNetwork * network,
                                                   unsigned int initialCycleSearchLength,
                                                   unsigned int extensionMode);

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
extern pAttractorInfo getAttractors_SAT_maxLength(BooleanNetwork * network,
                                                  unsigned int maxLength);

#endif
