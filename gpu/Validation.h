#ifndef GPU_VALIDATION_H
#define GPU_VALIDATION_H

#include "Configuration.h"


/**
 * Determines if the computed Betweenness Centrality scores are
 * valid for the given benchmark configuration.
 */
bool isValid(Configuration *config, float *bc);


#endif // gpu/Validation.h
