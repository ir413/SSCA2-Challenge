#ifndef GPU_BETWEENNESS_CENTRALITY_H
#define GPU_BETWEENNESS_CENTRALITY_H

#include "Configuration.h"
#include "Graph.h"


typedef struct
{
  int *list;
  unsigned int count;
} plist;


/**
 * Computes the Betweenness Centrality metric on the GPU.
 */
void computeBCGPU(Configuration *config, Graph *g, int *perm,  float *bc);

/**
 * Computes the Betweenness Centrality metric on the CPU.
 */
void computeBCCPU(Configuration *config, Graph *g, int *perm, float *bc);


#endif // gpu/BetweennessCentrality.h
