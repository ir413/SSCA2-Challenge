#ifndef GPU_BETWEENNESS_CENTRALITY_H
#define GPU_BETWEENNESS_CENTRALITY_H

#include "Configuration.h"
#include "Graph.h"


typedef struct
{
  int *list;
  int count;
  int degree;
} plist;


/**
 * Computes the Betweenness Centrality metric.
 */
void computeBetweennessCentrality(Configuration *config, Graph *g, double *bc);

/**
 * Computes the Betweenness Centrality on the CPU.
 */
void computeBCCPU(Configuration *config, Graph *g, int *perm, double *bc);


#endif // gpu/BetweennessCentrality.h
