#ifndef GPU_BETWEENNESS_CENTRALITY_H
#define GPU_BETWEENNESS_CENTRALITY_H

#include "Configuration.h"
#include "Graph.h"


/**
 * Computes the Betweenness Centrality metric.
 */
void computeBetweennessCentrality(Configuration *config, Graph *g, double *bc);


#endif // gpu/BetweennessCentrality.h
