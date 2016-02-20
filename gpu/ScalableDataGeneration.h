#ifndef GPU_SCALABLE_DATA_GENERATION_H
#define GPU_SCALABLE_DATA_GENERATION_H

#include "Configuration.h"


/**
 * Represents tuples of scalably generated data.
 */
typedef struct
{
  // Number of vertices.
  int n;
  // Number of edges.
  int m;
  // Edge lists.
  int *startVertex;
  int *endVertex;
  int *weight;
} GraphSDG;

/**
 * Generates tuples. 
 */
double generateScalableData(Configuration *config, int *permV, GraphSDG *tuples);

#endif // gpu/ScalableDataGeneration.h
