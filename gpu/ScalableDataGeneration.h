#ifndef GPU_SCALABLE_DATA_GENERATION_H
#define GPU_SCALABLE_DATA_GENERATION_H

#include <stdio.h>

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
void generateScalableData(Configuration *config, int *permV, GraphSDG *tuples);

/**
 * Generates a permutation of vertices.
 */
void generatePermutation(int n, int *permutation);

/**
 * Prints generated tuples.
 */
void printTuples(FILE *stream, GraphSDG *tuples); 

#endif // gpu/ScalableDataGeneration.h
