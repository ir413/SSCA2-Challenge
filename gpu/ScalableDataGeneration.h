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
} TuplesSDG;


/**
 * Allocates memory for the tuples.
 */
void allocate(Configuration *config, TuplesSDG *tuples);

/**
 * Frees the memory occupied by the tuples.
 */
void destroy(TuplesSDG *tuples);

/**
 * Generates tuples. 
 */
void generateScalableData(Configuration *config, TuplesSDG *tuples);

/**
 * Generates a permutation of vertices.
 */
void generatePermutation(int n, int *permutation);

/**
 * Prints generated tuples.
 */
void printTuples(FILE *stream, TuplesSDG *tuples); 

#endif // gpu/ScalableDataGeneration.h
