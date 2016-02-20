#ifndef GPU_SCALABLE_DATA_GENERATION_H
#define GPU_SCALABLE_DATA_GENERATION_H


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
} graphSDG;

/**
 * Generates tuples. 
 */
double generateScalableData(graphSDG *tuples);

#endif // gpu/ScalableDataGeneration.h
