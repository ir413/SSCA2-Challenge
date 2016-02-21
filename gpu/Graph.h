#ifndef GPU_GRAPH_H
#define GPU_GRAPH_H

#include "ScalableDataGeneration.h"


/**
 * Graph representation.
 * The grpah is stored in the Compressed Sparse Row (CSR) format.
 *
 * The representation consists of three arrays:
 * rowOffset - Points at the start of each vertex's
 *             adjecency list in the column and weight arrays.
 * column    - Concatenation of each vertex's adjecency list.
 * weight    - Weights corresponding to the column array. 
 */
typedef struct
{
  // Number of vertices.
  int n;
  // Number of edges.
  int m; 
  // Row offsets.
  int *rowOffset;
  // Adjecency lists.
  int *column;
  // Weights/
  int *weight; 
} Graph;

void constructGraph(GraphSDG *tuples, Graph *graph);

#endif // gpu/Graph.h
