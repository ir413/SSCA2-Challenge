#ifndef GPU_GRAPH_H
#define GPU_GRAPH_H

#include <stdio.h>

#include "Configuration.h"
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


/**
 * Allocates memory for the graph on the host.
 */
void allocateHost(Configuration *config, Graph **graph);

/**
 * Allocates memory for the graph in the unified memory.
 */
void allocateManaged(Configuration *config, Graph **graph);

/**
 * Allocates memory for the graph.
 */
void allocateManaged(Configuration *config, Graph **graph);

/**
 * Frees the host memory occupied by the graph.
 */
void destroyHost(Graph **graph);

/**
 * Frees the unified memory occupied by the graph.
 */
void destroyManaged(Graph **graph);

/**
 * Consturcts a graph from a list of tuples.
 */ 
void constructGraph(TuplesSDG *tuples, Graph *graph);

/**
 * Prints the graph.
 */
void printGraph(FILE *stream, Graph *graph);


#endif // gpu/Graph.h
