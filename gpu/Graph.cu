#include <assert.h>
#include <stdio.h>

#include "PrefixSums.h"

#include "Graph.h"

// Number of thread bloks.
#define THREAD_BLOCKS_COUNT 1
// Number of threads per block.
#define THREADS_PER_BLOCK 1024 


// TODO: Parallelize
void constructGraph(GraphSDG *tuples, Graph *graph)
{
  graph->n = tuples->n;
  graph->m = tuples->m;

  // Allocate memory for temporary arrays.
  int *pos = (int *) malloc(tuples->m * sizeof(int));
  assert( pos != NULL); 
  int *degree = (int *) calloc(graph->n, sizeof(int));
  assert(degree != NULL);

  // Compute the row offsets. 
  for (int i = 0; i < tuples->m; ++i)
  {
    int u = tuples->startVertex[i];
    pos[i] = degree[u]++;
  } 

  //TODO: use prefixSums(degree, graph->rowOffset, graph->n);
  graph->rowOffset[0] = 0;
  for (int i = 1; i < (graph->n + 1); ++i)
  {
    graph->rowOffset[i] = graph->rowOffset[i - 1] + degree[i - 1];
  }

  // Compute the adjecency lists and weights.
  for (int i = 0; i < tuples->m; ++i)
  {
    // Compute the index in the adjecency and weight lists.
    int u = tuples->startVertex[i];
    int columnIndex = graph->rowOffset[u] + pos[i];
    // Update the lists.
    graph->column[columnIndex] = tuples->endVertex[i];
    graph->weight[columnIndex] = tuples->weight[i];
  }

  // Clean up the temporary storage.
  free(degree);
  free(pos);
}

void printGraph(FILE *stream, Graph *graph)
{
  fprintf(stream, "Row offsets:\n");
  for (int i = 0; i < (graph->n + 1); ++i)
  {
    fprintf(stream, "%d ", graph->rowOffset[i]);
  }

  fprintf(stream, "\nAdjecency lists:\n");
  for (int i = 0; i < graph->m; ++i)
  {
    fprintf(stream, "%d ", graph->column[i]);
  }

  fprintf(stream, "\nWeights:\n");
  for (int i = 0; i < graph->m; ++i)
  {
    fprintf(stream, "%d ", graph->weight[i]);
  }
  fprintf(stream, "\n");
}

