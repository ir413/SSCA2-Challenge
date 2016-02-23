#include <assert.h>

#include "Graph.h"


void allocate(Configuration *config, Graph *graph)
{
  assert(config != NULL);
  assert(graph != NULL);

  graph->n = config->n;
  graph->m = config->m;

  graph->rowOffset = (int *) malloc((graph->n + 1) * sizeof(int));
  assert(graph->rowOffset != NULL);

  graph->column = (int *) malloc(graph->m * sizeof(int));
  assert(graph->column != NULL);

  graph->weight = (int *) malloc(graph->m * sizeof(int));
  assert(graph->weight != NULL);
}

void destroy(Graph *graph)
{
  assert(graph != NULL);

  free(graph->weight);
  free(graph->column);
  free(graph->rowOffset);
}

void constructGraph(GraphSDG *tuples, Graph *graph)
{
  assert(tuples != NULL);
  assert(graph != NULL);

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
  assert(stream != NULL);
  assert(graph != NULL);

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

