#include <assert.h>
#include <stdio.h>

#include "Configuration.h"
#include "Graph.h"
#include "ScalableDataGeneration.h"
#include "Timer.h"


/**
 * Entry point.
 */
int main(int argc, char **argv)
{
  if (argc != 2)
  {
    fprintf(stderr, "Usage: ./SSCA2 <SCALE>\n");
    return EXIT_FAILURE;
  }

  double elapsedTime;

  /* ----------------------------------------- */
  /* Initialization -- Untimed                 */
  /* ----------------------------------------- */
  fprintf(stderr, "Configuring the Benchmark...\n");

  Configuration config;
  configure(atoi(argv[1]), &config);
  fprintf(stderr, "N: %d M: %d\n", config.n, config.m);

  /* ----------------------------------------- */
  /* Scalable Data Generation -- Untimed       */
  /* ------------------------------------------*/
  fprintf(stderr, "Scalable Data Generation...\n");
  elapsedTime = getSeconds();

  // Consturct the tuples structure.
  GraphSDG tuples;
  // Allocate memory required for the tuple arrays.
  tuples.startVertex = (int *) malloc(config.m * sizeof(int));
  tuples.endVertex = (int *) malloc(config.m * sizeof(int));
  tuples.weight = (int *) malloc(config.m * sizeof(int));
  // Allocate memory for the temporary permV array. 
  int *permV = (int *) malloc(config.m * sizeof(int));

  // Consturct the tuples.
  generateScalableData(&config, permV, &tuples);

  // Free the memory used for the temporary permV array.
  free(permV);

  elapsedTime = getSeconds() - elapsedTime;
  fprintf(
      stderr,
      "Time taken for Scalable Data Generation is %9.6lf sec.\n",
      elapsedTime);

  /* ----------------------------------------- */
  /* Kernel 1 - Graph Construction             */
  /* ----------------------------------------- */
  fprintf(stderr, "Constructing the graph...\n");
  elapsedTime = getSeconds();

  // Consturct the graph structure. 
  Graph graph;
  // Allocate memory required for the graph.
  graph.rowOffset = (int *) malloc((config.n + 1) * sizeof(int));
  assert(graph.rowOffset != NULL);
  graph.column = (int *) malloc(config.m * sizeof(int));
  assert(graph.column != NULL);
  graph.weight = (int *) malloc(config.m * sizeof(int)); 
  assert(graph.weight != NULL);

  // Construct the graph.
  constructGraph(&tuples, &graph);

  elapsedTime = getSeconds() - elapsedTime;
  fprintf(stderr, "Time taken for Kernel 1 is %9.6lf sec.\n", elapsedTime);

  // Clean up the memory used to store generated data.
  free(tuples.weight);
  free(tuples.endVertex);
  free(tuples.startVertex);

  /* ---------------------------------------- */
  /* Kernel 2 - Find max edge weight          */
  /* ---------------------------------------- */
  // TODO

  // Clean up the memory used to store the graph.
  free(graph.weight);
  free(graph.column);
  free(graph.rowOffset);
}

