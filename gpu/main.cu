#include <assert.h>
#include <stdio.h>

#include "BetweennessCentrality.h"
#include "Configuration.h"
#include "Graph.h"
#include "ScalableDataGeneration.h"
#include "Timer.h"
#include "Validation.h"


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

  // Allocate memory required for the tuples.
  TuplesSDG tuples;
  allocate(&config, &tuples);

  // Generate the data.
  generateScalableData(&config, &tuples);

  elapsedTime = getSeconds() - elapsedTime;
  fprintf(
      stderr,
      "Time taken for Scalable Data Generation is %9.6lf sec.\n",
      elapsedTime);

  /* ----------------------------------------- */
  /* Kernel 1 - Graph Construction             */
  /* ----------------------------------------- */
  fprintf(stderr, "Kernel 1: Constructing the graph...\n");
  elapsedTime = getSeconds();

  // Allocate memory required for the graph.
  Graph graph;
  allocate(&config, &graph);

  // Construct the graph from the tuples.
  constructGraph(&tuples, &graph);

  elapsedTime = getSeconds() - elapsedTime;
  fprintf(stderr, "Time taken for Kernel 1 is %9.6lf sec.\n", elapsedTime);

  // Clean up the memory used to store generated data.
  destroy(&tuples);

  /* ---------------------------------------- */
  /* Kernel 2 - Find max edge weight          */
  /* ---------------------------------------- */
  // TODO

  /* ---------------------------------------- */
  /* Kernel 3 - Graph Extraction              */
  /* ---------------------------------------- */
  // TODO

  /* ---------------------------------------- */
  /* Kernel 4 - Betweenness Centrality        */
  /* ---------------------------------------- */
  fprintf(stderr, "Kernel 4: Computing Betweenness Centrality...\n");
  
  // Allocating bc array and generating vertex permutations is not timed.
  double *bc = (double *) calloc(config.n, sizeof(double));
  assert(bc != NULL);

  int *perm = (int *) malloc(config.n * sizeof(int));
  assert(perm != NULL);
  generatePermutation(config.n, perm);

  // Start timing.
  elapsedTime = getSeconds();
  
  // Compute the betweenes centrality metric.
  computeBCCPU(&config, &graph, perm, bc);
  
  elapsedTime = getSeconds() - elapsedTime;
  fprintf(stderr, "Time taken for Kernel 4 is %9.6lf sec.\n", elapsedTime);

  /* ---------------------------------------- */
  /* Validation                               */
  /* ---------------------------------------- */
  fprintf(stderr, "Validating the results...\n");

  if (isValid(&config, bc))
  {
    fprintf(stderr, "Kernel 4 validation successful!\n");
  }
  else
  {
    fprintf(stderr, "Kernel 4 validation failed!\n");
  }

  // Clean up.
  free(perm);
  free(bc);

  destroy(&graph);
}

