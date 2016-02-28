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
  if (argc != 3)
  {
    fprintf(stderr, "Usage: ./SSCA2 <SCALE> <GPU>\n");
    return EXIT_FAILURE;
  }

  // TODO: Refactor.
  // Determine whether to run CPU or GPU version.
  bool runGPU = false;
  if (atoi(argv[2]) == 1)
  {
    runGPU = true; 
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
  Graph *graph;

  // Allocate in unified memory if running on the GPU.
  if (runGPU)
  {
    allocateManaged(&config, &graph);
  }
  else
  {
    allocateHost(&config, &graph);
  }

  // Construct the graph from the tuples.
  constructGraph(&tuples, graph);

  elapsedTime = getSeconds() - elapsedTime;
  fprintf(stderr, "Time taken for Kernel 1 is %9.6lf sec.\n", elapsedTime);

  // Clean up the memory used to store generated data.
  destroy(&tuples);

  /* ---------------------------------------- */
  /* Kernel 4 - Betweenness Centrality        */
  /* ---------------------------------------- */
  fprintf(stderr, "Kernel 4: Computing Betweenness Centrality...\n");
  
  // Allocating bc array and generating vertex permutations is not timed.
  float *bc;
  int *perm;

  if (runGPU)
  {
    cudaMallocManaged(&bc, config.n * sizeof(float));
    cudaMemset(bc, 0.0, config.n);
    cudaMallocManaged(&perm, config.n * sizeof(int));
  }
  else
  {
    bc = (float *) calloc(config.n, sizeof(double)); // float!
    assert(bc != NULL);
    perm = (int *) malloc(config.n * sizeof(int));
    assert(perm != NULL);
  }

  // Permute the vertices.
  generatePermutation(config.n, perm);

  // Start timing.
  elapsedTime = getSeconds();
  
  // Compute the betweenes centrality metric.
  if (runGPU)
  {
    computeBCGPU(&config, graph, perm, bc);
  }
  else
  {
    computeBCCPU(&config, graph, perm, bc);
  }
  
  elapsedTime = getSeconds() - elapsedTime;
  fprintf(stderr, "Time taken for Kernel 4 is %9.6lf sec.\n", elapsedTime);

  fprintf(
      stderr,
      "TEPS score for Kernel 4 is %lf\n",
      7 * config.n * ((long) (1 << config.k4Approx)) / elapsedTime);

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
  if (runGPU)
  {
    cudaFree(perm);
    cudaFree(bc);
    destroyManaged(&graph);
  }
  else
  {
    free(perm);
    free(bc);
    destroyHost(&graph);
  }
}

