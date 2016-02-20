#ifndef GPU_CONFIGURATION_H
#define GPU_CONFIGURATION_H

/**
 * Benchmark configuration.
 * Constants used for data generation and subsequent kernels.
 */
struct Configuration {
  // Binary scaling heuristic parameter.
  int scale;
  // Number of vertices.
  int n;
  // Number of edges.
  int m;
  // Maximum integer weight.
  int maxIntWeight;
  // Maximum path length in subgraphs generated in Kernel 3.
  int subGraphPathLength;
  // Controls exactness of BC computation.
  int k4Approx;
  // R-MAT graph data generation parameters.
  double a;
  double b;
  double c;
  double d; 
};

void configure(int scale, struct Configuration *c);

#endif // gpu/Configuration.h
