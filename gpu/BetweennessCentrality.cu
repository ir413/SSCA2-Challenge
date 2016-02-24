#include <assert.h>
#include <stdio.h>

#include "BetweennessCentrality.h"

// Number of thread blocks.
#define BLOCKS_COUNT 1
// Maximal number of threads per block.
#define MAX_THREADS_PER_BLOCK 1024 

__global__ void initialize(
    int source,
    int n,
    int *d,
    int *sigma,
    double *delta,
    plist *p)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  if (idx < n)
  {
    if (idx == source)
    {
      d[idx] = 0;
      sigma[idx] = 1;
    }
    else
    {
      d[idx] = -1; 
      sigma[idx] = 0;
      delta[idx] = 0.0;
      p[idx].count = 0;
    }
  }
}

__global__ void vertexParallelBFS(
    int source,
    Graph *g,
    int *d,
    int *sigma,
    plist *p,
    int *stack,
    int *sPtr)
{
  // We use empty and level to implicitly represent the standard bfs queue.
  __shared__ bool empty;
  __shared__ int level;

  if (threadIdx.x == 0)
  {
    empty = false;
    level = 0;
  }

  __syncthreads();

  while (!empty)
  {
    __syncthreads();
    empty = true;
    __syncthreads();

    // TODO: unroll the first iteration.
    for (int v = threadIdx.x; v < g->n; v += blockDim.x)
    {
      // Check if vs neighbours are to be visited i.e. if v is in the queue.
      if (d[v] == level)
      {
        // Push v onto the stack.
        stack[level] = v;

        // Go through the successors of v.
        for (int j = g->rowOffset[v]; j < g->rowOffset[v + 1]; ++j)
        {
          // Skip edges whose weight is divisible by 8.
          if ((g->weight[j] & 7) == 0)
          {
            continue;
          }

          int w = g->column[j];

          // Not visited.
          if (d[w] == -1)
          {
            // Mark that there are more nodes to be explored.
            empty = false; 
            // No need for an atomic update.
            d[w] = d[v] + 1;
            // Set sigma.
            atomicExch(&sigma[w], sigma[v]);
            // Save the predecessor.
            p[w].list[p[w].count] = v;
            atomicInc(&(p[w].count), 1);
          }
          else if (d[w] == d[v] + 1)
          {
            atomicAdd(&sigma[w], sigma[v]);
            // Save the predecessor.
            p[w].list[p[w].count] = v;
            atomicInc(&(p[w].count), 1);
          }
        }
      }
    }

    if (threadIdx.x == 0)
    {
      level++;
    }

    __syncthreads();
  }

  if (threadIdx.x == 0)
  {
    *sPtr = level;
  }
}

__global__ void accumBC(
    int source,
    Graph *g,
    int *sigma,
    double *delta,
    plist *p,
    int *stack,
    double *bc)
{

}

void computeBCGPU(Configuration *config, Graph *g, int *perm, double *bc)
{
  // Declare the auxilary structures.
  int *d;
  int *sigma;
  double *delta;
  int *stack;
  plist *p;
  int *pListMem;
  int *sPtr;

  // Allocate temporary structures in global memory.
  cudaMallocManaged(&d, g->n * sizeof(int));
  cudaMallocManaged(&sigma, g->n * sizeof(int));
  cudaMallocManaged(&delta, g->n * sizeof(double));
  cudaMallocManaged(&stack, g->n * sizeof(int));
  cudaMallocManaged(&p, g->n * sizeof(plist));
  cudaMallocManaged(&pListMem, g->m * sizeof(int));
  cudaMallocManaged(&sPtr, sizeof(int));

  // --- TMP
  int *inDegree;
  int *numEdges;

  cudaMallocManaged(&inDegree, (g->n + 1) * sizeof(int));
  cudaMemset(inDegree, 0, g->n + 1);
  cudaMallocManaged(&numEdges, (g->n + 1) * sizeof(int));

  for (int i = 0; i < g->m; ++i)
  {
    inDegree[g->column[i]]++;
  }

  numEdges[0] = 0;
  for (int i = 1; i < (g->n + 1); ++i)
  {
    numEdges[i] = numEdges[i - 1] + inDegree[i - 1];
  }

  for (int i = 0; i < g->n; ++i)
  {
    p[i].list = pListMem + numEdges[i];
    p[i].degree = inDegree[i];
    p[i].count = 0;
  }

  cudaFree(inDegree);
  cudaFree(numEdges);
  // --- TMP

  for (int source = 0; source < g->n; ++source)
  {
    // Initialize the data structures.
    initialize<<<BLOCKS_COUNT, MAX_THREADS_PER_BLOCK>>>(
        source,
        g->n,
        d,
        sigma,
        delta,
        p);
    cudaDeviceSynchronize();

    // Run BFS.
    vertexParallelBFS<<<BLOCKS_COUNT, MAX_THREADS_PER_BLOCK>>>(
        source,
        g,
        d,
        sigma,
        p,
        stack,
        sPtr);
    cudaDeviceSynchronize();

    // Sum centrality scores.
    /*
    accumBC<<<BLOCKS_COUNT, MAX_THREADS_PER_BLOCK>>>(
        source,
        g,
        sigma,
        delta,
        p,
        stack,
        bc);
    cudaDeviceSynchronize();
    */

    // --- TMP
    while (*sPtr > 0)
    {
      (*sPtr)--;
      int w = stack[*sPtr];

      for (int k = 0; k < p[w].count; ++k)
      {
        int v = p[w].list[k];
        delta[v] = delta[v] + (
            ((double) sigma[v]) / ((double) sigma[w])) * (1.0 + delta[w]);
      }

      if (w != source)
      {
        bc[w] += delta[w];
      }
    }
    // --- TMP
  }

  // Clean up.
  cudaFree(sPtr);
  cudaFree(pListMem);
  cudaFree(p);
  cudaFree(stack);
  cudaFree(delta);
  cudaFree(sigma);
  cudaFree(d);
}

void computeBCCPU(Configuration *config, Graph *g, int *perm, double *bc)
{
  assert(config != NULL);
  assert(g != NULL);
  assert(perm != NULL);
  assert(bc != NULL);

  // Number of vertices to run BFS from.
  int rootsCount = 1 << config->k4Approx;

  // Predecessors of a vertex v on shortest paths from s.
  plist *p;
  int *pListMem;

  // Number of shortest paths.
  double *sigma;
  // Shortest path length between the pairs.
  int *d;
  // Dependency of vertices.
  double *delta;
 
  // inDegree array wasteful -> can use p[].degree 
  int *inDegree;
  int *numEdges;

  // The size of the predecessor list of each vertex is bounded
  // by its in-degree -> Compute the in-degree of every vertex.
  p = (plist *) calloc(g->n, sizeof(plist));

  inDegree = (int *) calloc(g->n + 1, sizeof(int));
  numEdges = (int *) malloc((g->n + 1) * sizeof(int));
  
  // Compute inDegree.
  for (int i = 0; i < g->m; ++i)
  {
    int v = g->column[i];
    inDegree[v]++;
  }

  // Prefix sums.
  numEdges[0] = 0;
  for (int i = 1; i < (g->n + 1); ++i)
  {
    numEdges[i] = numEdges[i - 1] + inDegree[i - 1];
  }

  // Allocate memory for plists.
  pListMem = (int *) malloc(g->m * sizeof(int));

  // Set the p pointers accordingly.
  for (int i = 0; i < g->n; ++i)
  {
    p[i].list = pListMem + numEdges[i];
    p[i].degree = inDegree[i];
    p[i].count = 0;
  }

  // Clean up temporary structures.
  free(inDegree);
  free(numEdges);

  // Allocate memory.
  sigma = (double *) malloc(g->n * sizeof(double));
  d = (int *) malloc(g->n * sizeof(int));
  delta = (double *) calloc(g->n, sizeof(double));

  int *stack = (int *) malloc(g->n * sizeof(int));
  int *queue = (int *) malloc(g->n * sizeof(int));

  for (int r = 0; r < g->n; ++r)
  {
    // Check if the required number of roots has been explored.
    if (rootsCount == 0)
    {
      break;
    }

    // Apply the permutation.
    int root = perm[r];

    // Skip vertices with no outgoing edges.
    if (g->rowOffset[root + 1] - g->rowOffset[root] == 0)
    {
      continue;
    }
    else
    {
      rootsCount--;
    }

    int sPtr = 0;
    int qHead = 0;
    int qTail = 0;

    // Clear the auxilary structures.
    for (int k = 0; k < g->n; ++k)
    {
      d[k] = -1;
      sigma[k] = 0.0;
      p[k].count = 0;
      delta[k] = 0.0;
    }

    sigma[root] = 1.0;
    d[root] = 0;

    queue[qTail] = root;
    qTail++;

    // While !empty(Q)
    while (qTail - qHead > 0)
    {
      // Dequeue v <- Q
      int v = queue[qHead];
      qHead++;
      // Push v -> Stack
      stack[sPtr] = v;
      sPtr++; 

      for (int j = g->rowOffset[v]; j < g->rowOffset[v + 1]; ++j)
      {
        // Skip edges whose weight is divisible by 8.
        if ((g->weight[j] & 7) == 0)
        {
          continue;
        }

        int w = g->column[j];

        if (v == w) {
          continue;
        }

        // Not visited.
        if (d[w] == -1)
        {
          // Enqueue w -> Q
          queue[qTail] = w;
          qTail++;
        
          // Distance to w is distance to v + 1.
          d[w] = d[v] + 1;

          sigma[w] = sigma[v];
          // Save the predecesor.
          p[w].list[p[w].count] = v;
          p[w].count++;
        }
        else if (d[w] == (d[v] + 1))
        {
          sigma[w] += sigma[v];
          // Save the predecesor.
          p[w].list[p[w].count] = v;
          p[w].count++;
        }
      }
    }

    // While !empty(Stack)
    while (sPtr > 0)
    {
      sPtr--;
      int w = stack[sPtr];

      for (int k = 0; k < p[w].count; ++k)
      {
        // v = pred of w
        int v = p[w].list[k];
        delta[v] = delta[v] + (sigma[v] / sigma[w]) * (1.0 + delta[w]);
      }

      if (w != root)
      {
        bc[w] += delta[w];
      }
    }
  }

  // Free the memory.
  free(queue);
  free(stack);
  free(delta);
  free(d);
  free(sigma);
  free(pListMem);
  free(p);
}

