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
    float *sigma,
    float *delta,
    plist *p)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  if (idx < n)
  {
    if (idx == source)
    {
      d[idx] = 0;
      sigma[idx] = 1.0;
    }
    else
    {
      d[idx] = -1; 
      sigma[idx] = 0.0;
      delta[idx] = 0.0;
    }
    p[idx].count = 0;
  }
}

__global__ void vertexParallelBFS(
    int source,
    Graph *g,
    int *d,
    float *sigma,
    plist *p,
    unsigned int *l)
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
          }
          
          if (d[w] == (d[v] + 1))
          {
            atomicAdd(&sigma[w], sigma[v]);
            // Save the predecessor.
            atomicExch(&(p[w].list[p[w].count]), v);
            atomicAdd(&(p[w].count), 1);
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

  __syncthreads();

  if (threadIdx.x == 0)
  {
    *l = level;
  }
}

__global__ void accumBC(
    int source,
    Graph *g,
    int *d,
    float *sigma,
    float *delta,
    plist *p,
    unsigned int l,
    float *bc)
{
  __shared__ int level;

  if (threadIdx.x == 0)
  {
    level = l;
  }

  __syncthreads();

  while (level > 0)
  {
    __syncthreads();

    for (int w = threadIdx.x; w < g->n; w += blockDim.x)
    {
      if (d[w] == level)
      {
        for (int k = 0; k < p[w].count; ++k)
        {
          int v = p[w].list[k];
       
          atomicAdd(&delta[v], (sigma[v] / sigma[w]) * (1.0 + delta[w]));
        }
      } 

      __syncthreads();
      if (w != source)
      {
        atomicAdd(&bc[w], delta[w]);
      }
    }

    /*
    __syncthreads();
    for (int v = threadIdx.x; v < g->n; v += blockDim.x)
    {
      if (v != source)
      {
        atomicAdd(&bc[v], delta[v]);
      }
    }
    */

    if (threadIdx.x == 0)
    {
      level--;
    }
    __syncthreads();
  }

  /*
  __syncthreads();

  for (int v = threadIdx.x; v < g->n; v += blockDim.x)
  {

    atomicAdd(&bc[v], delta[v]);
    
    //if (v != source)
    //{
    //  atomicAdd(&bc[v], delta[v]);
    //}
  }
  */
}

void computeBCGPU(Configuration *config, Graph *g, int *perm, float *bc)
{
  // Declare the auxilary structures.
  int *d;
  float *sigma;
  float *delta;
  int *stack;
  plist *p;
  int *pListMem;
  unsigned int *l;

  // Allocate temporary structures in global memory.
  cudaMallocManaged(&d, g->n * sizeof(int));
  cudaMallocManaged(&sigma, g->n * sizeof(float));
  cudaMallocManaged(&delta, g->n * sizeof(float));
  cudaMallocManaged(&stack, g->n * sizeof(int));
  cudaMallocManaged(&p, g->n * sizeof(plist));
  cudaMallocManaged(&pListMem, g->m * sizeof(int));
  cudaMallocManaged(&l, sizeof(unsigned int));

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
    p[i].count = 0;
  }

  cudaFree(inDegree);
  cudaFree(numEdges);
  // --- TMP

  for (int source = 0; source < g->n; ++source)
  {
    *l = 0;

    //printf("before %d \n", *l);

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
        l);
    cudaDeviceSynchronize();

    //printf("after %d \n", *l);
  
    // Sum centrality scores.
    accumBC<<<BLOCKS_COUNT, MAX_THREADS_PER_BLOCK>>>(
        source,
        g,
        d,
        sigma,
        delta,
        p,
        *l,
        bc);
    cudaDeviceSynchronize();

    printf("D:\n");
    for (int i = 0; i < g->n; ++i)
    {
      printf("%d, ", d[i]);
    }
    printf("\nSigma:\n");
    for (int i = 0; i < g->n; ++i)
    {
      printf("%f, ", sigma[i]);
    }
    printf("\n");

    /*
    printf("P:\n");
    for (int i = 0; i < g->n; ++i)
    {
      printf("p[%d].count = %d\n  ", i, p[i].count);
      for (int j = 0; j < p[i].count; ++j)
      {
        printf("%d, ", p[i].list[j]);
      }
      printf("\n");
    }
    */
   
    printf("Delta:\n");
    for (int i = 0; i < g->n; ++i)
    {
      printf("%f, ", delta[i]);
    }
    printf("\n");
  }

  // Clean up.
  cudaFree(l);
  cudaFree(pListMem);
  cudaFree(p);
  cudaFree(stack);
  cudaFree(delta);
  cudaFree(sigma);
  cudaFree(d);
}

void computeBCCPU(Configuration *config, Graph *g, int *perm, float *bc)
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
  float *sigma;
  // Shortest path length between the pairs.
  int *d;
  // Dependency of vertices.
  float *delta;
 
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
    p[i].count = 0;
  }

  // Clean up temporary structures.
  free(inDegree);
  free(numEdges);

  // Allocate memory.
  sigma = (float *) malloc(g->n * sizeof(float));
  d = (int *) malloc(g->n * sizeof(int));
  delta = (float *) calloc(g->n, sizeof(float));

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

        // Not visited.
        if (d[w] == -1)
        {
          // Enqueue w -> Q
          queue[qTail] = w;
          qTail++;
        
          // Distance to w is distance to v + 1.
          d[w] = d[v] + 1;
        }

        if (d[w] == (d[v] + 1))
        {
          sigma[w] += sigma[v];
          // Save the predecesor.
          p[w].list[p[w].count] = v;
          p[w].count++;
        }
      }
    }

    printf("Root: %d\n", root);
    printf("D:\n");
    for (int i = 0; i < g->n; ++i)
    {
      printf("%d, ", d[i]);
    }
    printf("\nSigma:\n");
    for (int i = 0; i < g->n; ++i)
    {
      printf("%f, ", sigma[i]);
    }
    printf("\n");
    /*
    printf("P:\n");
    for (int i = 0; i < g->n; ++i)
    {
      printf("p[%d].count = %d\n  ", i, p[i].count);
      for (int j = 0; j < p[i].count; ++j)
      {
        printf("%d, ", p[i].list[j]);
      }
      printf("\n");
    }
    */

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

    printf("Delta:\n");
    for (int i = 0; i < g->n; ++i)
    {
      printf("%f, ", delta[i]);
    }
    printf("\n");
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

