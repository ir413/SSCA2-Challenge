#include "sprng.h"

#include "BetweennessCentrality.h"

// Number of thread bloks.
#define THREAD_BLOCKS_COUNT 1
// Number of threads per block.
#define THREADS_PER_BLOCK 1024 

__global__ void computeBC(Graph *g, double *bc)
{

}

void computeBetweennessCentrality(Configuration *config, Graph *g, double *bc)
{
  // Zero bc.
  for (int i = 0; i < g->n; ++i)
  {
    bc[i] = 0.0;
  }

  // Stack of vertices in the order of non-decreasing
  // distance from s. Represents the BFS queue implicitly.
  //int *s;
  // Predecessoes of a vertex v on shortest paths from s.
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
  //s = (int *) malloc(g->n * sizeof(int));
  sigma = (double *) malloc(g->n * sizeof(double));
  d = (int *) malloc(g->n * sizeof(int));
  delta = (double *) calloc(g->n, sizeof(double));

  int *stack = (int *) malloc(g->n * sizeof(int));
  int *queue = (int *) malloc(g->n * sizeof(int));

  for (int root = 0; root < g->n; ++root)
  {
    int sPtr = 0;
    int qHead = 0;
    int qTail = 0;

    // Clear the auxilary structures.
    for (int k = 0; k < g->n; ++k)
    {
      d[k] = -1;
      sigma[k] = 0;
      p[k].count = 0;
    }

    sigma[root] = 1;
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
        else if (d[w] == (d[v] + 1))
        {
          sigma[w] += sigma[v];
          p[w].list[p[w].count] = v;
          p[w].count++;
        }
      }

    }

    // zero delta.
    for (int i = 0; i < g->n; ++i)
    {
      delta[i] = 0;
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

        delta[v] = delta[v] + sigma[v] * (1 + delta[w]) / sigma[w];
      }

      if (w != root)
      {
        bc[w] += delta[w];
      }
    }

  }


}
