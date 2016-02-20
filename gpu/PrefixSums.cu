#include "PrefixSums.h"

// Number of thread bloks.
#define THREAD_BLOCKS_COUNT 1
// Number of threads per block.
#define THREADS_PER_BLOCK 1024 

/**
 * Basic work-efficient parallel prefix sum implementation,
 * as described in "Parallel Prefix Sum (Scan) with CUDA" by Mark Harris.
 *
 * Pre: n is a power of 2
 */
__global__ void prefixSumsKernel(int *data, int *sums, int n)
{
  // Shared between the block.
  __shared__ int worklist[2 * THREADS_PER_BLOCK];

  int threadId = threadIdx.x;
  int offset = 1;

  // Load the data.
  worklist[2 * threadId] = data[2 * threadId];
  worklist[2 * threadId + 1] = data[2 * threadId + 1];

  // Perform the up-sweep phase.
  for (int d = n >> 1; d > 0; d >>= 1)
  {
    __syncthreads();

    if (threadId < d)
    {
      int ai = offset * (2 * threadId + 1) - 1;
      int bi = offset * (2 * threadId + 2) - 1;

      worklist[bi] += worklist[ai];
    }

    offset <<= 1;
  }

  // Clear the root.
  if (threadId == 0)
  {
    worklist[n - 1] = 0;
  } 

  // Perform the down-sweep phase.
  for (int d = 1; d < n; d <<= 1)
  {
    offset >>= 1;
    __syncthreads();

    if (threadId < d)
    {
      int ai = offset * (2 * threadId + 1) - 1;
      int bi = offset * (2 * threadId + 2) - 1;

      int t = worklist[ai];
      worklist[ai] = worklist[bi];
      worklist[bi] += t;
    }
  }

  __syncthreads();

  // Write the results.
  sums[2 * threadId] = worklist[2 * threadId];
  sums[2 * threadId + 1] = worklist[2 * threadId + 1];
}

void prefixSums(int *data, int *sums, int n)
{
  prefixSumsKernel<<<THREAD_BLOCKS_COUNT, THREADS_PER_BLOCK>>>(data, sums, n);
  cudaDeviceSynchronize();
}

