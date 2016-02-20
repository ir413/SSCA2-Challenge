#include <assert.h>
#include <stdio.h>

#include "PrefixSums.h"


int main(int argc, char **argv)
{
  int n = 2048;

  int *data;
  int *sums;

  // Allocate memory.
  cudaMallocManaged(&data, sizeof(int) * n);
  cudaMallocManaged(&sums, sizeof(int) * n);

  // Initialize data.
  for (int i = 0; i < n; ++i)
  {
    data[i] = 1; 
    sums[i] = 0;
  }

  prefixSums(data, sums, n);

  // Check the result.
  for (int i = 0; i < n; ++i)
  {
    assert(sums[i] == i);
  }

  // Clean up.
  cudaFree(sums);
  cudaFree(data);
}


