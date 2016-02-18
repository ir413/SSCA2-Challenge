#include <assert.h>
#include <stdio.h>

/**
 * Adds two vectors.
 */
__global__ void
vectorAdd(const float *a, const float *b, float *c, size_t size)
{
  size_t i = threadIdx.x;

  if (i < size)
  {
    c[i] = a[i] + b[i];
  }
}

/**
 * Entry point.
 */
int
main(int argc, char **argv)
{
  size_t size = 100;
  
  float *a;
  float *b;
  float *c;

  // Allocate memory.
  cudaMallocManaged(&a, sizeof(float) * size);
  cudaMallocManaged(&b, sizeof(float) * size);
  cudaMallocManaged(&c, sizeof(float) * size);

  // Initialize a & b.
  for (size_t i = 0; i < size; ++i)
  {
    a[i] = 1.0f;
    b[i] = 1.0f; 
  }

  vectorAdd<<<1, 100>>>(a, b, c, size);
  cudaDeviceSynchronize();

  // Check the result.
  for (size_t i = 0; i< size; ++i)
  {
    assert(abs(c[i] - 2.0f) < 0.0001f);
  }

  // Clean up.
  cudaFree(c);
  cudaFree(b);
  cudaFree(a);
}


