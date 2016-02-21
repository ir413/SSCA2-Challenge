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
}
