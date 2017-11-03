#ifndef _REDUCTION_KERNEL_CU_
#define _REDUCTION_KERNEL_CU_

// includes, kernels
#include <assert.h>
#include <cutil_inline.h>

// You can use any other block size you wish.
#define BLOCK_SIZE 256 

// **===-------- Lab 2.2 - Modify the body of this function -----------===**
// Reduction - Kernel 
__global__ void reduction_step(float *output, float *input, int nr_elm)
{
    __shared__  float scratch[BLOCK_SIZE*2];

    int tid = threadIdx.x;
    int stride;
    // Coalesced load to SHM
    int gidx = (blockIdx.x * BLOCK_SIZE * 2 + threadIdx.x);

    scratch[tid]   = input[gidx];
    scratch[tid+BLOCK_SIZE] = input[gidx+BLOCK_SIZE];

    __syncthreads();
    // build the sum tree (in-place)
    for (stride = BLOCK_SIZE; stride > 0; stride >>= 1) {
	if (threadIdx.x < stride) {
	    scratch[threadIdx.x] += scratch[threadIdx.x + stride];
	}
	__syncthreads();
    }
    // write results to global memory
    if(threadIdx.x == 0)
	output[0] = scratch[0];
}

// **===-------- LAB 2.2 - Modify the body of this function -----------===**
// You may need to make multiple kernel calls, make your own kernel
// function in this file, and then call them from here.
void reductionArray(float *outArray, float *inArray, int numElements)
{
    dim3 grid(1, 1, 1);
    dim3 block(BLOCK_SIZE, 1, 1);
    reduction_step<<<grid, block>>>(outArray, inArray, numElements); 
}
// **===-----------------------------------------------------------===**


#endif // _REDUCTION_KERNEL_CU_
