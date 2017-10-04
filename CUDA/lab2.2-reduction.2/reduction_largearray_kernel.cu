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
}

// **===-------- LAB 2.2 - Modify the body of this function -----------===**
// You may need to make multiple kernel calls, make your own kernel
// function in this file, and then call them from here.
void reductionArray(float *outArray, float *inArray, int numElements)
{
}
// **===-----------------------------------------------------------===**


#endif // _REDUCTION_KERNEL_CU_
