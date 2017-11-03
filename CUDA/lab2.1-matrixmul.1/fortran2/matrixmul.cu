#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <helper_cuda.h>
#include <helper_timer.h>

// includes, kernels
#include "matrixmul_kernel.cu"

#define ERROR_CHECK { cudaError_t err; \
  if ((err = cudaGetLastError()) != cudaSuccess) { \
    printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__);}}

////////////////////////////////////////////////////////////////////////////////
extern "C" void matmul_cuda(float *deviceM, float *deviceN, float *deviceP, int my, int mx, int nx) {
    StopWatchInterface *timer_compute = NULL;
    int block_size = 8;

    printf("  Setup kernel execution parameters.\n");
    dim3 block(block_size, block_size); // thread indexing : there will be block_size ^ 2 threads per block, indexed (threadIdx.x threadIdx.y)
    dim3 grid(nx/block.x, my/block.y);

    printf("  # of threads in a block: %d x %d (%d)\n",
        block.x, block.y, block.x * block.y);
    printf("  # of blocks in a grid  : %d x %d (%d)\n",
        grid.x, grid.y, grid.x * grid.y);

    // ================================================
    // Initialize the block and grid dimensions here
    // ================================================

    printf("  Executing the kernel...\n");
    // Start the timer_compute to calculate how much time we spent on it.
    sdkCreateTimer(&timer_compute);
    sdkStartTimer(&timer_compute);

    // Invoke the CUDA kernel here
    matrixMul<<<grid, block>>>(deviceP, deviceM, deviceN, my, mx, nx);

    // Make sure all threads have finished their jobs
    // before we stop the timer_compute.
    cudaThreadSynchronize();

    // Stop the timer_compute
    sdkStopTimer(&timer_compute);

    // check if kernel execution generated an error
    ERROR_CHECK
    getLastCudaError("Kernel execution failed");

    // Copy the results back from the host
    // ===================================================================

    // ================================================
    // Show timing information
    // ================================================

    printf("  GPU computation time  : %f (ms)\n",
        sdkGetTimerValue(&timer_compute));
    sdkDeleteTimer(&timer_compute);

}
