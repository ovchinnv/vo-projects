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
extern "C" void matmul_cuda(float *M, float *N, float *P, int my, int mx, int nx) {
    StopWatchInterface *timer_compute = NULL;
    StopWatchInterface *timer_memory = NULL;
    unsigned int mem_size_M, mem_size_N, mem_size_P;
    float * deviceM = NULL, * deviceN = NULL, * deviceP = NULL;
    int block_size = 8;

    mem_size_M=my*mx*sizeof(float);
    mem_size_N=mx*nx*sizeof(float);
    mem_size_P=my*nx*sizeof(float);

    sdkCreateTimer(&timer_memory);
    sdkStartTimer(&timer_memory);

    printf("  Allocate device memory.\n");
    checkCudaErrors(cudaMalloc(&deviceM, mem_size_M ));
    checkCudaErrors(cudaMalloc(&deviceN, mem_size_N ));

    printf("  Copy host memory to device.\n");
    checkCudaErrors(cudaMemcpy(deviceM, M, mem_size_M, cudaMemcpyHostToDevice ));
    checkCudaErrors(cudaMemcpy(deviceN, N, mem_size_N, cudaMemcpyHostToDevice ));

    printf("  Allocate device memory for results.\n");
    checkCudaErrors(cudaMalloc(&deviceP, mem_size_P ));

    // Clear device memory
    cudaMemset(deviceP, 0, mem_size_P);

    sdkStopTimer(&timer_memory);

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

    sdkStartTimer(&timer_memory);

    printf("  Copy result from device to host.\n");
    cudaMemcpy(P, deviceP, mem_size_P, cudaMemcpyDeviceToHost );

    sdkStopTimer(&timer_memory);

    // ================================================
    // Show timing information
    // ================================================

    printf("  GPU memory access time: %f (ms)\n",
        sdkGetTimerValue(&timer_memory));
    printf("  GPU computation time  : %f (ms)\n",
        sdkGetTimerValue(&timer_compute));
    printf("  GPU processing time   : %f (ms)\n",
        sdkGetTimerValue(&timer_compute) + sdkGetTimerValue(&timer_memory));
    sdkDeleteTimer(&timer_memory);
    sdkDeleteTimer(&timer_compute);

    checkCudaErrors(cudaFree(deviceM ));
    checkCudaErrors(cudaFree(deviceN ));
    checkCudaErrors(cudaFree(deviceP ));

}
