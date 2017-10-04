/*
 * Copyright 1993-2006 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOURCE CODE.
 *
 * U.S. Government End Users.  This source code is a "commercial item" as
 * that term is defined at 48 C.F.R. 2.101 (OCT 1995), consisting  of
 * "commercial computer software" and "commercial computer software
 * documentation" as such terms are used in 48 C.F.R. 12.212 (SEPT 1995)
 * and is provided to the U.S. Government only as a commercial end item.
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
 * source code with only those rights set forth herein.
 */

/* Matrix multiplication: P = M * N.
 * Host code.
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <cutil.h>

// includes, kernels
#include "matrixmul_kernel.cu"

#include "assist.h"

#define ERROR_CHECK { cudaError_t err; \
  if ((err = cudaGetLastError()) != cudaSuccess) { \
    printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__);}}

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
void runTest(int argc, char** argv);

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main(int argc, char** argv)
{
    runTest(argc, argv);
}

////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void
runTest(int argc, char** argv)
{
    bool if_quiet = true;
    unsigned int timer_compute = 0;
    unsigned int timer_memory = 0;
    int i, j;
    char *matrix_id = NULL, *input_fn = NULL, *gold_fn = NULL;
    float * deviceM = NULL, * deviceN = NULL, * deviceP = NULL;
    int Mw = 0, Mh = 0, Nw = 0, Nh = 0, Pw = 0, Ph = 0;
    int block_size = 0;

    if (argc == 2) {
        matrix_id = strdup(argv[1]);

    } else {
        fprintf(stderr, "Error: Wrong input parameter numbers.\n");
        fprintf(stderr, "Usage:\n"
                        "$> ./lab3.1-matrixmul <8, 128, 512, 3072, 4096>\n"
                        "Examples:\n"
                        "      $> ./lab3.1-matrixmul 128\n"
                        );
        exit(1);
    }

    // Note: Matrix width and height must be multiples of block size.
    if (!strcmp(matrix_id, "8")) {
        Mw = Mh = Nw = Nh = Pw = Ph = 8;
        block_size = BLOCK_SIZE; // thread number = block_size^2
        input_fn = strdup("matrix_8.bin");
        gold_fn = strdup("matrix_8.gold");
        if_quiet = false; // If not display matrix contents
    } else
    if (!strcmp(matrix_id, "128")) {
        Mw = Mh = Nw = Nh = Pw = Ph = 128;
        block_size = BLOCK_SIZE; // thread number = block_size^2
        input_fn = strdup("matrix_128.bin");
        gold_fn = strdup("matrix_128.gold");
        if_quiet = true; // If not display matrix contents
    } else
    if (!strcmp(matrix_id, "512")) {
        Mw = Mh = Nw = Nh = Pw = Ph = 512;
        block_size = BLOCK_SIZE; // thread number = block_size^2
        input_fn = strdup("matrix_512.bin");
        gold_fn = strdup("matrix_512.gold");
        if_quiet = true; // If not display matrix contents
    } else
    if (!strcmp(matrix_id, "3072")) {
        Mw = Mh = Nw = Nh = Pw = Ph = 3072;
        block_size = BLOCK_SIZE; // thread number = block_size^2
        input_fn = strdup("matrix_3072.bin");
        gold_fn = strdup("matrix_3072.gold");
        if_quiet = true; // If not display matrix contents
    } else
    if (!strcmp(matrix_id, "4096")) {
        Mw = Mh = Nw = Nh = Pw = Ph = 4096;
        block_size = BLOCK_SIZE; // thread number = block_size^2
        input_fn = strdup("matrix_4096.bin");
        gold_fn = strdup("matrix_4096.gold");
        if_quiet = true; // If not display matrix contents
    } else {
        printf("***Error on %s: %d: Undefined matrix ID.\n",
            __FILE__, __LINE__);
        printf("   You should add it to the source code.\n");
        printf("   Current available ID's are 8, 128, 512, 3072, 4096.\n");
        exit(1);
    }

    if (block_size > Mw) {
        printf("***Error on %s: %d: Block size %d is larger than matrix width %d.\n",
            __FILE__, __LINE__, block_size, Mw);
        printf("   You should define a smaller block size.\n");
        exit(1);
    }

    printf("Input matrix file name: %s\n", input_fn);

    // -----------------------------------------------------------------------
    // Setup host side
    // -----------------------------------------------------------------------

    printf("Setup host side environment and launch kernel:\n");

    // allocate host memory for matrices M and N
    printf("  Allocate host memory for matrices M and N.\n");
    printf("    M: %d x %d\n", Mw, Mh);
    printf("    N: %d x %d\n", Nw, Nh);
    unsigned int size_M = Mw * Mh;
    unsigned int mem_size_M = sizeof(float) * size_M;
    float* hostM = (float*) malloc(mem_size_M);
    unsigned int size_N = Nw * (Nh);
    unsigned int mem_size_N = sizeof(float) * size_N;
    float* hostN = (float*) malloc(mem_size_N);

    // allocate memory for the result on host side
    printf("  Allocate memory for the result on host side.\n");
    unsigned int size_P = Pw * Ph;
    unsigned int mem_size_P = sizeof(float) * size_P;
    float* hostP = (float*) malloc(mem_size_P);

    // Initialize the input matrices.
    printf("  Initialize the input matrices.\n");
    unsigned int * matrix = ReadMatrixFile(input_fn, Pw, Ph, if_quiet);
    for (i = 0; i < Mw; i++)
        for (j = 0; j < Nw; j++)
	        hostM[i * Mw + j] = hostN[i * Mw + j] = (float) matrix[i*Mw + j];
    free(matrix); matrix = NULL;

    // ===================================================================
    //  Allocate device memory for the input matrices.
    //  Copy memory from the host memory to the device memory.
    // ===================================================================

    CUT_SAFE_CALL(cutCreateTimer(&timer_memory));
    CUT_SAFE_CALL(cutStartTimer(timer_memory));

    printf("  Allocate device memory.\n");

    CUDA_SAFE_CALL(cudaMalloc((void**) &deviceM, mem_size_M));
    CUDA_SAFE_CALL(cudaMalloc((void**) &deviceN, mem_size_N));

    printf("  Copy host memory data to device.\n");

    CUDA_SAFE_CALL(cudaMemcpy(deviceM, hostM, mem_size_M,
        cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(deviceN, hostN, mem_size_N,
        cudaMemcpyHostToDevice));

    printf("  Allocate device memory for results.\n");

    CUDA_SAFE_CALL(cudaMalloc((void**) &deviceP, mem_size_P));
    cudaMemset(deviceP, 0, mem_size_P);

    CUT_SAFE_CALL(cutStopTimer(timer_memory));

    // ===================================================================
    // Initialize the thread block and kernel grid dimensions
    // and invoke the CUDA kernel.
    // You may assume that each matrix dimension is a multiple
    // of the defined constant block_size.
    // ===================================================================

    printf("  Setup kernel execution parameters.\n");

    // Different ways of declarations
    #if 1
    dim3 block;
    dim3 grid;
    grid.x = Pw/block_size/WORK_SIZE;  // this must imply that each thread computes WORK_SIZE elements of P
    grid.y = Pw/block_size;
    block.x = block_size;
    block.y = block_size;
    #else
    dim3 block(block_size, block_size);
    dim3 grid(Pw/block.x/WORK_SIZE, Pw/block.y);
    #endif

    printf("  # of threads in a block: %d x %d (%d)\n",
        block.x, block.y, block.x * block.y);
    printf("  # of blocks in a grid  : %d x %d (%d)\n",
        grid.x, grid.y, grid.x * grid.y);

    // ================================================
    // Initialize the block and grid dimensions here
    // ================================================

    printf("  Executing the kernel...\n");

    printf("  Optimization parameters:\n");
    printf("    Block size: %d\n", BLOCK_SIZE);
    printf("    Unrolling factor: %d\n", UNROLL);
    printf("    Working size: %d\n", WORK_SIZE);
    #ifdef SPILL
    printf("    Register spilling: %s\n", "enabled");
    #else
    printf("    Register spilling: %s\n", "disabled");
    #endif
    #ifdef PREFETCH
    printf("    Data prefetch: %s\n", "enabled");
    #else
    printf("    Data prefetch: %s\n", "disabled");
    #endif

    // Start the timer_compute to calculate how much time we spent on it.
    CUT_SAFE_CALL(cutCreateTimer(&timer_compute));
    CUT_SAFE_CALL(cutStartTimer(timer_compute));

    // Invoke the CUDA kernel here
    matrixMul<<<grid, block>>> (deviceP, deviceM, deviceN, Mw, Nw);

    // Make sure all threads have finished their jobs
    // before we stop the timer_compute.
    cudaThreadSynchronize();

    // Stop the timer_compute
    CUT_SAFE_CALL(cutStopTimer(timer_compute));

    // check if kernel execution generated an error
    ERROR_CHECK
    CUT_CHECK_ERROR("Kernel execution failed");

    // ===================================================================
    // Copy the results back from the host
    // ===================================================================

    printf("  Copy result from device to host.\n");

    CUT_SAFE_CALL(cutStartTimer(timer_memory));
    cudaMemcpy(hostP, deviceP, mem_size_P, cudaMemcpyDeviceToHost);
    CUT_SAFE_CALL(cutStopTimer(timer_memory));

    // ================================================
    // Show timing information
    // ================================================

    printf("  GPU memory access time: %f (ms)\n",
        cutGetTimerValue(timer_memory));
    printf("  GPU computation time  : %f (ms)\n",
        cutGetTimerValue(timer_compute));
    printf("  GPU processing time   : %f (ms)\n",
        cutGetTimerValue(timer_compute) + cutGetTimerValue(timer_memory));
    CUT_SAFE_CALL(cutDeleteTimer(timer_memory));
    CUT_SAFE_CALL(cutDeleteTimer(timer_compute));

    // ================================================
    // Do comparison
    // ================================================

    // Full result check when input matrix is <= 512x512
    //if (0) {
    if (Mw * Nw > 512*512) {
        printf("\nInput matrix size is too big. Skip computing reference.\n");

    } else {
        printf("\nCheck results with those computed by CPU.\n");
        printf ("  Computing reference solution.\n");
        CUT_SAFE_CALL(cutCreateTimer(&timer_compute));
        CUT_SAFE_CALL(cutStartTimer(timer_compute));

        float* reference = (float*) malloc(mem_size_P);
        computeGold(reference, hostM, hostN, Mh, Mw, Nw);
        CUT_SAFE_CALL(cutStopTimer(timer_compute));

        printf("  CPU Processing time   : %f (ms)\n",
            cutGetTimerValue(timer_compute));
        CUT_SAFE_CALL(cutDeleteTimer(timer_compute));

        printf("  CPU checksum: %g\n", CheckSum(reference, Mw, Nw));

        matrix = (unsigned int *) malloc(Pw * Ph * sizeof(unsigned int));
        for (i = 0; i < Ph; i++)
            for (j = 0; j < Pw; j++)
                matrix[i*Pw + j] = (unsigned int) reference[i*Pw + j];

        WriteMatrixFile("lab3.1-matrixmul.gold", matrix, Pw, Ph, 1);
        free(matrix); matrix = NULL;
        free(reference);
    }

    printf("  GPU checksum: %g\n", CheckSum(hostP, Mw, Nw));

    /* Write matrix C to output binary file */
    matrix = (unsigned int *) malloc (Pw * Ph * sizeof(unsigned int));
    for (i = 0; i < Ph; i++)
        for (j = 0; j < Pw; j++)
            matrix[i*Pw + j] = (unsigned int) hostP[i*Pw + j];
    WriteMatrixFile("lab3.1-matrixmul.bin", matrix, Pw, Ph, 1);
    free (matrix); matrix = NULL;

    if (Mw >= 3072 && Mh >= 3072) {
        CompareMatrixFile("lab3.1-matrixmul.bin", gold_fn, Pw, Ph, if_quiet);
    } else {
        CompareMatrixFile("lab3.1-matrixmul.bin", "lab3.1-matrixmul.gold",
            Pw, Ph, if_quiet);
    }

    // clean up memory
    free(hostM); free(hostN); free(hostP);
    free(input_fn); free(gold_fn);

    // ===================================================================
    // Free the device memory
    // ===================================================================

    CUDA_SAFE_CALL(cudaFree(deviceM));
    CUDA_SAFE_CALL(cudaFree(deviceN));
    CUDA_SAFE_CALL(cudaFree(deviceP));
}


