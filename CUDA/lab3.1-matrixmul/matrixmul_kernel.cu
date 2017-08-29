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
 * IMPLIED AwRRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL AwRRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED AwRRANTIES OF
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

/* Matrix multiplication: C = A * B.
 * Device code.
 */

#ifndef _MATRIXMUL_KERNEL_H_
#define _MATRIXMUL_KERNEL_H_

#include <stdio.h>
#include "matrixmul.h"

#define CHECK_BANK_CONFLICTS 0
#if CHECK_BANK_CONFLICTS
# define AS(i, j) CUT_BANK_CHECKER(((float *)&As[0][0]), (BLOCK_SIZE * i + j))
# define BS(i, j) CUT_BANK_CHECKER(((float *)&Bs[0][0]), (BLOCK_SIZE * i + j))
#else
# define AS(i, j) As[i][j]
# define BS(i, j) Bs[i][j]
#endif

////////////////////////////////////////////////////////////////////////////////
//! Simple test kernel for device functionality
//! @param g_idata  input data in global memory
//! @param g_odata  output data in global memory
////////////////////////////////////////////////////////////////////////////////

    __global__ void
matrixMul(float *C, float *A, float *B, int Aw, int Bw)
{
    // Block index
    #define bx blockIdx.x
    #define by blockIdx.y

    // Thread index
    #define tx threadIdx.x
    #define ty threadIdx.y
    // Index of the first element of A loaded by this thread by the block
    int aEnd = Aw * ty;
    aEnd += tx;
    int a = Aw * BLOCK_SIZE * by + aEnd;
    int b = BLOCK_SIZE * bx * WORK_SIZE;

    #ifdef SPILL
    // Create a shared-memory buffer to spill a register value
    // into shared memory, hopefully reducing the total required
    // register count.
    __shared__ int c[BLOCK_SIZE][BLOCK_SIZE];
    c[tx][ty] = a + b;
    #else
    int c =  a + b;
    #endif

    // Index of the first sub-matrix of B processed by the block
    b = b + aEnd;

    // Index of the last sub-matrix of A processed by the block
    aEnd = a + Aw;

    // Step size used to iterate through the sub-matrices of A
    #define aStep BLOCK_SIZE

    // Step size used to iterate through the sub-matrices of B
    #define bStep BLOCK_SIZE * Bw

    // Initialize result(s) to 0.
    float Csub = 0;
    #if (WORK_SIZE == 2 || WORK_SIZE == 4)
    float Dsub = 0;
    # if WORK_SIZE == 4
    float Esub = 0;
    float Fsub = 0;
    # endif // WORK_SIZE == 4
    #endif

    // Initial prefetch.  Issues loads to main memory and store
    // in temporary variables which will later be stored to shared memory
    #ifdef PREFETCH
    float fa = A[a];
    float fb = B[b];
    # if (WORK_SIZE == 2 || WORK_SIZE == 4)
    float fb2 = B[b + BLOCK_SIZE];
    #  if WORK_SIZE == 4
    float fb3 = B[b + BLOCK_SIZE * 2];
    float fb4 = B[b + BLOCK_SIZE * 3];
    #  endif // WORK_SIZE = 4
    # endif
    #endif

    // Loop over all the sub-matrices of A and B
    // required to compute the block sub-matrix
    while (a < aEnd) {
        int i;

        // Declaration of the shared memory array As used to
        // store the sub-matrix of A
        __shared__ float As[BLOCK_SIZE][BLOCK_SIZE];

        // Declaration of the shared memory array Bs used to
        // store the sub-matrices of B
        __shared__ float Bs[BLOCK_SIZE][BLOCK_SIZE * WORK_SIZE];

        #ifdef PREFETCH
        // If performing prefetching, the values are already loaded
        // from memory, and the temporary variables holding the loaded
        // values are stored to shared memory.
        As[ty][tx] = fa;
        Bs[ty][tx] = fb;
        # if (WORK_SIZE == 2 || WORK_SIZE == 4)
        Bs[ty][tx + BLOCK_SIZE] = fb2;
        #  if WORK_SIZE == 4
        Bs[ty][tx + BLOCK_SIZE * 2] = fb3;
        Bs[ty][tx + BLOCK_SIZE * 3] = fb4;
        #  endif // WORK_SIZE = 4
        # endif
        #else
        // Load the matrices from device memory
        // directly to shared memory
        As[ty][tx] = A[a];
        Bs[ty][tx] = B[b];
        # if (WORK_SIZE == 2 || WORK_SIZE == 4)
        Bs[ty][tx + BLOCK_SIZE] = B[b + BLOCK_SIZE];
        #  if WORK_SIZE == 4
        Bs[ty][tx + BLOCK_SIZE * 2] = B[b + BLOCK_SIZE * 2];
        Bs[ty][tx + BLOCK_SIZE * 3] = B[b + BLOCK_SIZE * 3];
        #  endif // WORK_SIZE = 4
        # endif
        #endif // PREFETCH
               // Update for next loop
        a += aStep;
        b += bStep;

        // Synchronize to make sure the shared memory
        // tiles are ready
        __syncthreads();

        #ifdef PREFETCH
        // If prefetching, issue the loads for the next tiles preemptively.
        // The loads will complete and be stored into these temporary
        // variables while the current shared memory tiles
        // are being operated on.
        fa = A[a];
        fb = B[b];
        # if (WORK_SIZE == 2 || WORK_SIZE == 4)
        fb2 = B[b + BLOCK_SIZE];
        #  if WORK_SIZE == 4
        fb3 = B[b + BLOCK_SIZE * 2];
        fb4 = B[b + BLOCK_SIZE * 3];
        #  endif
        # endif
        #endif

        // Multiply the two matrices together.
//#pragma unroll UNROLL
NVCC_MACRO
        for (i = 0; i < BLOCK_SIZE; i++) {
            Csub += As[ty][i] *Bs[i][tx];
            # if (WORK_SIZE == 2 || WORK_SIZE == 4)
            Dsub += As[ty][i] *Bs[i][tx + BLOCK_SIZE];
            #  if WORK_SIZE == 4
            Esub += As[ty][i] *Bs[i][tx + BLOCK_SIZE * 2];
            Fsub += As[ty][i] *Bs[i][tx + BLOCK_SIZE * 3];
            #  endif
            # endif
        }

        // Synchronize to make sure that the preceding
        // computation is done before overwriting new
        // shared memory sub-matrices of A and B in the next iteration
        __syncthreads();
    }

    #ifdef SPILL
    // If we spilled the output index at the beginning, load it back
    // from the shared memory array.
    // Output the result(s) for each thread.
    C[c[tx][ty]] = Csub;
    # if (WORK_SIZE == 2 || WORK_SIZE == 4)
    C[c[tx][ty] + BLOCK_SIZE] = Dsub;
    #  if WORK_SIZE == 4
    C[c[tx][ty] + BLOCK_SIZE * 2] = Esub;
    C[c[tx][ty] + BLOCK_SIZE * 3] = Fsub;
    #  endif
    # endif

    #else
    // Output the final result(s) for each thread.
    C[c] = Csub;
    # if (WORK_SIZE == 2 || WORK_SIZE == 4)
    C[c + BLOCK_SIZE] = Dsub;
    #  if WORK_SIZE == 4
    C[c + BLOCK_SIZE * 2] = Esub;
    C[c + BLOCK_SIZE * 3] = Fsub;
    #  endif
    # endif
    #endif
}

#endif // #ifndef _MATRIXMUL_KERNEL_H_

