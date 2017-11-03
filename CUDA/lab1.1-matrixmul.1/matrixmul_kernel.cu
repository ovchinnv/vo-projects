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
 * Device code.
 */

#ifndef _MATRIXMUL_KERNEL_H_
#define _MATRIXMUL_KERNEL_H_

//#include "matrixmul.h"

////////////////////////////////////////////////////////////////////////////////
//! Simple test kernel for device functionality
//! @param g_idata  input data in global memory
//! @param g_odata  output data in global memory
////////////////////////////////////////////////////////////////////////////////

    __global__ void matrixMul( // device kernel
    float* P, const float* M, const float* N,
    const int Mh, const int Mw, const int Nw) {

    float Psub = 0;
    int i, j;
    unsigned int im, in, ip;
    // ===================================================================
    // Begin solution part 5
    // Determine the output index of each thread.
    // Compute the dot product of one row of M and one column of N
    // for each thread.
    // Write the computed value to matrix P at the correct index.
    // ===================================================================

    //
    // find index of next row of M (remember that columns vary fast ; i.e. M[ irow * ncol + icol ] = M(irow, icol)
    // recall that each thread corresponds to a computed entry
    i=blockIdx.x * blockDim.x + threadIdx.x ; // take this as the row index of M
    j=blockIdx.y * blockDim.y + threadIdx.y ; // take this as the column index of N
    // compute P(i,j) = M(i,:) * N(:,j)
    // For each index from [0, Width of M)
    for (int k = 0; k < Mw; k++) {
        // Multiply the corresponding elements of M and N, and accumulate
        // into partial sum Psub.
        im=i*Mw + k;
        in=k*Nw + j;
        Psub += M[im] * N[in]; // note integer index multiplications are cheap
    }
    // P has size( Mh x Nw )
    ip=i*Nw + j;
    P[ip] = Psub;
    // End of solution part 5 ============================================
}

#endif // #ifndef _MATRIXMUL_KERNEL_H_
