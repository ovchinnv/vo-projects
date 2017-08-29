/* Matrix multiplication: P = M * N.
 * Device code.
 */

#ifndef _MATRIXMUL_KERNEL_H_
#define _MATRIXMUL_KERNEL_H_

#include <stdio.h>
//#include "matrixmul.h"


////////////////////////////////////////////////////////////////////////////////
//! Matrix multiplication on the device: P = M * N
//! Mw is M's width and Nw is N's width
////////////////////////////////////////////////////////////////////////////////
    __global__ void matrixMul( float* P, float* M, float* N, int Mh, int Mw, int Nw) {
#define bx blockIdx.x
#define by blockIdx.y
#define tx threadIdx.x
#define ty threadIdx.y

    const int irow=bx*BLOCK_SIZE_X+tx; // row of P
    const int icol=by*BLOCK_SIZE_Y+ty; // col of P

    __shared__ float Ms[BLOCK_SIZE_X * BLOCK_SIZE_Y];
    __shared__ float Ns[BLOCK_SIZE_X * BLOCK_SIZE_Y];

    int j, k;
    int i1, i2;
    float Pval=0;

    // ===================================================================
    // load submatrices
#define IO(i,j) (Mw * ( i ) + ( j ))
#define II(i,j) (BLOCK_SIZE_Y * ( i ) + ( j ))
    for ( j=0 ; j < Mw ; j+=BLOCK_SIZE_Y ) {
      Ms[II(tx,ty)] = M[IO(irow, j + ty )];
      Ns[II(tx,ty)] = N[IO(j+tx, icol)];
      __syncthreads();

// note : no penalty benefit to unrolling, or integer index computations
#define _UNROLL 8
      for (k = 0; k < BLOCK_SIZE_Y; k+=_UNROLL) {
         i1=II(tx,k);
         i2=II(k,ty);
         Pval+=Ms[i1]*Ns[i2];
//         Pval+=Ms[II(tx,k)]*Ns[II(k,ty)];
#if _UNROLL >= 2
         i1++;
         i2+=BLOCK_SIZE_Y;
         Pval+=Ms[i1]*Ns[i2];
#endif
#if _UNROLL >= 4
         i1++;
         i2+=BLOCK_SIZE_Y;
         Pval+=Ms[i1]*Ns[i2];
         i1++;
         i2+=BLOCK_SIZE_Y;
         Pval+=Ms[i1]*Ns[i2];
#endif
#if _UNROLL >= 8
         i1++;
         i2+=BLOCK_SIZE_Y;
         Pval+=Ms[i1]*Ns[i2];
         i1++;
         i2+=BLOCK_SIZE_Y;
         Pval+=Ms[i1]*Ns[i2];
         i1++;
         i2+=BLOCK_SIZE_Y;
         Pval+=Ms[i1]*Ns[i2];
         i1++;
         i2+=BLOCK_SIZE_Y;
         Pval+=Ms[i1]*Ns[i2];
#endif
        }

        __syncthreads(); // unclear to me why needed, but is needed

    }

    // ===================================================================
    // Code segment 3
    // Store the data back to global memory
    // ===================================================================
 P[IO(irow,icol)]=Pval;

#undef IO
#undef II
}

#endif // #ifndef _MATRIXMUL_KERNEL_H_
