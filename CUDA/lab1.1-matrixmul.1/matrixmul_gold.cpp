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

#include <stdio.h>
#include "matrixmul.h"

////////////////////////////////////////////////////////////////////////////////
//! Compute reference data set
//! P = M * N
//! @param P          reference data, computed but preallocated
//! @param M          matrix M as provided to device
//! @param N          matrix N as provided to device
//! @param Mh         height of matrix M
//! @param Nw         width of matrix N
////////////////////////////////////////////////////////////////////////////////
void
computeGold(float* P, const float* M, const float* N, int Mh, int Mw, int Nw)
{
  int i, j, k;
  float sum, a, b;

  for (i = 0; i < Mh; i++) // rows
    for (j = 0; j < Nw; j++) // columns
      {
	    sum = 0;
	    for (k = 0; k < Mw; k++)
	    {
	        a = M[i * Mw + k]; // skip i rows, then take k column : M(i,k) ; note that column varies fast
	        b = N[k * Nw + j]; // skip k rowsm=, take jth col :     N(k,j)
            //printf ("A[%d] * B[%d]\n", i * Mw + k, k * Nw + j);
	        sum += a * b;
	    }
	    P[i * Nw + j] = (float)sum;
      }
}
