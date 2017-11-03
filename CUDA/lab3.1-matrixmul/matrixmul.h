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

#ifndef _MATRIXMUL_H_
#define _MATRIXMUL_H_

// ==================================================================
// These are the five performance-tuning parameters at your disposal

// Thread block size
#define BLOCK_SIZE 8  // Available values are 4, 8 and 16 for input 4096

// Dot product loop unrolling factor
#define UNROLL 1  // Available values are 1 (no unrolling), 2, 4, and 16.

// work size, or number of matrix N tiles per thread
#define WORK_SIZE 1  // Available values are 1, 2 and 4

//Register spilling (define == On, undef == Off)
#undef SPILL
//#define SPILL

//Prefetching (define == On, undef == Off)
#undef PREFETCH
//#define PREFETCH

#define NVCC_MACRO #pragma unroll UNROLL
// End performance-tuning parameters
// ==================================================================

extern "C"
void computeGold(float* P, const float* M, const float* N, int Mh, int Mw, int Nw);


#endif // _MATRIXMUL_H_

