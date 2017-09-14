#ifdef __TEX
 texture<__CTYPE> tex_reduce; // may not work for types other than float
#endif

#define _BDIM 128 // threadblock size
#define _VPT 32    // values per thread

#include <stdio.h>
#include "kernel.cu"

#define _NBLK(N,n) (N)/(n) + ( (N) % (n) > 0 )


extern "C" void reduction_c(__CTYPE *A, int n, __CTYPE *val) {


  __CTYPE *devA, *devAout ;
  unsigned int numblk = _NBLK ( n , 2*_BDIM*_VPT ) ;
  __CTYPE redA[numblk];

// if (n<10000) 
//for (int i=0;i<n;i++) { printf("%12.5f\n",A[i]);}

 printf("N :%5d\n",n);

  cudaMalloc( &devA, n*sizeof(__CTYPE) ); // device copy of A
  cudaMemcpy( devA, A, n*sizeof(__CTYPE), cudaMemcpyHostToDevice);

  cudaMalloc( &devAout, numblk*sizeof(__CTYPE) ); // device copy of reduced array data

  dim3 block ( _BDIM, 1, 1 ); // thread indices
  dim3 grid  ( numblk, 1, 1 ); // block indices
 printf("threads/block :%5d\n",_BDIM);
 printf("blocks :%5d\n",numblk);
 printf("values/thread :%5d\n",_VPT);

// launch kernel

#ifdef __TEX
  cudaBindTexture(NULL, tex_reduce, devA, n*sizeof(__CTYPE));
#endif
 for (int i=0 ; i<1000;i++){
  reduction_cuda<<<grid,block>>>(devA, devAout, n);
//  reduction_cuda_simple<<<grid,block>>>(devA, devAout, n);
// reduction_cuda_simple2<<<grid,block>>>(devA, devAout, n);
// reduction_cuda_simple3<<<grid,block>>>(devA, devAout, n);
 }
#ifdef __TEX
  cudaUnbindTexture(tex_reduce);
#endif

// copy partially reduced array
  cudaMemcpy( redA, devAout, numblk*sizeof(__CTYPE), cudaMemcpyDeviceToHost);
//
// compute final value
//
#define _FMT "%12.5f\n"
#define _FMT "%5d\n"
  *val=0.f ; for (int i=0;i<numblk;i++) { *val+=redA[i] ;}//  printf(_FMT,redA[i]);}
  printf(_FMT,*val);
// free device memory
  cudaFree(devA);
  cudaFree(devAout);

}