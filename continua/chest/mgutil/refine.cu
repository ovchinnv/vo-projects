 __global__ void Refine_Cuda_3D(__CUFLOAT *fine,
#ifndef __MGTEX
  __CUFLOAT *coarse,
#endif
 const __CINT i3f, const __CINT i3c, const __CINT nnx, const __CINT nny, const __CINT nnz) {

#ifdef __MGTEX
//use texture addressing
#define coarse(_IND) tex1Dfetch(texcoarse,_IND)
#else
#define coarse(_IND) coarse[_IND]
#endif

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z


}