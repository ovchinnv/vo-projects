 __global__ void Coarsen_Cuda_3D(
#ifndef __MGTEX
  __CUFLOAT *fine,
#endif
  __CUFLOAT *coarse, const __CINT i3, const __CINT nnx, const __CINT nny, const __CINT nnz, const __CINT ibc) {

#ifdef __MGTEX
//use texture addressing
#define fine(_IND) tex1Dfetch(texfine,_IND)
#else
#define fine(_IND) fine[_IND]
#endif

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define ixc (blockIdx.x * blockDim.x + tx )
#define iyc (blockIdx.y * blockDim.y + ty )
#define izc (blockIdx.z * blockDim.z + tz )

#define ixf (ixc<<1)
#define iyf (iyc<<1)
#define izf (izc<<1)

// coarse index accounts for possible BC ; nnx are coarse lengths
#define IDC(i,j,k) ( i3 - 1 + (k+ibc)*(nnx+2*ibc)*(nny+2*ibc) + (j+ibc)*(nnx+2*ibc)   + (i+ibc) )
// fine index accounts for possible BC points ; index is assumed to start at 0 ; fine lengths are twice the coarse lengths ; hence taking 2 out
#define IDF(i,j,k) (          (k+ibc)*(4*(nnx+ibc)*(nny+ibc)) + (j+ibc)*2*(nnx+ibc)   + (i+ibc) )

#define coef 0.1250000000000

 coarse[ IDC(ixc,iyc,izc) ] = coef * ( fine( IDF(ixf, iyf, izf))   + fine( IDF(ixf+1, iyf, izf)  ) + fine( IDF(ixf+1, iyf+1, izf)  ) + fine( IDF(ixf, iyf+1, izf)) +
                                       fine( IDF(ixf, iyf, izf+1)) + fine( IDF(ixf+1, iyf, izf+1)) + fine( IDF(ixf+1, iyf+1, izf+1)) + fine( IDF(ixf, iyf+1, izf+1)) ) ;

}
#undef coef
#undef tx
#undef ty
#undef tz
#undef IDC
#undef IDF
#undef fine
#undef ixc
#undef iyc
#undef izc
#undef ixf
#undef iyf
#undef izf

