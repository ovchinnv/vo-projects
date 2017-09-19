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
//
// load coarse points into local tiles, similarly in spirit to GS kernel
//
#define __VOLATILE
#define tilesize ( ( _BRFNE_X + 2 ) * ( _BRFNE_Y + 2 ) )
 __VOLATILE __shared__ float  clocal[ 2 * tilesize ];
 float *cfront = clocal ;
 float *cback  = cfront + tilesize ;

 unsigned int ixc = (blockIdx.x * _BRFNE_X) + tx;
 unsigned int iyc = (blockIdx.y * _BRFNE_Y) + ty;
#define ixf 2*ixc
#define iyf 2*iyc

//define index functions
#define IDC(i,j,k)  ( i3c - 1 + (k)*(nnx+2)*(nny+2)   + (j)*(nnx+2)   + (i) )
//note that the nnx/y/z are the coarse lengths; hence factors of two below
#define IDF(i,j,k)  ( i3f - 1 + (k)*4*(nnx+1)*(nny+2) + (j)*2*(nnx+1) + (i) )
//local memory:
#define IOL(i,j)  ((j)*(_BRFNE_X+2) + (i))
//
// unsigned int ind, indl;
// load first coarse row
 cfront[IOL(tx+1,ty+1)]=coarse(IDC(ixc+tx+1,iyc+1,0));
//
// loop over all z-slices
//
 for (unsigned int k=1 ; k<nnz-1 ; k++) {
// load next coarse row
  cback[IOL(tx+1,ty+1)]=coarse(IDC(ixc+1,iyc+1,k));
// load boundary points

// populate coarse ghost points by interpolation

// compute fine grid data by interpolation
  __syncthreads() ; // make sure shared data is up-to-date


// swap cfront and cback
  cfront+=tilesize;
  if (cfront!=cback) {
   cfront=cback;
   cback+=tilesize;
  } else {
   cback-=tilesize; // move cback pointer up
  }

 }

}


#undef tx
#undef ty
#undef __VOLATILE
#undef IOL
#undef IDC
#undef IDF
#undef ixf
#undef iyf

