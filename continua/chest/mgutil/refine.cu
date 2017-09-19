 __global__ void Refine_Cuda_3D(__CUFLOAT *fine, 
#ifndef __MGTEX
 __CUFLOAT *coarse, 
#endif
 const __CINT i3f, const __CINT i3c, const __CINT nnx, const __CINT nny, const __CINT nnz) {
// interpolation coefficient in 3D is 1/64
#define coef 0.015625000000000000
#define cb   0.666666666666666666
//constants
#define three 3.0
#define nine 9.0
#define twentyseven 27.0


#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z
#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z
#define nbx gridDim.x
#define nby gridDim.y
#define nbz gridDim.z
#define ntx _BRFNE_X
#define nty _BRFNE_Y
#define ntz _BRFNE_Z

//
// load coarse points into local tiles, similarly in spirit to GS kernel
//
#define __VOLATILE
#define tilesize ( ( _BRFNE_X + 2 ) * ( _BRFNE_Y + 2 ) )
 __VOLATILE __shared__ float  clocal[ 2 * tilesize ];
 float *cfront = clocal ;
 float *cback  = cfront + tilesize ;
 float *cswap ;
 float cwsf, cesf, cwnf, cenf, cwsb, cesb, cwnb, cenb ; // cell values places in registers

 unsigned int ixc = (bx * ntx) + tx;
 unsigned int iyc = (by * nty) + ty;
 bool qedgeblock;
#define ixf (ixc<<1)
#define iyf (iyc<<1)

//define index functions
#define IDC(i,j,k)  ( i3c - 1 + (k)*(nnx+2)*(nny+2)   + (j)*(nnx+2)   + (i) )
//note that the nnx/y/z are the coarse lengths; hence factors of two below
#define IDF(i,j,k)  ( i3f - 1 + (k)*4*(nnx+1)*(nny+1) + (j)*2*(nnx+1) + (i) )
//local memory:
#define IOL(i,j)  ((j)*(_BRFNE_X+2) + (i))
//
#define cfront(i,j) cfront[IOL(i,j)]
#define cback(i,j) cback[IOL(i,j)]
#define fine(i,j,k) fine[IDF(i,j,k)]

#ifdef __MGTEX
#define coarse(i,j,k) tex1Dfetch(texcoarse,IDC(i,j,k));
#else 
#define coarse(i,j,k) coarse[IDC(i,j,k)]
#endif
// unsigned int ind, indl;
// loop over all z-slices, including first row (k=-1)
//
 for (int k=-1 ; k<=nnz ; k++) { 
// load next coarse row
  cback(tx+1,ty+1)=coarse(ixc+1,iyc+1,k+1);
// load additional boundary points
  if (ty < 2) { // BC in y-direction
   cback( tx+1, (nty+1)*ty) = coarse(ixc+1, iyc + (nty+1)*ty, k+1);
  }
  if (tx < 2) { // BC in x-direction
   cback( (ntx+1)*(tx%2), ty+1) = coarse(ixc + (ntx+1)*(tx%2), iyc+1, k+1);
  }

// populate corners either by interpolation (boundary corners) or by extra reads from device
  qedgeblock = ( (bx==0 || bx == nbx-1) && (by==0 || by == nby-1)  ) ;
  if  (qedgeblock) __syncthreads(); // make sure tiles are populated prior to interpolation

  if ( ((tx==0) || (tx==ntx-1) || (ixc==nnx-1) ) && ((ty==0) || (ty==nty-1) || (iyc==nny-1))) {
#define dx  (1-2*(bool)(tx))
#define dy  (1-2*(bool)(tx))
   cback(tx+1-dx,ty+1-dy)=coarse(ixc+1-dx, iyc+1-dy, k+1) ; // read from device (possibly wrong values)
#undef dx
#undef dy
   if (qedgeblock) { // obtain values at grid edge by interpolation
#define dx  (1-2*(bool)(ixc))
#define dy  (1-2*(bool)(iyc))
    cback(tx+1-dx,ty+1-dy) = cback(tx+1-dx,ty+1) + cback(tx+1,ty+1-dy) - cback(tx+1,ty+1) ; // treat four possible corner points togetherOR
//    cback(0,0) = cback(1,0) + cback(0,1) - cback(1,1)
   }
  }
// explicitly treat eight corners here (k=0,nnz) case
  if (k==0) {
   __syncthreads(); // need this because further interpolation below :
   if ( ((ixc==0) || (ixc==nnx-1) ) && ((iyc==0) || (iyc==nny-1)) ) {
    cfront(tx+1-dx,ty+1-dy) = cb * ( cfront(tx+1-dx,ty+1) + cfront(tx+1,ty+1-dy) + cback(tx+1-dx,ty+1-dy) ) - cback(tx+1,ty+1) ; // treat four possible corner points together
   }
  } else if (k==nnz) {
   __syncthreads(); // need this because further interpolation below :
   if ( ((ixc==0) || (ixc==nnx-1) ) && ((iyc==0) || (iyc==nny-1)) ) {
    cback(tx+1-dx,ty+1-dy) = cb * ( cback(tx+1-dx,ty+1) + cback(tx+1,ty+1-dy) + cfront(tx+1-dx,ty+1-dy) ) - cfront(tx+1,ty+1) ; // treat four possible corner points together
#undef dx
#undef dy
   }
  }
//
// compute fine grid data by interpolation
  __syncthreads() ; // make sure shared data is up-to-date
//
  cwsf=cfront(tx,ty) ;
  cesf=cfront(tx+1,ty) ;
  cenf=cfront(tx+1,ty+1) ;
  cwnf=cfront(tx,ty+1) ;
//
  cwsb=cback(tx,ty) ;
  cesb=cback(tx+1,ty) ;
  cenb=cback(tx+1,ty+1) ;
  cwnb=cback(tx,ty+1) ;
//
// swap cfront and cback local shmem pointers
//
  cswap =cfront;
  cfront=cback;
  cback =cswap;
  if (k<0) continue; // skip to next iteration
//
// perform the actual interpolation, uploading directly to device :
#define izf 2*k
  fine(ixf+1,iyf+1,izf+1) += (cwsf + (cesf + cwnf + cwsb)*three + (cenf + cesb + cwnb)*nine + cenb*twentyseven);
  fine(ixf,  iyf+1,izf+1) += (cesf + (cwsf + cenf + cesb)*three + (cwnf + cwsb + cenb)*nine + cwnb*twentyseven);
  fine(ixf+1,iyf,  izf+1) += (cwnf + (cenf + cwsf + cwnb)*three + (cesf + cenb + cwsb)*nine + cesb*twentyseven);
  fine(ixf,  iyf,  izf+1) += (cenf + (cwnf + cesf + cenb)*three + (cwsf + cwnb + cesb)*nine + cwsb*twentyseven);
//
  fine(ixf+1,iyf+1,izf)   += (cwsb + (cesb + cwnb + cwsf)*three + (cenb + cesf + cwnf)*nine + cenf*twentyseven);
  fine(ixf,  iyf+1,izf)   += (cesb + (cwsb + cenb + cesf)*three + (cwnb + cwsf + cenf)*nine + cwnf*twentyseven);
  fine(ixf+1,iyf,  izf)   += (cwnb + (cenb + cwsb + cwnf)*three + (cesb + cenf + cwsf)*nine + cesf*twentyseven);
  fine(ixf,  iyf,  izf)   += (cenb + (cwnb + cesb + cenf)*three + (cwsb + cwnf + cesf)*nine + cwsf*twentyseven);
//
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
#undef izf
#undef coarse
#undef fine
#undef tilesize
#undef dx
#undef dy
#undef nbx
#undef nby
#undef nbz
#undef ntx
#undef nty
#undef ntz
