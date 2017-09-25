// NOTE : cannot use textures because input and output are stored in the same array
 __global__ void Refine_Cuda_3D(__CUFLOAT *fine, 
//#ifndef __MGTEX
 __CUFLOAT *coarse, 
//#endif
 const __CINT i3f, const __CINT i3c, const __CINT nnx, const __CINT nny, const __CINT nnz) {
// interpolation coefficient in 3D is 1/64
#define coef 0.0156250000000
#define cb   0.6666666666666
//constants

#define three        3.0
#define nine         9.0
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
#define tilesize ( ( ntx + 2 ) * ( nty + 2 ) )
 __VOLATILE __shared__ float  clocal[ 2 * tilesize ];
 float *cfront = clocal ;
 float *cback  = cfront + tilesize ;
 float *cswap ;
// float cwsf, cesf, cwnf, cenf, cwsb, cesb, cwnb, cenb ; // cell values places in registers
 float cwsf, cesf, cwnf, cenf, cwsb, cesb, cwnb, cenb ; // cell values places in registers

 unsigned int ixc = (bx * ntx) + tx;
 unsigned int iyc = (by * nty) + ty;
 bool qcorthread;
#define ixf (ixc<<1)
#define iyf (iyc<<1)

//define index functions
#define IDC(i,j,k)  ( i3c - 1 + (k)*(nnx+2)*(nny+2)   + (j)*(nnx+2)   + (i) )
//note that the nnx/y/z are the coarse lengths; hence factors of two below
#define IDF(i,j,k)  ( i3f - 1 + (k)*4*(nnx+1)*(nny+1) + (j)*2*(nnx+1) + (i) )
//local memory:
#define IOL(i,j)  ((j)*(ntx+2) + (i))
//
#define cfront(i,j) cfront[IOL(i,j)]
#define cback(i,j) cback[IOL(i,j)]
#define fine(i,j,k) fine[IDF(i,j,k)]

//#ifdef __MGTEX
//#define coarse(i,j,k) tex1Dfetch(texcoarse,IDC(i,j,k))
//#else
#define coarse(i,j,k) coarse[IDC(i,j,k)]
//#endif
// unsigned int ind, indl;
// loop over all z-slices, including first row (k=-1)
//
 for (int k=-1 ; k<=nnz ; k++) {
//  __syncthreads(); // should not be needed
// load next coarse row
  cback(tx+1,ty+1)=coef*coarse(ixc+1,iyc+1,k+1);
// load additional boundary points
  if ( (ty==0) || ((ty==nty-1)&&(iyc==nny-1)) ) { // BC in y ; both boundaries at far edge only
#define dyl  (1-2*(bool)(ty))
#define dxl  (1-2*(bool)(tx))
   cback( tx+1, ty+1-dyl ) = coef*coarse(ixc+1, iyc+1-dyl, k+1) ;
  }
  if ( (tx==0) || ((tx==ntx-1)&&(ixc==nnx-1)) ) { // BC in x
   cback( tx+1-dxl, ty+1 ) = coef*coarse(ixc+1-dxl, iyc+1, k+1) ;
// corners :
   if ( (ty==0) || ((ty==nty-1)&&(iyc==nny-1)) ) {
    cback(tx+1-dxl,ty+1-dyl) = coef*coarse(ixc+1-dxl, iyc+1-dyl, k+1) ;
   }
  }
// flag corner threads
  qcorthread = ( (ixc==0) || (ixc==nnx-1) ) && ( (iyc==0) || (iyc==nny-1) ) ;
#define dxg  (1-2*(bool)(ixc))
#define dyg  (1-2*(bool)(iyc))
// populate boundary corners by interpolation
  if (k==0) {
    __syncthreads();
// interpolate k=0 corner "bc" arrays
// y-const.
    if ( (ixc==0) || (ixc==nnx-1) ) {
     cfront(tx+1-dxg,ty+1) = cback(tx+1-dxg,ty+1) + cfront(tx+1,ty+1) - cback(tx+1,ty+1) ; // treat four possible corner points together
     if ( (ty==0) || ((ty==nty-1)&&(iyc==nny-1)) ) { // BC in y
      cfront(tx+1-dxg,ty+1-dyl) = cback(tx+1-dxg,ty+1-dyl) + cfront(tx+1,ty+1-dyl) - cback(tx+1,ty+1-dyl) ; // treat four possible corner points together 
     }
    }
// x-const
    if ( (iyc==0) || (iyc==nny-1) ) {
     cfront(tx+1,ty+1-dyg) = cback(tx+1,ty+1-dyg) + cfront(tx+1,ty+1) - cback(tx+1,ty+1) ; // treat four possible corner points together
     if ( (tx==0) || ((tx==ntx-1)&&(ixc==nnx-1)) ) { // BC in x
      cfront(tx+1-dxl,ty+1-dyg) = cback(tx+1-dxl,ty+1-dyg) + cfront(tx+1-dxl,ty+1) - cback(tx+1-dxl,ty+1) ; // treat four possible corner points together
// NOTE : corners should be taken care of below because they are the outermost corners
     }
    }
  } else if (k==nnz) {
    __syncthreads();
// interpolate k=nnz corner "bc" arrays
// y-const.
    if ( (ixc==0) || (ixc==nnx-1) ) {
     cback(tx+1-dxg,ty+1) = cfront(tx+1-dxg,ty+1) + cback(tx+1,ty+1) - cfront(tx+1,ty+1) ; // treat four possible corner points together
     if ( (ty==0) || ((ty==nty-1)&&(iyc==nny-1)) ) { // BC in y
      cback(tx+1-dxg,ty+1-dyl) = cfront(tx+1-dxg,ty+1-dyl) + cback(tx+1,ty+1-dyl) - cfront(tx+1,ty+1-dyl) ; // treat four possible corner points together 
     }
    }
// x-const
    if ( (iyc==0) || (iyc==nny-1) ) {
     cback(tx+1,ty+1-dyg) = cfront(tx+1,ty+1-dyg) + cback(tx+1,ty+1) - cfront(tx+1,ty+1) ; // treat four possible corner points together
     if ( (tx==0) || ((tx==ntx-1)&&(ixc==nnx-1)) ) { // BC in x
      cback(tx+1-dxl,ty+1-dyg) = cfront(tx+1-dxl,ty+1-dyg) + cback(tx+1-dxl,ty+1) - cfront(tx+1-dxl,ty+1) ; // treat four possible corner points together
     }
    }
  }

// populate boundary corners by interpolation
  if ( ( (bx==0) || (bx==nbx-1) ) && ( (by==0) || (by==nby-1) ) ) { // whether this is a corner block
   __syncthreads();
//
   if  ( qcorthread ) {
    cback(tx+1-dxg,ty+1-dyg) = cback(tx+1-dxg,ty+1) + cback(tx+1,ty+1-dyg) - cback(tx+1,ty+1) ; // treat four possible corner points together
   }
// explicitly treat four corners here (k=0) case
   if (k==0) {
    __syncthreads();
    if ( qcorthread ) {
     cfront(tx+1-dxg,ty+1-dyg) = cb * ( cfront(tx+1-dxg,ty+1) + cfront(tx+1,ty+1-dyg) + cback(tx+1-dxg,ty+1-dyg) ) - cback(tx+1,ty+1) ; // treat four possible corner points together
    }
   } else if (k==nnz) {
// four corners at k=nnz
    __syncthreads();
    if (qcorthread) {
     cback(tx+1-dxg,ty+1-dyg) = cb * ( cback(tx+1-dxg,ty+1) + cback(tx+1,ty+1-dyg) + cfront(tx+1-dxg,ty+1-dyg) ) - cfront(tx+1,ty+1) ; // treat four possible corner points together
#undef dx
#undef dy
    }
   }
  }
//
// compute fine grid data by interpolation
//
  if (k>-1) { // only for positive k
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
// perform the actual interpolation, uploading directly to device :
#define izf 2*k
//
   fine(ixf+1,iyf+1,izf)   += (cwsb + (cesb + cwnb + cwsf)*three + (cenb + cesf + cwnf)*nine + cenf*twentyseven);
   fine(ixf,  iyf+1,izf)   += (cesb + (cwsb + cenb + cesf)*three + (cwnb + cwsf + cenf)*nine + cwnf*twentyseven);
   fine(ixf+1,iyf,  izf)   += (cwnb + (cenb + cwsb + cwnf)*three + (cesb + cenf + cwsf)*nine + cesf*twentyseven);
   fine(ixf,  iyf,  izf)   += (cenb + (cwnb + cesb + cenf)*three + (cwsb + cwnf + cesf)*nine + cwsf*twentyseven);
//
   fine(ixf+1,iyf+1,izf+1) += (cwsf + (cesf + cwnf + cwsb)*three + (cenf + cesb + cwnb)*nine + cenb*twentyseven);
   fine(ixf,  iyf+1,izf+1) += (cesf + (cwsf + cenf + cesb)*three + (cwnf + cwsb + cenb)*nine + cwnb*twentyseven);
   fine(ixf+1,iyf,  izf+1) += (cwnf + (cenf + cwsf + cwnb)*three + (cesf + cenb + cwsb)*nine + cesb*twentyseven);
   fine(ixf,  iyf,  izf+1) += (cenf + (cwnf + cesf + cenb)*three + (cwsf + cwnb + cesb)*nine + cwsb*twentyseven);
//
//additional far boundary points
   if (ixc==nnx-1) {
    cwsf=cfront(tx+1,ty) ;
    cesf=cfront(tx+2,ty) ;
    cenf=cfront(tx+2,ty+1) ;
    cwnf=cfront(tx+1,ty+1) ;
//
    cwsb=cback(tx+1,ty) ;
    cesb=cback(tx+2,ty) ;
    cenb=cback(tx+2,ty+1) ;
    cwnb=cback(tx+1,ty+1) ;

    fine(ixf+3,iyf+1,izf)   += (cwsb + (cesb + cwnb + cwsf)*three + (cenb + cesf + cwnf)*nine + cenf*twentyseven); // optional bc point ; can be omitted
    fine(ixf+2,iyf+1,izf)   += (cesb + (cwsb + cenb + cesf)*three + (cwnb + cwsf + cenf)*nine + cwnf*twentyseven);
    fine(ixf+3,iyf,  izf)   += (cwnb + (cenb + cwsb + cwnf)*three + (cesb + cenf + cwsf)*nine + cesf*twentyseven); // omit
    fine(ixf+2,iyf,  izf)   += (cenb + (cwnb + cesb + cenf)*three + (cwsb + cwnf + cesf)*nine + cwsf*twentyseven);
//
    fine(ixf+3,iyf+1,izf+1) += (cwsf + (cesf + cwnf + cwsb)*three + (cenf + cesb + cwnb)*nine + cenb*twentyseven); // omit
    fine(ixf+2,iyf+1,izf+1) += (cesf + (cwsf + cenf + cesb)*three + (cwnf + cwsb + cenb)*nine + cwnb*twentyseven);
    fine(ixf+3,iyf,  izf+1) += (cwnf + (cenf + cwsf + cwnb)*three + (cesf + cenb + cwsb)*nine + cesb*twentyseven); // omit
    fine(ixf+2,iyf,  izf+1) += (cenf + (cwnf + cesf + cenb)*three + (cwsf + cwnb + cesb)*nine + cwsb*twentyseven);
   } //ixc
//
   if (iyc==nny-1) {
    cwsf=cfront(tx,ty+1) ;
    cesf=cfront(tx+1,ty+1) ;
    cenf=cfront(tx+1,ty+2) ;
    cwnf=cfront(tx,ty+2) ;
//
    cwsb=cback(tx,ty+1) ;
    cesb=cback(tx+1,ty+1) ;
    cenb=cback(tx+1,ty+2) ;
    cwnb=cback(tx,ty+2) ;
//
    fine(ixf+1,iyf+3,izf)   += (cwsb + (cesb + cwnb + cwsf)*three + (cenb + cesf + cwnf)*nine + cenf*twentyseven); // omit optional bc point
    fine(ixf,  iyf+3,izf)   += (cesb + (cwsb + cenb + cesf)*three + (cwnb + cwsf + cenf)*nine + cwnf*twentyseven); //omit
    fine(ixf+1,iyf+2,izf)   += (cwnb + (cenb + cwsb + cwnf)*three + (cesb + cenf + cwsf)*nine + cesf*twentyseven);
    fine(ixf,  iyf+2,izf)   += (cenb + (cwnb + cesb + cenf)*three + (cwsb + cwnf + cesf)*nine + cwsf*twentyseven);
//
    fine(ixf+1,iyf+3,izf+1) += (cwsf + (cesf + cwnf + cwsb)*three + (cenf + cesb + cwnb)*nine + cenb*twentyseven); //omit
    fine(ixf,  iyf+3,izf+1) += (cesf + (cwsf + cenf + cesb)*three + (cwnf + cwsb + cenb)*nine + cwnb*twentyseven); //omit
    fine(ixf+1,iyf+2,izf+1) += (cwnf + (cenf + cwsf + cwnb)*three + (cesf + cenb + cwsb)*nine + cesb*twentyseven);
    fine(ixf,  iyf+2,izf+1) += (cenf + (cwnf + cesf + cenb)*three + (cwsf + cwnb + cesb)*nine + cwsb*twentyseven);
// x/y corner
    if (ixc==nnx-1) {
     cwsf=cfront(tx+1,ty+1) ;
     cesf=cfront(tx+2,ty+1) ;
     cenf=cfront(tx+2,ty+2) ;
     cwnf=cfront(tx+1,ty+2) ;
//
     cwsb=cback(tx+1,ty+1) ;
     cesb=cback(tx+2,ty+1) ;
     cenb=cback(tx+2,ty+2) ;
     cwnb=cback(tx+1,ty+2) ;
//
     fine(ixf+3,iyf+3,izf)   += (cwsb + (cesb + cwnb + cwsf)*three + (cenb + cesf + cwnf)*nine + cenf*twentyseven); // omit
     fine(ixf+2,iyf+3,izf)   += (cesb + (cwsb + cenb + cesf)*three + (cwnb + cwsf + cenf)*nine + cwnf*twentyseven); // omit
     fine(ixf+3,iyf+2,izf)   += (cwnb + (cenb + cwsb + cwnf)*three + (cesb + cenf + cwsf)*nine + cesf*twentyseven); // omit
     fine(ixf+2,iyf+2,izf)   += (cenb + (cwnb + cesb + cenf)*three + (cwsb + cwnf + cesf)*nine + cwsf*twentyseven);
//
     fine(ixf+3,iyf+3,izf+1) += (cwsf + (cesf + cwnf + cwsb)*three + (cenf + cesb + cwnb)*nine + cenb*twentyseven); // omit
     fine(ixf+2,iyf+3,izf+1) += (cesf + (cwsf + cenf + cesb)*three + (cwnf + cwsb + cenb)*nine + cwnb*twentyseven); // omit
     fine(ixf+3,iyf+2,izf+1) += (cwnf + (cenf + cwsf + cwnb)*three + (cesf + cenb + cwsb)*nine + cesb*twentyseven); // omit
     fine(ixf+2,iyf+2,izf+1) += (cenf + (cwnf + cesf + cenb)*three + (cwsf + cwnb + cesb)*nine + cwsb*twentyseven);
    } //ixc
   } //iyc
  } // k>=0
//
// swap cfront and cback local shmem pointers
//
  cswap =cfront;
  cfront=cback;
  cback =cswap;
 } //k
} //refine

#undef tx
#undef ty
#undef bx
#undef by
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
