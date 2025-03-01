 __global__ void Residual_Cuda_3D(__CUFLOAT *devres, 
#ifndef __MGTEX
 __CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, 
#endif
 __CUFLOAT *devoodx, __CUFLOAT *devoody, __CUFLOAT *devoodz,
 const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz, const __CINT qmaxres, const __CINT qresnorm,
 __CUFLOAT *maxres, __CINT *imax) {
//based on the GS kernel
//NOTE : I am loading two adjacent tiles (in red/black fashion, as for GS, but this is not necessary
// because there is no particular update order that needs to be enforced ; might rewrite kernel to use one tile

#ifdef __DSHMEM
#define ntx blockDim.x
#define nty blockDim.y
#else
#define ntx _BGSMAX_X
#define nty _BGSMAX_Y
#endif
// thread indices into device memory:
#define tx threadIdx.x
#define ty threadIdx.y

#ifdef __MGTEX
//use texture addressing
#define devp(_IND) tex1Dfetch(texp,_IND)
#define deveps(_IND) tex1Dfetch(texeps,_IND)
#define devrhs(_IND) tex1Dfetch(texrhs,_IND)
#define devkappa(_IND) tex1Dfetch(texkappa,_IND)
#else
#define devp(_IND) devp[_IND]
#define deveps(_IND) deveps[_IND]
#define devrhs(_IND) devrhs[_IND]
#define devkappa(_IND) devkappa[_IND]
#endif

#define sizep      (( _SX * ntx + 2 ) * ( _SY * nty + 2 ))
#define sizeeps    sizep
#define sizeres    (( _SX * ntx     ) * ( _SY * nty     ))
#define sizedxcen  (_SX * ntx + 1)
#define sizedxcor  (_SX * ntx)
#define sizedycen  (_SY * nty + 1)
#define sizedycor  (_SY * nty)

// 2D local sh/mem arrays
#define __VOLATILE
#ifndef __DSHMEM
 __VOLATILE __shared__ float p       [ sizep ];
 __VOLATILE __shared__ float eps     [ sizeeps ];
 __VOLATILE __shared__ float res     [ sizeres ]; // for computing maxumum residual
 __VOLATILE __shared__ __CINT resind [ sizeres ]; // for computing index location of maxumum residual
 __shared__ float oodxcen [ sizedxcen ];
 __shared__ float oodxcor [ sizedxcor ];
 __shared__ float oodycen [ sizedycen ];
 __shared__ float oodycor [ sizedycor ];

#define p(i)       p[i]
#define eps(i)     eps[i]
#define res(i)     res[i]
#define resind(i)  resind[i]
#define oodxcen(i) oodxcen[i]
#define oodycen(i) oodycen[i]
#define oodxcor(i) oodxcor[i]
#define oodycor(i) oodycor[i]

#else
#define __VOLATILE
 extern __VOLATILE __shared__ float local[];
// array offsets below
/* does this waste registers ?
 float * p =   local ;
 float * eps = p + sizep ;
 float * res = eps + sizeeps ;
 float * resind = res + sizeres ;
 float * oodxcen = resind + sizeres ;
 float * oodxcor = oodxcen + sizedxcen ;
 float * oodycen = oodxcor + sizedxcor ;
 float * oodycor = oodycen + sizedycen ;
/*/
#define p(i)        local[i]
#define eps(i)      local[i+sizep]
#define res(i)      local[i+sizep+sizeeps]
#define resind(i)   local[i+sizep+sizeeps+sizeres]
#define oodxcen(i)  local[i+sizep+sizeeps+sizeres+sizeres]
#define oodxcor(i)  local[i+sizep+sizeeps+sizeres+sizeres+sizedxcen]
#define oodycen(i)  local[i+sizep+sizeeps+sizeres+sizeres+sizedxcen+sizedxcor]
#define oodycor(i)  local[i+sizep+sizeeps+sizeres+sizeres+sizedxcen+sizedxcor+sizedycen]
//*/
#endif


// registers -- 16
// NOTE : tiles are shfted from "front" to "back"
 float pback[2], pfront[2], pcur[2]; // solution z-points
 float eback[2], efront[2], ecur[2]; // epsilon z-points
 float rhs[2], kappa[2] ;
// 1D local shmem metrics
//registers
 float oodzcor, oodzcenfront, oodzcenback ;
// NOTE that the code was written with the assumption the _SX=2, _SY=1 ; other values may not work yet
// registers 15 + 2=17
 unsigned int ix = _SX*(blockIdx.x * ntx) + tx;
 unsigned int iy = _SY*(blockIdx.y * nty) + ty;
//registers
 float w, e, s, n, b, f, o ;// FD stencil coefficients
//registers
 unsigned int ind, indl ;// indices
// index maps : NOTE that they are zero-based (not one-based as in FORTRAN)
// global memory:
#define IOG(i,j,k)  ( i3b - 1 + (k)*(nx*ny)       + (j)*(nx)   + (i) )
#define IIG(i,j,k)  ( i3  - 1 + (k)*(nx-2)*(ny-2) + (j)*(nx-2) + (i) )
//local memory:
#define IOL(i,j)  ((j)*(_SX*ntx+2) + (i))
#define IIL(i,j)  ((j)*(_SX*ntx)   + (i))

 ind=IOG(ix+tx+1,iy+1,0) ;
 pfront[0]=devp(ind)   ; efront[0]=deveps(ind);       // red

 ind++;
 pfront[1]=devp(ind)   ; efront[1]=deveps(ind);       // black

 ind=IOG(ix+tx+1,iy+1,1) ;
 pcur[0]=devp(ind) ; ecur[0]=deveps(ind);

 ind++;
 pcur[1]=devp(ind) ; ecur[1]=deveps(ind);

// read metrics
// x-metrics
 if (ty < _SX) {
  oodxcen(tx+ntx*ty+(ty>0) ) = devoodx[i1 - 1 + ix + ntx*ty+(ty>0)];
  oodxcor(tx+ntx*ty)         = devoodx[i1 - 1 + (nx-1) + ix + ntx*ty];
 }
// NOTE that if we know that tx=ty we can use a single if clause to read both metrics
// y-metrics
 if (tx < _SY) {
  oodycen(ty+nty*tx+(tx>0)) = devoody[j1 - 1 + iy + nty*tx+(tx>0)];
  oodycor(ty+nty*tx)        = devoody[j1 - 1 + ny-1 + iy + nty*tx];
 }
//
// z - metric
 oodzcenfront=devoodz[k1 - 1];
 __syncthreads(); // need dxcor/dycor updated
//
 if ( !( (bool)ty || (bool)tx ) ) { // first (0th) thread reads
// deal with y-metric first, since it will usually require global mem read.
  if (_SY<2) { // read from gmem
   oodycen(nty) = devoody[j1 - 1 + (iy - ty) + nty];
  } else {
   oodycen(nty) = 2.0 / ( 1.0/oodycor(nty-1)+1.0/oodycor(nty-(_SY<2)) ); // ( _SY<2 ) avoids static memory check warnings
  }
   oodxcen(ntx) = devoodx[i1 - 1 + (ix - tx) + ntx];
 }
//
// initialize local maximum residual arrays
 if (qmaxres) {
  indl=IIL(2*tx,ty);
  res(indl)  =-INFINITY; resind(indl)=-1;
  res(++indl)=-INFINITY; resind(indl)=-1;
 }
//
// loop over all z-slices
 for (unsigned int k=1 ; k<nz-1 ; k++) {
  __syncthreads(); // This is necessary but I do not know why, since there are no thread-dependend memory access reads between here and the next syncthreads
// might not be optimal due to striped access:
// populate slice
  ind=IOL(2*tx+1,ty+1);
  p   (ind) = pcur[0];
  eps (ind) = ecur[0];
  ind++;
  p   (ind) = pcur[1];
  eps (ind) = ecur[1];
// update z-metric (registers) ; might convert to shmem
  oodzcenback=devoodz[k1 - 1 + k];
  oodzcor=devoodz[k1 - 1 + k - 1 + (nz-1) ];
// load ghost points
  if (ty < 4) { // BC in y direction
   p  (IOL(tx+1+ntx*(ty%2), (ty/2)*(nty+1))) = devp  (IOG(ix+1+ntx*(ty%2), iy-ty+(ty/2)*(nty+1), k)); // subtract ty to get first tile coordinate
   eps(IOL(tx+1+ntx*(ty%2), (ty/2)*(nty+1))) = deveps(IOG(ix+1+ntx*(ty%2), iy-ty+(ty/2)*(nty+1), k));
  }
  if (tx < 2) { // BC in x-direction
   p  (IOL((_SX*ntx+1)*(tx%2), ty+1)) = devp  (IOG(ix-tx+(_SX*ntx+1)*(tx%2), iy+1, k));
   eps(IOL((_SX*ntx+1)*(tx%2), ty+1)) = deveps(IOG(ix-tx+(_SX*ntx+1)*(tx%2), iy+1, k));
  }

// load next p & eps slices
  ind=IOG(ix+tx+1,iy+1,k+1) ;
  pback[0] = devp(ind) ; eback[0] = deveps(ind);
  ind++;
  pback[1] = devp(ind) ; eback[1] = deveps(ind);

// compute residual for two adjacent points
// NOTE that the two residual points are not interdependent (unlike in the p update case)
// thus, we can define DX in arbitrary order (i.e. not in alternating fashion, as required for GS/RB update)
#//define _DX (ty%2)
#define _DX 0
  ind=IIG(ix+tx+_DX,iy,k-1) ; // inner global index
  indl=IOL(2*tx+1+_DX,ty+1) ; // outer local index -- used for p
  __syncthreads();

  //prefetch
  rhs[_DX]    =devrhs(ind);
  kappa[_DX]  =devkappa(ind);
  rhs[1-_DX]  =devrhs(ind+1-_DX);
  kappa[1-_DX]=devkappa(ind+1-_DX);

  w = ( eps(indl-1) + ecur[_DX] ) * ( oodxcor(2*tx+_DX) * oodxcen(2*tx+_DX) );
  e = ( eps(indl+1) + ecur[_DX] ) * ( oodxcor(2*tx+_DX) * oodxcen(2*tx+1+_DX) );
  s = ( eps(indl-(_SX*ntx+2)) + ecur[_DX] ) * ( oodycor(ty) * oodycen(ty) );
  n = ( eps(indl+(_SX*ntx+2)) + ecur[_DX] ) * ( oodycor(ty) * oodycen(ty+1) );
  f = ( efront[_DX] + ecur[_DX] )              * ( oodzcor * oodzcenfront );
  b = ( eback[_DX]  + ecur[_DX] )              * ( oodzcor * oodzcenback );
  o = - ( w + e + s + n + b + f ) + 2.0 * kappa[_DX];
// note that o (and other metrics) are twice what they should be ; this is accounted for explicitly below
// reuse register e
#define residual e
//RHS o-weighted
// e =        - 0.5 * ( o * ( pcur[_DX] - devrhs[ind] ) + w*p(indl-1) + e*p(indl+1) + s*p(indl-(_SX*ntx+2)) + n*p(indl+(_SX*ntx+2)) + f*pfront[_DX] + b*pback[_DX] );
//RHS not weighted
  residual = rhs[_DX] - 0.5 * ( o * pcur[_DX] + w*p(indl-1) + e*p(indl+1) + s*p(indl-(_SX*ntx+2)) + n*p(indl+(_SX*ntx+2)) + f*pfront[_DX] + b*pback[_DX] );
  devres[ind-i3+1] = residual; // upload to device ; note that the residual array always starts at zero offset ; hence the subtraction of i3
// -- update max residual and index ; note that the index stored is the inner global index (per mg_resout.src)
#define __UPDATE_RES \
  if (qmaxres) { \
   if (qresnorm) { \
    residual=fabs(2.0f*residual/(o)); \
   } else { \
    residual=fabs(residual); \
   } \
   if (res(IIL(2*tx+_DX,ty)) < residual) { \
    res   (IIL(2*tx+_DX,ty)) = residual ;\
    resind(IIL(2*tx+_DX,ty)) = ind - i3 + 1 ; \
   } \
  }

  __UPDATE_RES

#undef _DX
#//define _DX (1-ty%2)
#define _DX 1

  ind=IIG(ix+tx+_DX,iy,k-1) ; // inner global index
  indl=IOL(2*tx+1+_DX,ty+1) ;

  w = ( eps(indl-1) + ecur[_DX] ) * ( oodxcor(2*tx+_DX) * oodxcen(2*tx+_DX) );
  e = ( eps(indl+1) + ecur[_DX] ) * ( oodxcor(2*tx+_DX) * oodxcen(2*tx+1+_DX) );
  s = ( eps(indl-(_SX*ntx+2)) + ecur[_DX] ) * ( oodycor(ty) * oodycen(ty) );
  n = ( eps(indl+(_SX*ntx+2)) + ecur[_DX] ) * ( oodycor(ty) * oodycen(ty+1) );
  f = ( efront[_DX] + ecur[_DX] )              * ( oodzcor * oodzcenfront );
  b = ( eback[_DX]  + ecur[_DX] )              * ( oodzcor * oodzcenback );
  o = - ( w + e + s + n + b + f ) + 2.0 * kappa[_DX];
// note that o (and other metrics) are twice what they should be ; this is accounted for explicitly below
// reuse register e
#define residual e
//RHS o-weighted
// e =        - 0.5 * ( o * ( pcur[_DX] - devrhs[ind] ) + w*p(indl-1) + e*p(indl+1) + s*p(indl-(_SX*ntx+2)) + n*p(indl+(_SX*ntx+2)) + f*pfront[_DX] + b*pback[_DX] );
//RHS not weighted
  residual = rhs[_DX] - 0.5 * ( o * pcur[_DX] + w*p(indl-1) + e*p(indl+1) + s*p(indl-(_SX*ntx+2)) + n*p(indl+(_SX*ntx+2)) + f*pfront[_DX] + b*pback[_DX] );
  devres[ind-i3+1] = residual; // upload to device ; note that the residual array always starts at zero offset ; hence the subtraction of i3
// -- update max residual
  __UPDATE_RES
// move forward in z-direction
  pfront[0]=pcur[0];  efront[0]=ecur[0];
  pcur[0]  =pback[0]; ecur[0]  =eback[0];
  pfront[1]=pcur[1];  efront[1]=ecur[1];
  pcur[1]  =pback[1]; ecur[1]  =eback[1];
  oodzcenfront=oodzcenback;
 } // for
// compute maximum residual over the tile
 if (qmaxres) {
// note that I need to remap 2D indices to 1D in order to use the code below
// reuse register
#define step indl
#undef __UPDATE_RES
#define __UPDATE_RES(dx,dy) \
   if (res(IIL(tx,ty)) < res   (IIL(tx+dx,ty+dy)) ) { \
       res(IIL(tx,ty)) = res   (IIL(tx+dx,ty+dy)) ;\
    resind(IIL(tx,ty)) = resind(IIL(tx+dx,ty+dy)) ; \
   }
// NOTE : block size must be a power of two
// otherwise the loop below will not work
  for (step=ntx ; step > 0 ; step>>=1) {
   if (tx < step) {
    __UPDATE_RES(step,0);
   }
   __syncthreads();
  }
  for (step=nty/2 ; step > 0 ; step>>=1) {
   if (ty < step) {
    __UPDATE_RES(0,step);
   }
   __syncthreads();
  }
//
  if (tx==0 & ty==0 ) {
   maxres[ blockIdx.y * gridDim.x + blockIdx.x ]=res(0);
   imax  [ blockIdx.y * gridDim.x + blockIdx.x ]=resind(0);
  } // 0th thread
 } // maxres
} // residual

#undef IOG
#undef IIG
#undef IOL
#undef IIL
#undef __VOLATILE
#undef _DX
#undef __UPDATE_RES
#undef step
#undef residual
#undef tx
#undef ty
#undef bx
#undef by
#undef ntx
#undef nty
#undef p
#undef eps
#undef res
#undef oodxcen
#undef oodxcor
#undef oodycen
#undef oodycor

