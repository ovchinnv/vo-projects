 __global__ void Residual_Cuda_3D(__CUFLOAT *res, 
#ifndef __MGTEX
 __CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, 
#endif
 __CUFLOAT *devoodx, __CUFLOAT *devoody, __CUFLOAT *devoodz,
 const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz) {

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

// this kernel is based on the GS kernel, except that there is no red/black partition
// 2D local sh/mem arrays
 __shared__ float p    [ ( _BSIZE_X + 2 ) * ( _BSIZE_Y + 2 ) ];
 __shared__ float eps  [ ( _BSIZE_X + 2 ) * ( _BSIZE_Y + 2 ) ];
// registers -- 12
// NOTE : tiles are shfted from "front" to "back"
 float pback[1], pfront[1], pcur[1]; // solution z-points ; 2 for red/black
 float eback[1], efront[1], ecur[1]; // epsilon z-points
// 1D local shmem metrics
 __shared__ float oodxcen [ _BSIZE_X + 1 ];
 __shared__ float oodxcor [ _BSIZE_X ];
 __shared__ float oodycen [ _BSIZE_Y + 1 ];
 __shared__ float oodycor [ _BSIZE_Y ];
//registers 12 + 3 = 15
 float oodzcor, oodzcenfront, oodzcenback ;
// thread indices into device memory:
#define tx threadIdx.x
#define ty threadIdx.y
// NOTE that the code was written with the assumption the _SX=2, _SY=1 ; other values may not work yet
// registers 15 + 2=17
 int ix = (blockIdx.x * _BSIZE_X) + tx;
 int iy = (blockIdx.y * _BSIZE_Y) + ty;
//registers 17 + 7 = 22
 float w, e, s, n, b, f, o ;// FD stencil coefficients
//registers 22 + 1 = 23
 int ind, indl ;// indices
// index maps : NOTE that they are zero-based (not one-based as in FORTRAN)
// global memory:
#define IOG(i,j,k)  ( i3b - 1 + (k)*(nx*ny)       + (j)*(nx)   + (i) )
#define IIG(i,j,k)  ( i3  - 1 + (k)*(nx-2)*(ny-2) + (j)*(nx-2) + (i) )
//local memory:
#define IOL(i,j)  ((j)*( _BSIZE_X+2) + (i))
#define IIL(i,j)  ((j)*( _BSIZE_X)   + (i))

 ind=IOG(ix+tx+1,iy+1,0) ;
 pfront[0]=devp(ind)   ; efront[0]=deveps(ind);

 ind=IOG(ix+tx+1,iy+1,1) ;
 pcur[0]=devp(ind) ; ecur[0]=deveps(ind);

// read metrics
//x
 oodxcen[tx] = devoodx[i1 - 1 + ix];
 oodxcor[tx] = devoodx[i1 - 1 + (nx-1) + ix];
//y
 oodycen[ty] = devoody[j1 - 1 + iy];
 oodycor[ty] = devoody[j1 - 1 + (ny-1) + iy];
//
// z - metric
 oodzcenfront=devoodz[k1-1];
 if ( !( (bool)ty || (bool)tx ) ) { // first (0th) thread reads
   oodycen[_BSIZE_Y] = devoody[j1 - 1 + (iy - ty) + _BSIZE_Y];
   oodxcen[_BSIZE_X] = devoodx[i1 - 1 + (ix - tx) + _BSIZE_X];
 }
//
// loop over all z-slices
 for (int k=1 ; k<nz-1 ; k++) {
 __syncthreads(); // This is necessary but I do not know why, since there are no thread-dependend memory access reads between here and the next syncthreads
// populate slice
  ind=IOL(tx+1,ty+1);
  p   [ind] = pcur[0];
  eps [ind] = ecur[0];
// update z-metric (registers) ; might convert to shmem
  oodzcenback=devoodz[k1 - 1 + k];
  oodzcor=devoodz[k1 - 1 + k - 1 + (nz-1)];
//
  if (ty < 2) { // BC in y direction
   p  [IOL(tx+1, (ty%2)*(_BSIZE_Y+1))] = devp  (IOG(ix+1, iy-ty+(ty%2)*(_BSIZE_Y+1), k)); // subtract ty to get first tile coordinate
   eps[IOL(tx+1, (ty%2)*(_BSIZE_Y+1))] = deveps(IOG(ix+1, iy-ty+(ty%2)*(_BSIZE_Y+1), k));
  }
  if (tx < 2) {
   p  [IOL((_BSIZE_X+1)*(tx%2), ty+1)] = devp  (IOG(ix-tx+(_BSIZE_X+1)*(tx%2), iy+1, k));
   eps[IOL((_BSIZE_X+1)*(tx%2), ty+1)] = deveps(IOG(ix-tx+(_BSIZE_X+1)*(tx%2), iy+1, k));
  }

//
// load next p & eps slices
  ind=IOG(ix+1,iy+1,k+1) ;
  pback[0] = devp(ind) ; eback[0] = deveps(ind);
// compute residual
  ind=IIG(ix,iy,k-1) ;
  indl=IOL(tx+1,ty+1) ;
  __syncthreads();
  w = ( eps[indl-1] + ecur[0] ) * ( oodxcor[tx] * oodxcen[tx] );
  e = ( eps[indl+1] + ecur[0] ) * ( oodxcor[tx] * oodxcen[tx+1] );
  s = ( eps[indl-(_BSIZE_X+2)] + ecur[0] ) * ( oodycor[ty] * oodycen[ty] );
  n = ( eps[indl+(_BSIZE_X+2)] + ecur[0] ) * ( oodycor[ty] * oodycen[ty+1] );
  f = ( efront[0] + ecur[0] )              * ( oodzcor * oodzcenfront );
  b = ( eback[0]  + ecur[0] )              * ( oodzcor * oodzcenback );
  o = - 0.5 * ( w + e + s + n + b + f ) + devkappa(ind);
  o = 0.5 / o;
  res[ind] = - o * ( w*p[indl-1] + e*p[indl+1] + s*p[indl-(_BSIZE_X+2)] + n*p[indl+(_BSIZE_X+2)] + f*pfront[0] + b*pback[0]
//                              ) + devrhs[ind]; // RHS o-weighted
                              - 2.0 * devrhs(ind) ); //RHS not weighted
// move forward in z-direction
  pfront[0]=pcur[0];  efront[0]=ecur[0];
  pcur[0]  =pback[0]; ecur[0]  =eback[0];
  oodzcenfront=oodzcenback;
 }
}

#undef IOG
#undef IIG
#undef IOL
#undef IIL
