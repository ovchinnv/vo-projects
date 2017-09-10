
// Red-Black Gauss Seidel CUDA kernel
 __global__ void Gauss_Seidel_Cuda_3D(__CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, __CUFLOAT *devoodx, __CUFLOAT *devoody, __CUFLOAT *devoodz,
                        const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz, const __CFLOAT dt) {

// NOTE : for GS Red/Black read twice as many locations as there are threads, since the update will invilve two passes with half the memory accessed
// 2D local sh/mem arrays
 __shared__ float p    [ ( _SX * _BSIZE_X + 2 ) * ( _SY * _BSIZE_Y + 2 ) ];
 __shared__ float rhs  [ ( _SX * _BSIZE_X     ) * ( _SY * _BSIZE_Y     ) ];
 __shared__ float eps  [ ( _SX * _BSIZE_X + 2 ) * ( _SY * _BSIZE_Y + 2 ) ];
 __shared__ float kappa[ ( _SX * _BSIZE_X     ) * ( _SY * _BSIZE_Y     ) ];
// registers -- 12
// NOTE : tiles are shfted from "front" to "back"
 float pback[2], pfront[2], pcur[2]; // solution z-points ; 2 for red/black
 float eback[2], efront[2], ecur[2]; // epsilon z-points
// 1D local shmem metrics
 __shared__ float oodxcen [ _SX * _BSIZE_X + 1 ];
 __shared__ float oodxcor [ _SX * _BSIZE_X ];
 __shared__ float oodycen [ _SY * _BSIZE_Y + 1 ];
 __shared__ float oodycor [ _SY * _BSIZE_Y ];
//registers 12 + 3 = 15
 float oodzcor, oodzcenfront, oodzcenback ;
// thread indices into device memory:
#define tx threadIdx.x
#define ty threadIdx.y
// NOTE that the code was written with the assumption the _SX=2, _SY=1 ; other values may not work yet
// registers 15 + 2=17
 int ix = _SX*(blockIdx.x * _BSIZE_X) + tx;
 int iy = _SY*(blockIdx.y * _BSIZE_Y) + ty;
//registers 17 + 7 = 22
 float w, e, s, n, b, f, o ;// FD stencil coefficients
//registers 22 + 1 = 23
 int ind ;// index
// index maps : NOTE that they are zero-based (not one-based as in FORTRAN)
// global memory:
#define IOG(i,j,k)  ( i3b - 1 + (k)*(nx*ny)       + (j)*(nx)   + (i) )
#define IIG(i,j,k)  ( i3  - 1 + (k)*(nx-2)*(ny-2) + (j)*(nx-2) + (i) )
//local memory:
#define IOL(i,j)  ((j)*(_SX*_BSIZE_X+2) + (i))
#define IIL(i,j)  ((j)*(_SX*_BSIZE_X)   + (i))

 ind=IOG(ix+tx+1,iy+1,0) ;
 pfront[0]=devp[ind]   ; efront[0]=deveps[ind];       // red

 ind++; ;
 pfront[1]=devp[ind]   ; efront[1]=deveps[ind];       // black

 pcur[0]=devp[IOG(ix+tx+1,iy+1,1)] ; ecur[0]=deveps[IOG(ix+tx+1,iy+1,1)];
 pcur[1]=devp[IOG(ix+tx+2,iy+1,1)] ; ecur[1]=deveps[IOG(ix+tx+2,iy+1,1)];

// read metrics
// note that it may be better to read all metrics into registers
// rather than into shared memory (or copy to registers after reading into shared memory)
// the answer depends on the penalty of accessing the same global data by different threads
// and on the amount of available memories
// x-metrics
 if (ty < _SX) {
// center metric ; NOTE that the center metric is larger, of size 2*_BSIZE_X+1 ; we read two blocks of size _BSIZE_X from each end ; skip one in the middle; then compute it
//  dxcen[ix]         = devdx[i1 - 1 + ix] ; // first half of tile
//  dxcen[ix+_BSIZE+1]= devdx[i1 - 1 + ix +_BSIZE+1] ; // second half of tile
  oodxcen[tx+_BSIZE_X*ty+(ty>0) ] = devoodx[i1 - 1 + ix + _BSIZE_X*ty+(ty>0)];
// corner metric
//  dxcor[ix] = devdx[i1 - 1 + nx-1 + ix ];                   // first  x-half of tile
//  dxcor[ix+_BSIZE_X] = devdx[i1 - 1 + nx-1 + ix +_BSIZE_X]; // second x-half of tile
// in one shot, assuming threIdx.y<2 :
  oodxcor[tx+_BSIZE_X*ty] = devoodx[i1 - 1 + (nx-1) + ix + _BSIZE_X*ty];
 }
// NOTE that if we know that tx=ty we can use a single if clause to read both metrics
// y-metrics
 if (tx < _SY) {
// center
  oodycen[ty+_BSIZE_Y*tx+(tx>0)] = devoody[j1 - 1 + iy + _BSIZE_Y*tx+(tx>0)];
// corner
  oodycor[ty+_BSIZE_Y*tx]     = devoody[j1 - 1 + ny-1 + iy + _BSIZE_Y*tx];
 }
//
 __syncthreads(); // need dxcor/dycor updated
// z - metric
 oodzcenfront=devoodz[0];
// single missing elements in dxcen & dycen (assuming _SX = 2)
// NOTE : it might be better to prefetch z metrics into a shmem array;
 if ( !( (bool)ty || (bool)tx ) ) { // first (0th) thread reads
// deal with y-metric first, since it will usually require global mem read.
  if (_SY<2) { // read from gmem
   oodycen[_BSIZE_Y] = devoody[j1 - 1 + (iy - ty) + _BSIZE_Y];
  } else {
   oodycen[_BSIZE_Y] = 2.0f / ( 1.0f/oodycor[_BSIZE_Y-1]+1.0f/oodycor[_BSIZE_Y-(_SY<2)] ); // ( _SY<2 ) avoids static memory check warnings
  }
//
  if (_SX<2) {  // read from gmem
   oodxcen[_BSIZE_X] = devoodx[i1 - 1 + (ix - tx) + _BSIZE_X];
  } else {
   oodxcen[_BSIZE_X] = 2.0f / ( 1.0f/oodxcor[_BSIZE_X-1]+1.0f/oodxcor[_BSIZE_X-(_SX<2)] );
  }
 }
//
// loop over all z-slices
 for (int k=1 ; k<nz-1 ; k++) {
// might not be optimal due to striped access:
  p   [ IOL(2*tx+1,ty+1) ] = pcur[0];
  p   [ IOL(2*tx+2,ty+1) ] = pcur[1];
  eps [ IOL(2*tx+1,ty+1) ] = ecur[0];
  eps [ IOL(2*tx+2,ty+1) ] = ecur[1];
// update z-metric (registers) ; might convert to shmem
  oodzcenback=devoodz[k];
  oodzcor=devoodz[k-1+nz-1];
// load ghost points
// there is enough work for 4*tx and 2*ty threads for y and x B/C, respectively
// naive approach :
//  if (ty<1) { // x-row of g/p
// near boundary
//   p[IOL(tx+1,         ty)] = devp[IOG(ix+1,         iy,k); // first block
//   p[IOL(tx+1+_BSIZE_X,ty)] = devp[IOG(ix+1+_BSIZE_X,iy,k); // second
// far boundary
//   p[IOL(tx+1,ty+_BSIZE_Y+1)]          = devp[IOG(ix+1,         iy+_BSIZE_Y+1,k); // red
//   p[IOL(tx+1+_BSIZE_X,ty+_BSIZE_Y+1)] = devp[IOG(ix+1+_BSIZE_X,iy+_BSIZE_Y+1,k); // black
//
// epsilon
//   eps[IOL(tx+1,         ty)] = deveps[IOG(ix+1,         iy,k); // first block
//   eps[IOL(tx+1+_BSIZE_X,ty)] = deveps[IOG(ix+1+_BSIZE_X,iy,k); // second
//
//   eps[IOL(tx+1,         ty+_BSIZE_Y+1)] = deveps[IOG(ix+1,         iy+_BSIZE_Y+1,k); // red
//   eps[IOL(tx+1+_BSIZE_X,ty+_BSIZE_Y+1)] = deveps[IOG(ix+1+_BSIZE_X,iy+_BSIZE_Y+1,k); // black
//  }
// same for y-rows -- omitted here
// ...
  if (ty < 4) { // BC in y direction
   p  [IOL(tx+1+_BSIZE_X*(ty%2), (ty/2)*(_BSIZE_Y+1))] = devp  [IOG(ix+1+_BSIZE_X*(ty%2), iy-ty+(ty/2)*(_BSIZE_Y+1), k)]; // subtract ty to get first tile coordinate
   eps[IOL(tx+1+_BSIZE_X*(ty%2), (ty/2)*(_BSIZE_Y+1))] = deveps[IOG(ix+1+_BSIZE_X*(ty%2), iy-ty+(ty/2)*(_BSIZE_Y+1), k)];
  }
  if (tx < 2) {
   p  [IOL((_SX*_BSIZE_X+1)*(tx%2), ty+1)] = devp  [IOG(ix-tx+(_SX*_BSIZE_X+1)*(tx%2), iy+1, k)];
   eps[IOL((_SX*_BSIZE_X+1)*(tx%2), ty+1)] = deveps[IOG(ix-tx+(_SX*_BSIZE_X+1)*(tx%2), iy+1, k)];
  }

//
//  rhs [ IIL(2*tx  ,ty) ] = devrhs [ IIG(ix,  iy,k-1) ];
//  rhs [ IIL(2*tx+1,ty) ] = devrhs [ IIG(ix+1,iy,k-1) ];
  rhs   [ IIL(tx,ty) ]          = devrhs   [ IIG(ix,         iy,k-1) ];
  rhs   [ IIL(tx+_BSIZE_X,ty) ] = devrhs   [ IIG(ix+_BSIZE_X,iy,k-1) ];
//
  kappa [ IIL(tx,         ty) ] = devkappa [ IIG(ix,         iy,k-1) ];
  kappa [ IIL(tx+_BSIZE_X,ty) ] = devkappa [ IIG(ix+_BSIZE_X,iy,k-1) ];
// load top p & eps slices
  ind=IOG(ix+tx+1,iy+1,k+1) ;
  pback[0] = devp[ind] ; eback[0] = deveps[ind];    // red
  ind++;;
  pback[1] = devp[ind] ; eback[1] = deveps[ind];    // black

// compute new value for p
// (1) red
  __syncthreads();
  w = ( eps[IOL(2*tx,   ty+1)] + ecur[0] ) * ( oodxcor[2*tx] * oodxcen[2*tx] );
  e = ( eps[IOL(2*tx+2, ty+1)] + ecur[0] ) * ( oodxcor[2*tx] * oodxcen[2*tx+1] );
  s = ( eps[IOL(2*tx+1, ty)]   + ecur[0] ) * ( oodycor[ty] * oodycen[ty] );
  n = ( eps[IOL(2*tx+1, ty+2)] + ecur[0] ) * ( oodycor[ty] * oodycen[ty+1] );
  f = ( efront[0] + ecur[0] )              * ( oodzcor * oodzcenfront );
  b = ( eback[0]  + ecur[0] )              * ( oodzcor * oodzcenback );
  o = - 0.5 * ( w + e + s + n + b + f ) + kappa[IIL(2*tx,ty)];
  o = 0.5 / o;
// NOTE : RHS is not scaled by o in the CUDA kernels (at least not yet) ; this is a difference from CPU implementation
//  pcur[0] -= dt * ( pcur[0] + o * ( w*p[IOL(2*tx,ty+1)] + e*p[IOL(2*tx+2,ty+1)] + s*p[IOL(2*tx+1,ty)] + n*p[IOL(2*tx+1,ty+2)] + f*pfront[0] + b*pback[0] - 2.0f*rhs[IIL(2*tx,ty)] ) );
  pcur[0] -= dt * ( pcur[0] + o * ( w*p[IOL(2*tx,ty+1)] + e*p[IOL(2*tx+2,ty+1)] + s*p[IOL(2*tx+1,ty)] + n*p[IOL(2*tx+1,ty+2)] + f*pfront[0] + b*pback[0]  ) - rhs[IIL(2*tx,ty)] );
  p[IOL(2*tx+1,ty+1)] = pcur[0]; // upload to shmem
// also upload to device :
  ind=IOG(ix+tx+1,iy+1,k) ;
  devp[ind]=pcur[0];
// (2) black
  __syncthreads();
  w = ( eps[IOL(2*tx+1, ty+1)] + ecur[1] ) * ( oodxcor[2*tx+1] * oodxcen[2*tx+1] );
  e = ( eps[IOL(2*tx+3, ty+1)] + ecur[1] ) * ( oodxcor[2*tx+1] * oodxcen[2*tx+2] );
  s = ( eps[IOL(2*tx+2, ty)]   + ecur[1] ) * ( oodycor[ty]   * oodycen[ty] );
  n = ( eps[IOL(2*tx+2, ty+2)] + ecur[1] ) * ( oodycor[ty]   * oodycen[ty+1] );
  f = ( efront[1] + ecur[1] )              * ( oodzcor * oodzcenfront );
  b = ( eback[1]  + ecur[1] )              * ( oodzcor * oodzcenback );
  o = - 0.5 * ( w + e + s + n + b + f ) + kappa[IIL(2*tx+1,ty)];
  o = 0.5 / o;
//  pcur[1] -= dt * ( pcur[1] + d * ( w*p[IOL(2*tx+1,ty+1)] + e*p[IOL(2*tx+3,ty+1)] + s*p[IOL(2*tx+2,ty)] + n * p[IOL(2*tx+2,ty+2)] + b*pbot[1] + f*ptop[1] - rhs[IIL(2*tx+1,ty)] ))
// upload to device
  ind++;
//  devp[ind]= pcur[1] - dt * ( pcur[1] + o * ( w*p[IOL(2*tx+1,ty+1)] + e*p[IOL(2*tx+3,ty+1)] + s*p[IOL(2*tx+2,ty)] + n*p[IOL(2*tx+2,ty+2)] + f*pfront[1] + b*pback[1]
//                                             - 2.0f*rhs[IIL(2*tx+1,ty)] ) );
  devp[ind]= pcur[1] - dt * ( pcur[1] + o * ( w*p[IOL(2*tx+1,ty+1)] + e*p[IOL(2*tx+3,ty+1)] + s*p[IOL(2*tx+2,ty)] + n*p[IOL(2*tx+2,ty+2)] + f*pfront[1] + b*pback[1])
                                             - rhs[IIL(2*tx+1,ty)] );
// no need to upload to shmem
//  devp[IOG(ix+2,iy+1,k)]=pcur[1];
// move forward in z-direction
  pfront[0]=pcur[0];  efront[0]=ecur[0];
  pcur[0]  =pback[0]; ecur[0]=eback[0];
  pfront[1]=pcur[1];  efront[1]=ecur[1];
  pcur[1]  =pback[1]; ecur[1]=eback[1];
  oodzcenfront=oodzcenback;
 }
}

#undef IOG
#undef IIG
#undef IOL
#undef IIL
