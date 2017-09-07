
// Red-Black Gauss Seidel CUDA kernel
 __global__ void Gauss_Seidel_Cuda_3D(__CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, __CUFLOAT *devdx, __CUFLOAT *devdy, __CUFLOAT *devdz,
                        const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz, const __CFLOAT dt) {

// NOTE : for GS Red/Black read twice as many locations as there are threads, since the update will invilve two passes with half the memory accessed
// 2D local sh/mem arrays
 __shared__ float p    [ ( _SX * _BSIZE_X + 2 ) * ( _SY * _BSIZE_Y + 2 ) ];
 __shared__ float rhs  [ ( _SX * _BSIZE_X     ) * ( _SY * _BSIZE_Y     ) ];
 __shared__ float eps  [ ( _SX * _BSIZE_X + 2 ) * ( _SY * _BSIZE_Y + 2 ) ];
 __shared__ float kappa[ ( _SX * _BSIZE_X     ) * ( _SY * _BSIZE_Y     ) ];
// registers
 float ptop[2], pbot[2], pcur[2]; // solution z-points ; 2 for red/black
 float etop[2], ebot[2], ecur[2]; // epsilon z-points
// 1D local shmem metrics
 __shared__ float dxcen [ _SX * _BSIZE_X + 1 ];
 __shared__ float dxcor [ _SX * _BSIZE_X ];
 __shared__ float dycen [ _SY * _BSIZE_Y + 1 ];
 __shared__ float dycor [ _SY * _BSIZE_Y ];
//registers
 float dzcor, dzcentop, dzcenbot ;
// thread indices into device memory:
#define tx threadIdx.x
#define ty threadIdx.y
// NOTE that the code was written with the assumption the _SX=2, _SY=1 ; other values may not work yet
 int ix = _SX*(blockIdx.x * _BSIZE_X) + tx;
 int iy = _SY*(blockIdx.y * _BSIZE_Y) + ty;
//registers
 float w, e, s, n, b, f, d ;// FD stencil coefficients
//
 int ind1, ind2 ;// indices
// index maps : NOTE that they are zero-based (not one-based as in FORTRAN)
// global memory:
#define IOG(i,j,k)  ( i3b - 1 + (k)*(nx*ny)       + (j)*(nx)   + (i) )
#define IIG(i,j,k)  ( i3  - 1 + (k)*(nx-2)*(ny-2) + (j)*(nx-2) + (i) )
//local memory:
#define IOL(i,j)  ((j)*(_SX*_BSIZE_X+2) + (i))
#define IIL(i,j)  ((j)*(_SX*_BSIZE_X)   + (i))

 ind1=IOG(ix+tx+1,iy+1,0) ;
 ind2=ind1+1 ;

 pbot[0]=devp[ind1]   ; ebot[0]=deveps[ind1];       // red
 pbot[1]=devp[ind2]   ; ebot[1]=deveps[ind2];     // black

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
  dxcen[tx+(_BSIZE_X+1)*ty] = devdx[i1 - 1 + ix + (_BSIZE_X+1)*ty];
// corner metric
//  dxcor[ix] = devdx[i1 - 1 + nx-1 + ix ];                   // first  x-half of tile
//  dxcor[ix+_BSIZE_X] = devdx[i1 - 1 + nx-1 + ix +_BSIZE_X]; // second x-half of tile
// in one shot, assuming threIdx.y<2 :
  dxcor[tx+_BSIZE_X*ty] = devdx[i1 - 1 + nx-1 + ix + _BSIZE_X*ty];
 }
// NOTE that if we know that tx=ty we can use a single if clause to read both metrics
// y-metrics
 if (tx < _SY) {
// center
  dycen[ty+(_BSIZE_Y+1)*tx] = devdy[j1 - 1 + iy + (_BSIZE_Y+1)*tx];
// corner
  dycor[ty+_BSIZE_Y*tx]     = devdy[j1 - 1 + ny-1 + iy + _BSIZE_Y*tx];
 }
//
 __syncthreads(); // need dxcor/dycor updated
// z - metric
 dzcenbot=devdz[1];
// single missing element in dxcen & dycen
// NOTE : it might be better to prefetch z metrics into a shmem array;
 if ( !( (bool)ty || (bool)tx ) ) { // first (0th) thread reads
  for (int i=1 ; i < _SX ; i++) { // should work for any _SX ; could be slow ; 
   dxcen[i*_BSIZE_X] = 0.5 * ( dxcor[i*_BSIZE_X]+dxcor[i*_BSIZE_X+1] );
  }
  for (int j=1 ; j < _SY ; j++) {
   dycen[j*_BSIZE_Y] = 0.5 * ( dycor[j*_BSIZE_Y]+dycor[j*_BSIZE_Y+1]);
  }
 }

// loop over all z-slices
 for (int k=1 ; k<nz-1 ; k++) {
// might not be optimal due to striped access:
  p   [ IOL(2*tx+1,ty+1) ] = pcur[0];
  p   [ IOL(2*tx+2,ty+1) ] = pcur[1];
  eps [ IOL(2*tx+1,ty+1) ] = ecur[0];
  eps [ IOL(2*tx+2,ty+1) ] = ecur[1];
// update z-metric (registers) ; might convert to shmem
  dzcentop=devdz[k];
  dzcor=devdz[k-1+nz-1];
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
  ind1=IOG(ix+tx+1,iy+1,k+1) ;
  ind2=ind1+1 ;
  ptop[0] = devp[ind1] ; ebot[0] = deveps[ind1];    // red
  ptop[1] = devp[ind2] ; ebot[1] = deveps[ind2];    // black

// compute new value for p
// (1) red
  __syncthreads();
  w = ( eps[IOL(2*tx,   ty+1)] + ecur[0] ) / ( dxcor[2*tx] * dxcen[2*tx] );
  e = ( eps[IOL(2*tx+2, ty+1)] + ecur[0] ) / ( dxcor[2*tx] * dxcen[2*tx+1] );
  s = ( eps[IOL(2*tx+1, ty)]   + ecur[0] ) / ( dycor[ty] * dycen[ty] );
  n = ( eps[IOL(2*tx+1, ty+2)] + ecur[0] ) / ( dycor[ty] * dycen[ty+1] );
  b = ( ebot[0] + ecur[0] )              / ( dzcor * dzcenbot );
  f = ( etop[0] + ecur[0] )              / ( dzcor * dzcentop );
  d = - 0.5 * ( w + e + s + n + b + f ) + kappa[IIL(2*tx,ty)];
  d = 0.5 / d;
  pcur[0] -= dt * ( pcur[0] + d * ( w*p[IOL(2*tx,ty+1)] + e*p[IOL(2*tx+2,ty+1)] + s*p[IOL(2*tx+1,ty)] + n * p[IOL(2*tx+1,ty+2)] + b*pbot[0] + f*ptop[0] - rhs[IIL(2*tx,ty)] ) );
  p[IOL(2*tx+1,ty+1)] = pcur[0]; // upload to shmem
  __syncthreads();
  w = ( eps[IOL(2*tx+1, ty+1)] + ecur[1] ) / ( dxcor[2*tx+1] * dxcen[2*tx+1] );
  e = ( eps[IOL(2*tx+3, ty+1)] + ecur[1] ) / ( dxcor[2*tx+1] * dxcen[2*tx+2] );
  s = ( eps[IOL(2*tx+2, ty)]   + ecur[1] ) / ( dycor[ty]   * dycen[ty] );
  n = ( eps[IOL(2*tx+2, ty+2)] + ecur[1] ) / ( dycor[ty]   * dycen[ty+1] );
  b = ( ebot[1] + ecur[1] )              / ( dzcor * dzcenbot );
  f = ( etop[1] + ecur[1] )              / ( dzcor * dzcentop );
  d = - 0.5 * ( w + e + s + n + b + f ) + kappa[IIL(2*tx+1,ty)];
  d = 0.5 / d;
//  pcur[1] -= dt * ( pcur[1] + d * ( w*p[IOL(2*tx+1,ty+1)] + e*p[IOL(2*tx+3,ty+1)] + s*p[IOL(2*tx+2,ty)] + n * p[IOL(2*tx+2,ty+2)] + b*pbot[1] + f*ptop[1] - rhs[IIL(2*tx+1,ty)] ))
// directly to device :
  ind1=IOG(ix+tx+1,iy+1,k) ;
  ind2=ind1+1 ;
  devp[ind2]= pcur[1] - dt * ( pcur[1] + d * ( w*p[IOL(2*tx+1,ty+1)] + e*p[IOL(2*tx+3,ty+1)] + s*p[IOL(2*tx+2,ty)] + n*p[IOL(2*tx+2,ty+2)] + b*pbot[1] + f*ptop[1] - rhs[IIL(2*tx+1,ty)] ) );
// no need to upload to shmem
// upload to device
  devp[ind1]=pcur[0];
//  devp[IOG(ix+2,iy+1,k)]=pcur[1];
// move forward in z-direction
  pbot[0]=pcur[0]; ebot[0]=ecur[0];
  pcur[0]=ptop[0]; ecur[0]=etop[0];
  dzcenbot=dzcentop;
 }
}


