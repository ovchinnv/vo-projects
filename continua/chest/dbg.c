 __global__ void Gauss_Seidel_Cuda_3D(__CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, __CUFLOAT *devoodx, __CUFLOAT *devoody, __CUFLOAT *devoodz,
                        const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz, const __CFLOAT dt,
                        const int which, const bool qpinitzero) {
 volatile __shared__ float p [ (( _SX * _BGSMAX_X + 2 ) * ( _SY * _BGSMAX_Y + 2 )) ];
 __shared__ float rhs [ (( _SX * _BGSMAX_X ) * ( _SY * _BGSMAX_Y )) ];
 __shared__ float eps [ (( _SX * _BGSMAX_X + 2 ) * ( _SY * _BGSMAX_Y + 2 )) ];
 __shared__ float kappa[ (( _SX * _BGSMAX_X ) * ( _SY * _BGSMAX_Y )) ];
 __shared__ float oodxcen [ (_SX * _BGSMAX_X + 1) ];
 __shared__ float oodxcor [ (_SX * _BGSMAX_X) ];
 __shared__ float oodycen [ (_SY * _BGSMAX_Y + 1) ];
 __shared__ float oodycor [ (_SY * _BGSMAX_Y) ];
 float pback[2], pfront[2], pcur[2];
 float eback[2], efront[2], ecur[2];
 float oodzcor, oodzcenfront, oodzcenback ;
 unsigned int ix = _SX*(blockIdx.x * _BGSMAX_X) + threadIdx.x;
 unsigned int iy = _SY*(blockIdx.y * _BGSMAX_Y) + threadIdx.y;
 float w, e, s, n, b, f, o ;
 unsigned int ind, indl ;
 ind=( i3b - 1 + (0)*(nx*ny) + (iy+1)*(nx) + (ix+threadIdx.x+1) ) ;
 pfront[0]=devp[ind];
 efront[0]=deveps[ind];
 ind++;
 pfront[1]=devp[ind];
 efront[1]=deveps[ind];
 ind=( i3b - 1 + (1)*(nx*ny) + (iy+1)*(nx) + (ix+threadIdx.x+1) ) ;
 pcur[0]=devp[ind] ; ecur[0]=deveps[ind];
 ind++;
 pcur[1]=devp[ind] ; ecur[1]=deveps[ind];
 if (threadIdx.y < _SX) {
  oodxcen[threadIdx.x+_BGSMAX_X*threadIdx.y+(threadIdx.y>0)] = devoodx[i1 - 1 + ix + _BGSMAX_X*threadIdx.y+(threadIdx.y>0)];
  oodxcor[threadIdx.x+_BGSMAX_X*threadIdx.y] = devoodx[i1 - 1 + (nx-1) + ix + _BGSMAX_X*threadIdx.y];
 }
 if (threadIdx.x < _SY) {
  oodycen[threadIdx.y+_BGSMAX_Y*threadIdx.x+(threadIdx.x>0)] = devoody[j1 - 1 + iy + _BGSMAX_Y*threadIdx.x+(threadIdx.x>0)];
  oodycor[threadIdx.y+_BGSMAX_Y*threadIdx.x] = devoody[j1 - 1 + ny-1 + iy + _BGSMAX_Y*threadIdx.x];
 }
 oodzcenfront=devoodz[k1 - 1];
 __syncthreads();
 if ( !( (bool)threadIdx.y || (bool)threadIdx.x ) ) {
  if (_SY<2) {
   oodycen[_BGSMAX_Y] = devoody[j1 - 1 + (iy - threadIdx.y) + _BGSMAX_Y];
  } else {
   oodycen[_BGSMAX_Y] = 2.0 / ( 1.0/oodycor[_BGSMAX_Y-1]+1.0/oodycor[_BGSMAX_Y-(_SY<2)] );
  }
   oodxcen[_BGSMAX_X] = devoodx[i1 - 1 + (ix - threadIdx.x) + _BGSMAX_X];
 }
 for (unsigned int k=1 ; k<nz-1 ; k++) {
  __syncthreads();
  ind=((threadIdx.y+1)*(_SX*_BGSMAX_X+2) + (2*threadIdx.x+1));
  p[ind] = pcur[0];
  eps[ind] = ecur[0];
  ind++;
  p[ind] = pcur[1];
  eps[ind] = ecur[1];
  oodzcenback=devoodz[k1 - 1 + k];
  oodzcor=devoodz[k1 - 1 + k - 1 + (nz-1) ];
  if (threadIdx.y < 4) {
   p[(((threadIdx.y/2)*(_BGSMAX_Y+1))*(_SX*_BGSMAX_X+2) + (threadIdx.x+1+_BGSMAX_X*(threadIdx.y%2)))] = devp [( i3b - 1 + (k)*(nx*ny) + (iy-threadIdx.y+(threadIdx.y/2)*(_BGSMAX_Y+1))*(nx) + (ix+1+_BGSMAX_X*(threadIdx.y%2)) )];
   eps[(((threadIdx.y/2)*(_BGSMAX_Y+1))*(_SX*_BGSMAX_X+2) + (threadIdx.x+1+_BGSMAX_X*(threadIdx.y%2)))] = deveps[( i3b - 1 + (k)*(nx*ny) + (iy-threadIdx.y+(threadIdx.y/2)*(_BGSMAX_Y+1))*(nx) + (ix+1+_BGSMAX_X*(threadIdx.y%2)) )];
  }
  if (threadIdx.x < 2) {
   p[((threadIdx.y+1)*(_SX*_BGSMAX_X+2) + ((_SX*_BGSMAX_X+1)*(threadIdx.x%2)))] = devp [( i3b - 1 + (k)*(nx*ny) + (iy+1)*(nx) + (ix-threadIdx.x+(_SX*_BGSMAX_X+1)*(threadIdx.x%2)) )];
   eps[((threadIdx.y+1)*(_SX*_BGSMAX_X+2) + ((_SX*_BGSMAX_X+1)*(threadIdx.x%2)))] = deveps[( i3b - 1 + (k)*(nx*ny) + (iy+1)*(nx) + (ix-threadIdx.x+(_SX*_BGSMAX_X+1)*(threadIdx.x%2)) )];
  }
  rhs[((threadIdx.y)*(_SX*_BGSMAX_X) + (threadIdx.x))] = devrhs[( i3 - 1 + (k-1)*(nx-2)*(ny-2) + (iy)*(nx-2) + (ix) )];
  rhs[((threadIdx.y)*(_SX*_BGSMAX_X) + (threadIdx.x+_BGSMAX_X))] = devrhs[( i3 - 1 + (k-1)*(nx-2)*(ny-2) + (iy)*(nx-2) + (ix+_BGSMAX_X) )];
  kappa[((threadIdx.y)*(_SX*_BGSMAX_X) + (threadIdx.x))] = devkappa[( i3 - 1 + (k-1)*(nx-2)*(ny-2) + (iy)*(nx-2) + (ix) )];
  kappa[((threadIdx.y)*(_SX*_BGSMAX_X) + (threadIdx.x+_BGSMAX_X))] = devkappa[( i3 - 1 + (k-1)*(nx-2)*(ny-2) + (iy)*(nx-2) + (ix+_BGSMAX_X) )];
  ind=( i3b - 1 + (k+1)*(nx*ny) + (iy+1)*(nx) + (ix+threadIdx.x+1) ) ;
  pback[0] = devp[ind] ; eback[0] = deveps[ind];
  ind++;
  pback[1] = devp[ind] ; eback[1] = deveps[ind];
  ind=( i3b - 1 + (k)*(nx*ny) + (iy+1)*(nx) + (ix+threadIdx.x+1+(threadIdx.y%2)) ) ;
  indl=((threadIdx.y+1)*(_SX*_BGSMAX_X+2) + (2*threadIdx.x+1+(threadIdx.y%2))) ;
  __syncthreads();
  w = ( eps[indl-1] + ecur[(threadIdx.y%2)] ) * ( oodxcor[2*threadIdx.x+(threadIdx.y%2)] * oodxcen[2*threadIdx.x+(threadIdx.y%2)] );
  e = ( eps[indl+1] + ecur[(threadIdx.y%2)] ) * ( oodxcor[2*threadIdx.x+(threadIdx.y%2)] * oodxcen[2*threadIdx.x+1+(threadIdx.y%2)] );
  s = ( eps[indl-(_SX*_BGSMAX_X+2)] + ecur[(threadIdx.y%2)] ) * ( oodycor[threadIdx.y] * oodycen[threadIdx.y] );
  n = ( eps[indl+(_SX*_BGSMAX_X+2)] + ecur[(threadIdx.y%2)] ) * ( oodycor[threadIdx.y] * oodycen[threadIdx.y+1] );
  f = ( efront[(threadIdx.y%2)] + ecur[(threadIdx.y%2)] ) * ( oodzcor * oodzcenfront );
  b = ( eback[(threadIdx.y%2)] + ecur[(threadIdx.y%2)] ) * ( oodzcor * oodzcenback );
  o = - 0.5 * ( w + e + s + n + b + f ) + kappa[((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x+(threadIdx.y%2)))];
  o = 0.5 / o;
  pcur[(threadIdx.y%2)] -= dt * ( pcur[(threadIdx.y%2)] + o * ( w*p[indl-1] + e*p[indl+1] + s*p[indl-(_SX*_BGSMAX_X+2)] + n*p[indl+(_SX*_BGSMAX_X+2)] + f*pfront[(threadIdx.y%2)] + b*pback[(threadIdx.y%2)]
                              - 2.0 * rhs[((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x+(threadIdx.y%2)))] ) );
  devp[ind]=pcur[(threadIdx.y%2)];
   p[indl] = pcur[(threadIdx.y%2)];
  ind=( i3b - 1 + (k)*(nx*ny) + (iy+1)*(nx) + (ix+threadIdx.x+1+(1-threadIdx.y%2)) ) ;
  indl=((threadIdx.y+1)*(_SX*_BGSMAX_X+2) + (2*threadIdx.x+1+(1-threadIdx.y%2))) ;
  __syncthreads();
  w = ( eps[indl-1] + ecur[(1-threadIdx.y%2)] ) * ( oodxcor[2*threadIdx.x+(1-threadIdx.y%2)] * oodxcen[2*threadIdx.x+(1-threadIdx.y%2)] );
  e = ( eps[indl+1] + ecur[(1-threadIdx.y%2)] ) * ( oodxcor[2*threadIdx.x+(1-threadIdx.y%2)] * oodxcen[2*threadIdx.x+1+(1-threadIdx.y%2)] );
  s = ( eps[indl-(_SX*_BGSMAX_X+2)] + ecur[(1-threadIdx.y%2)] ) * ( oodycor[threadIdx.y] * oodycen[threadIdx.y] );
  n = ( eps[indl+(_SX*_BGSMAX_X+2)] + ecur[(1-threadIdx.y%2)] ) * ( oodycor[threadIdx.y] * oodycen[threadIdx.y+1] );
  f = ( efront[(1-threadIdx.y%2)] + ecur[(1-threadIdx.y%2)] ) * ( oodzcor * oodzcenfront );
  b = ( eback[(1-threadIdx.y%2)] + ecur[(1-threadIdx.y%2)] ) * ( oodzcor * oodzcenback );
  o = - 0.5 * ( w + e + s + n + b + f ) + kappa[((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x+(1-threadIdx.y%2)))];
  o = 0.5 / o;
  pcur[(1-threadIdx.y%2)] -= dt * ( pcur[(1-threadIdx.y%2)] + o * ( w*p[indl-1] + e*p[indl+1] + s*p[indl-(_SX*_BGSMAX_X+2)] + n*p[indl+(_SX*_BGSMAX_X+2)] + f*pfront[(1-threadIdx.y%2)] + b*pback[(1-threadIdx.y%2)]
                              - 2.0 * rhs[((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x+(1-threadIdx.y%2)))] ) );
  devp[ind]=pcur[(1-threadIdx.y%2)];
  pfront[0]=pcur[0]; efront[0]=ecur[0];
  pcur[0] =pback[0]; ecur[0] =eback[0];
  pfront[1]=pcur[1]; efront[1]=ecur[1];
  pcur[1] =pback[1]; ecur[1] =eback[1];
  oodzcenfront=oodzcenback;
 }
}
 __global__ void Apply_BC_Cuda_3D(__CUFLOAT *devu, __CUFLOAT *devg, const __CINT i3b, const __CINT ind2, const __CINT nx, const __CINT ny, const __CINT nz,
                                   const __CINT boundary, const __CINT bctype, const __CUFLOAT wgt) {
      int ix = blockIdx.x*_BCDIM_X+threadIdx.x ;
      int iy = blockIdx.y*_BCDIM_Y+threadIdx.y ;
      int iz = blockIdx.z*_BCDIM_Z+threadIdx.z ;
      int ib, ie, jb, je, kb, ke ;
      int di, dj, dk ;
      int dir ;
      ib=1; ie=nx-2;
      jb=1; je=ny-2;
      kb=1; ke=nz-2;
      if (boundary==1) {
         ib--; ie=ib; di=+1; dj=0; dk=0; dir=di;
        } else if (boundary==2) {
         ie++; ib=ie; di=-1; dj=0; dk=0; dir=di;
        } else if (boundary==3) {
         jb--; je=jb; dj=+1; di=0; dk=0; dir=dj;
        } else if (boundary==4) {
         je++; jb=je; dj=-1; di=0; dk=0; dir=dj;
        } else if (boundary==5) {
         kb--; ke=kb; dk=+1; di=0; dj=0; dir=dk;
        } else if (boundary==6) {
         ke++; kb=ke; dk=-1; di=0; dj=0; dir=dk;
        } else {
         ke=0;
      }
      if (bctype==1 || bctype==4) {
         devu[( i3b - 1 + (iz+kb)*(nx*ny) + (iy+jb)*(nx) + (ix+ib) )] = devu[( i3b - 1 + (iz+kb+dk)*(nx*ny) + (iy+jb+dj)*(nx) + (ix+ib+di) )] + wgt * (devg[( ind2 - 1 + (iz)*((ie-ib+1)*(je-jb+1)) + (iy)*(ie-ib+1) + (ix) )] - devu[( i3b - 1 + (iz+kb+dk)*(nx*ny) + (iy+jb+dj)*(nx) + (ix+ib+di) )]);
        } else if (bctype==2 || bctype==5) {
         devu[( i3b - 1 + (iz+kb)*(nx*ny) + (iy+jb)*(nx) + (ix+ib) )] = devu[( i3b - 1 + (iz+kb+dk)*(nx*ny) + (iy+jb+dj)*(nx) + (ix+ib+di) )] - dir * wgt * devg[( ind2 - 1 + (iz)*((ie-ib+1)*(je-jb+1)) + (iy)*(ie-ib+1) + (ix) )];
        } else if (bctype==3 ) {
         dk*=(nz-2); dj*=(ny-2); di*=(nx-2);
         devu[( i3b - 1 + (iz+kb)*(nx*ny) + (iy+jb)*(nx) + (ix+ib) )] = devu[( i3b - 1 + (iz+kb+dk)*(nx*ny) + (iy+jb+dj)*(nx) + (ix+ib+di) )] ;
      }
 }
 __global__ void Apply_BC_Cuda_3Dx2(__CUFLOAT *devu, __CUFLOAT *devg, __CUFLOAT *devg2, const __CINT i3b, const __CINT ind2, const __CINT nx, const __CINT ny, const __CINT nz,
                                   const __CINT boundary, const __CINT boundary2, const __CINT bctype, const __CINT bctype2, const __CFLOAT wgt, const __CUFLOAT wgt2, const bool qpinitzero){
      int ix = blockIdx.x*_BCDIM_X+threadIdx.x ;
      int iy = blockIdx.y*_BCDIM_Y+threadIdx.y ;
      int iz = blockIdx.z*_BCDIM_Z+threadIdx.z ;
      int ib, ie, jb, je, kb, ke ;
      int di, dj, dk ;
      int dir ;
      ib=1; ie=nx-2;
      jb=1; je=ny-2;
      kb=1; ke=nz-2;
      if (boundary==1) {
         ib--; ie=ib; di=+1; dj=0; dk=0; dir=di;
        } else if (boundary==2) {
         ie++; ib=ie; di=-1; dj=0; dk=0; dir=di;
        } else if (boundary==3) {
         jb--; je=jb; dj=+1; di=0; dk=0; dir=dj;
        } else if (boundary==4) {
         je++; jb=je; dj=-1; di=0; dk=0; dir=dj;
        } else if (boundary==5) {
         kb--; ke=kb; dk=+1; di=0; dj=0; dir=dk;
        } else if (boundary==6) {
         ke++; kb=ke; dk=-1; di=0; dj=0; dir=dk;
      }
      if (bctype==1 || bctype==4) {
         devu[( i3b - 1 + (iz+kb)*(nx*ny) + (iy+jb)*(nx) + (ix+ib) )] = (1.0 - wgt2) * ((qpinitzero) ? 0.0 : devu[( i3b - 1 + (iz+kb+dk)*(nx*ny) + (iy+jb+dj)*(nx) + (ix+ib+di) )]) + wgt2 * devg[( ind2 - 1 + (iz)*((ie-ib+1)*(je-jb+1)) + (iy)*(ie-ib+1) + (ix) )];
        } else if (bctype==2 || bctype==5) {
         devu[( i3b - 1 + (iz+kb)*(nx*ny) + (iy+jb)*(nx) + (ix+ib) )] = ((qpinitzero) ? 0.0 : devu[( i3b - 1 + (iz+kb+dk)*(nx*ny) + (iy+jb+dj)*(nx) + (ix+ib+di) )]) - dir * wgt * devg[( ind2 - 1 + (iz)*((ie-ib+1)*(je-jb+1)) + (iy)*(ie-ib+1) + (ix) )];
        } else if (bctype==3 ) {
         dk*=(nz-2); dj*=(ny-2); di*=(nx-2);
         devu[( i3b - 1 + (iz+kb)*(nx*ny) + (iy+jb)*(nx) + (ix+ib) )] = ((qpinitzero) ? 0.0 : devu[( i3b - 1 + (iz+kb+dk)*(nx*ny) + (iy+jb+dj)*(nx) + (ix+ib+di) )]);
      }
      ib=1; ie=nx-2;
      jb=1; je=ny-2;
      kb=1; ke=nz-2;
      if (boundary2==1) {
         ib--; ie=ib; di=+1; dj=0; dk=0; dir=di;
        } else if (boundary2==2) {
         ie++; ib=ie; di=-1; dj=0; dk=0; dir=di;
        } else if (boundary2==3) {
         jb--; je=jb; dj=+1; di=0; dk=0; dir=dj;
        } else if (boundary2==4) {
         je++; jb=je; dj=-1; di=0; dk=0; dir=dj;
        } else if (boundary2==5) {
         kb--; ke=kb; dk=+1; di=0; dj=0; dir=dk;
        } else if (boundary2==6) {
         ke++; kb=ke; dk=-1; di=0; dj=0; dir=dk;
        } else {
         devg[( ind2 - 1 + (iz)*((ie-ib+1)*(je-jb+1)) + (iy)*(ie-ib+1) + (ix) )]=0.f;
      }
      if (bctype2==1 || bctype2==4) {
         devu[( i3b - 1 + (iz+kb)*(nx*ny) + (iy+jb)*(nx) + (ix+ib) )] = (1.0 - wgt2) * ((qpinitzero) ? 0.0 : devu[( i3b - 1 + (iz+kb+dk)*(nx*ny) + (iy+jb+dj)*(nx) + (ix+ib+di) )]) + wgt2 * devg2[( ind2 - 1 + (iz)*((ie-ib+1)*(je-jb+1)) + (iy)*(ie-ib+1) + (ix) )];
        } else if (bctype2==2 || bctype2==5) {
         devu[( i3b - 1 + (iz+kb)*(nx*ny) + (iy+jb)*(nx) + (ix+ib) )] = ((qpinitzero) ? 0.0 : devu[( i3b - 1 + (iz+kb+dk)*(nx*ny) + (iy+jb+dj)*(nx) + (ix+ib+di) )]) - dir * wgt2 * devg2[( ind2 - 1 + (iz)*((ie-ib+1)*(je-jb+1)) + (iy)*(ie-ib+1) + (ix) )];
        } else if (bctype2==3 ) {
         dk*=(nz-2); dj*=(ny-2); di*=(nx-2);
         devu[( i3b - 1 + (iz+kb)*(nx*ny) + (iy+jb)*(nx) + (ix+ib) )] = ((qpinitzero) ? 0.0 : devu[( i3b - 1 + (iz+kb+dk)*(nx*ny) + (iy+jb+dj)*(nx) + (ix+ib+di) )]);
      }
 }
 __global__ void Residual_Cuda_3D(__CUFLOAT *devres,
 __CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa,
 __CUFLOAT *devoodx, __CUFLOAT *devoody, __CUFLOAT *devoodz,
 const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz, const __CINT qmaxres, const __CINT qresnorm,
 __CUFLOAT *maxres, __CINT *imax) {
 __shared__ float p [ (( _SX * _BGSMAX_X + 2 ) * ( _SY * _BGSMAX_Y + 2 )) ];
 __shared__ float eps [ (( _SX * _BGSMAX_X + 2 ) * ( _SY * _BGSMAX_Y + 2 )) ];
 __shared__ float res [ (( _SX * _BGSMAX_X ) * ( _SY * _BGSMAX_Y )) ];
 __shared__ __CINT resind [ (( _SX * _BGSMAX_X ) * ( _SY * _BGSMAX_Y )) ];
 __shared__ float oodxcen [ (_SX * _BGSMAX_X + 1) ];
 __shared__ float oodxcor [ (_SX * _BGSMAX_X) ];
 __shared__ float oodycen [ (_SY * _BGSMAX_Y + 1) ];
 __shared__ float oodycor [ (_SY * _BGSMAX_Y) ];
 float pback[2], pfront[2], pcur[2];
 float eback[2], efront[2], ecur[2];
 float rhs[2], kappa[2] ;
 float oodzcor, oodzcenfront, oodzcenback ;
 unsigned int ix = _SX*(blockIdx.x * _BGSMAX_X) + threadIdx.x;
 unsigned int iy = _SY*(blockIdx.y * _BGSMAX_Y) + threadIdx.y;
 float w, e, s, n, b, f, o ;
 unsigned int ind, indl ;
 ind=( i3b - 1 + (0)*(nx*ny) + (iy+1)*(nx) + (ix+threadIdx.x+1) ) ;
 pfront[0]=devp[ind] ; efront[0]=deveps[ind];
 ind++;
 pfront[1]=devp[ind] ; efront[1]=deveps[ind];
 ind=( i3b - 1 + (1)*(nx*ny) + (iy+1)*(nx) + (ix+threadIdx.x+1) ) ;
 pcur[0]=devp[ind] ; ecur[0]=deveps[ind];
 ind++;
 pcur[1]=devp[ind] ; ecur[1]=deveps[ind];
 if (threadIdx.y < _SX) {
  oodxcen[threadIdx.x+_BGSMAX_X*threadIdx.y+(threadIdx.y>0)] = devoodx[i1 - 1 + ix + _BGSMAX_X*threadIdx.y+(threadIdx.y>0)];
  oodxcor[threadIdx.x+_BGSMAX_X*threadIdx.y] = devoodx[i1 - 1 + (nx-1) + ix + _BGSMAX_X*threadIdx.y];
 }
 if (threadIdx.x < _SY) {
  oodycen[threadIdx.y+_BGSMAX_Y*threadIdx.x+(threadIdx.x>0)] = devoody[j1 - 1 + iy + _BGSMAX_Y*threadIdx.x+(threadIdx.x>0)];
  oodycor[threadIdx.y+_BGSMAX_Y*threadIdx.x] = devoody[j1 - 1 + ny-1 + iy + _BGSMAX_Y*threadIdx.x];
 }
 oodzcenfront=devoodz[k1 - 1];
 __syncthreads();
 if ( !( (bool)threadIdx.y || (bool)threadIdx.x ) ) {
  if (_SY<2) {
   oodycen[_BGSMAX_Y] = devoody[j1 - 1 + (iy - threadIdx.y) + _BGSMAX_Y];
  } else {
   oodycen[_BGSMAX_Y] = 2.0 / ( 1.0/oodycor[_BGSMAX_Y-1]+1.0/oodycor[_BGSMAX_Y-(_SY<2)] );
  }
   oodxcen[_BGSMAX_X] = devoodx[i1 - 1 + (ix - threadIdx.x) + _BGSMAX_X];
 }
 if (qmaxres) {
  indl=((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x));
  res[indl] =-INFINITY; resind[indl]=-1;
  res[++indl]=-INFINITY; resind[indl]=-1;
 }
 for (unsigned int k=1 ; k<nz-1 ; k++) {
  __syncthreads();
  ind=((threadIdx.y+1)*(_SX*_BGSMAX_X+2) + (2*threadIdx.x+1));
  p[ind] = pcur[0];
  eps[ind] = ecur[0];
  ind++;
  p[ind] = pcur[1];
  eps[ind] = ecur[1];
  oodzcenback=devoodz[k1 - 1 + k];
  oodzcor=devoodz[k1 - 1 + k - 1 + (nz-1) ];
  if (threadIdx.y < 4) {
   p[(((threadIdx.y/2)*(_BGSMAX_Y+1))*(_SX*_BGSMAX_X+2) + (threadIdx.x+1+_BGSMAX_X*(threadIdx.y%2)))] = devp[( i3b - 1 + (k)*(nx*ny) + (iy-threadIdx.y+(threadIdx.y/2)*(_BGSMAX_Y+1))*(nx) + (ix+1+_BGSMAX_X*(threadIdx.y%2)) )];
   eps[(((threadIdx.y/2)*(_BGSMAX_Y+1))*(_SX*_BGSMAX_X+2) + (threadIdx.x+1+_BGSMAX_X*(threadIdx.y%2)))] = deveps[( i3b - 1 + (k)*(nx*ny) + (iy-threadIdx.y+(threadIdx.y/2)*(_BGSMAX_Y+1))*(nx) + (ix+1+_BGSMAX_X*(threadIdx.y%2)) )];
  }
  if (threadIdx.x < 2) {
   p[((threadIdx.y+1)*(_SX*_BGSMAX_X+2) + ((_SX*_BGSMAX_X+1)*(threadIdx.x%2)))] = devp[( i3b - 1 + (k)*(nx*ny) + (iy+1)*(nx) + (ix-threadIdx.x+(_SX*_BGSMAX_X+1)*(threadIdx.x%2)) )];
   eps[((threadIdx.y+1)*(_SX*_BGSMAX_X+2) + ((_SX*_BGSMAX_X+1)*(threadIdx.x%2)))] = deveps[( i3b - 1 + (k)*(nx*ny) + (iy+1)*(nx) + (ix-threadIdx.x+(_SX*_BGSMAX_X+1)*(threadIdx.x%2)) )];
  }
  ind=( i3b - 1 + (k+1)*(nx*ny) + (iy+1)*(nx) + (ix+threadIdx.x+1) ) ;
  pback[0] = devp[ind] ; eback[0] = deveps[ind];
  ind++;
  pback[1] = devp[ind] ; eback[1] = deveps[ind];
  ind=( i3 - 1 + (k-1)*(nx-2)*(ny-2) + (iy)*(nx-2) + (ix+threadIdx.x+0) ) ;
  indl=((threadIdx.y+1)*(_SX*_BGSMAX_X+2) + (2*threadIdx.x+1+0)) ;
  __syncthreads();
  rhs[0] =devrhs[ind];
  kappa[0] =devkappa[ind];
  rhs[1-0] =devrhs[ind+1-0];
  kappa[1-0]=devkappa[ind+1-0];
  w = ( eps[indl-1] + ecur[0] ) * ( oodxcor[2*threadIdx.x+0] * oodxcen[2*threadIdx.x+0] );
  e = ( eps[indl+1] + ecur[0] ) * ( oodxcor[2*threadIdx.x+0] * oodxcen[2*threadIdx.x+1+0] );
  s = ( eps[indl-(_SX*_BGSMAX_X+2)] + ecur[0] ) * ( oodycor[threadIdx.y] * oodycen[threadIdx.y] );
  n = ( eps[indl+(_SX*_BGSMAX_X+2)] + ecur[0] ) * ( oodycor[threadIdx.y] * oodycen[threadIdx.y+1] );
  f = ( efront[0] + ecur[0] ) * ( oodzcor * oodzcenfront );
  b = ( eback[0] + ecur[0] ) * ( oodzcor * oodzcenback );
  o = - ( w + e + s + n + b + f ) + 2.0 * kappa[0];
  e = rhs[0] - 0.5 * ( o * pcur[0] + w*p[indl-1] + e*p[indl+1] + s*p[indl-(_SX*_BGSMAX_X+2)] + n*p[indl+(_SX*_BGSMAX_X+2)] + f*pfront[0] + b*pback[0] );
  devres[ind-i3+1] = e;
  if (qmaxres) { if (qresnorm) { e=fabs(2.0f*e/(o)); } else { e=fabs(e); } if (res[((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x+0))] < e) { res[((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x+0))] = e ; resind[((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x+0))] = ind - i3 + 1 ; } }
  ind=( i3 - 1 + (k-1)*(nx-2)*(ny-2) + (iy)*(nx-2) + (ix+threadIdx.x+1) ) ;
  indl=((threadIdx.y+1)*(_SX*_BGSMAX_X+2) + (2*threadIdx.x+1+1)) ;
  w = ( eps[indl-1] + ecur[1] ) * ( oodxcor[2*threadIdx.x+1] * oodxcen[2*threadIdx.x+1] );
  e = ( eps[indl+1] + ecur[1] ) * ( oodxcor[2*threadIdx.x+1] * oodxcen[2*threadIdx.x+1+1] );
  s = ( eps[indl-(_SX*_BGSMAX_X+2)] + ecur[1] ) * ( oodycor[threadIdx.y] * oodycen[threadIdx.y] );
  n = ( eps[indl+(_SX*_BGSMAX_X+2)] + ecur[1] ) * ( oodycor[threadIdx.y] * oodycen[threadIdx.y+1] );
  f = ( efront[1] + ecur[1] ) * ( oodzcor * oodzcenfront );
  b = ( eback[1] + ecur[1] ) * ( oodzcor * oodzcenback );
  o = - ( w + e + s + n + b + f ) + 2.0 * kappa[1];
  e = rhs[1] - 0.5 * ( o * pcur[1] + w*p[indl-1] + e*p[indl+1] + s*p[indl-(_SX*_BGSMAX_X+2)] + n*p[indl+(_SX*_BGSMAX_X+2)] + f*pfront[1] + b*pback[1] );
  devres[ind-i3+1] = e;
  if (qmaxres) { if (qresnorm) { e=fabs(2.0f*e/(o)); } else { e=fabs(e); } if (res[((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x+1))] < e) { res[((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x+1))] = e ; resind[((threadIdx.y)*(_SX*_BGSMAX_X) + (2*threadIdx.x+1))] = ind - i3 + 1 ; } }
  pfront[0]=pcur[0]; efront[0]=ecur[0];
  pcur[0] =pback[0]; ecur[0] =eback[0];
  pfront[1]=pcur[1]; efront[1]=ecur[1];
  pcur[1] =pback[1]; ecur[1] =eback[1];
  oodzcenfront=oodzcenback;
 }
 if (qmaxres) {
  for (indl=_BGSMAX_X ; indl > 0 ; indl>>=1) {
   if (threadIdx.x < indl) {
    if (res[((threadIdx.y)*(_SX*_BGSMAX_X) + (threadIdx.x))] < res[((threadIdx.y+0)*(_SX*_BGSMAX_X) + (threadIdx.x+indl))] ) { res[((threadIdx.y)*(_SX*_BGSMAX_X) + (threadIdx.x))] = res[((threadIdx.y+0)*(_SX*_BGSMAX_X) + (threadIdx.x+indl))] ; resind[((threadIdx.y)*(_SX*_BGSMAX_X) + (threadIdx.x))] = resind[((threadIdx.y+0)*(_SX*_BGSMAX_X) + (threadIdx.x+indl))] ; };
   }
   __syncthreads();
  }
  for (indl=_BGSMAX_Y/2 ; indl > 0 ; indl>>=1) {
   if (threadIdx.y < indl) {
    if (res[((threadIdx.y)*(_SX*_BGSMAX_X) + (threadIdx.x))] < res[((threadIdx.y+indl)*(_SX*_BGSMAX_X) + (threadIdx.x+0))] ) { res[((threadIdx.y)*(_SX*_BGSMAX_X) + (threadIdx.x))] = res[((threadIdx.y+indl)*(_SX*_BGSMAX_X) + (threadIdx.x+0))] ; resind[((threadIdx.y)*(_SX*_BGSMAX_X) + (threadIdx.x))] = resind[((threadIdx.y+indl)*(_SX*_BGSMAX_X) + (threadIdx.x+0))] ; };
   }
   __syncthreads();
  }
  if (threadIdx.x==0 & threadIdx.y==0 ) {
   maxres[ blockIdx.y * gridDim.x + blockIdx.x ]=res[0];
   imax [ blockIdx.y * gridDim.x + blockIdx.x ]=resind[0];
  }
 }
}
 __global__ void Coarsen_Cuda_3D(
  __CUFLOAT *fine,
  __CUFLOAT *coarse, const __CINT i3, const __CINT nnx, const __CINT nny, const __CINT nnz, const __CINT ibc) {
 coarse[ ( i3 - 1 + ((blockIdx.z * blockDim.z + threadIdx.z )+ibc)*(nnx+2*ibc)*(nny+2*ibc) + ((blockIdx.y * blockDim.y + threadIdx.y )+ibc)*(nnx+2*ibc) + ((blockIdx.x * blockDim.x + threadIdx.x )+ibc) ) ] = 0.1250000000000 * ( fine[( (((blockIdx.z * blockDim.z + threadIdx.z )<<1)+ibc)*(4*(nnx+ibc)*(nny+ibc)) + (((blockIdx.y * blockDim.y + threadIdx.y )<<1)+ibc)*2*(nnx+ibc) + (((blockIdx.x * blockDim.x + threadIdx.x )<<1)+ibc) )] + fine[( (((blockIdx.z * blockDim.z + threadIdx.z )<<1)+ibc)*(4*(nnx+ibc)*(nny+ibc)) + (((blockIdx.y * blockDim.y + threadIdx.y )<<1)+ibc)*2*(nnx+ibc) + (((blockIdx.x * blockDim.x + threadIdx.x )<<1)+1 +ibc) )] + fine[( (((blockIdx.z * blockDim.z + threadIdx.z )<<1)+ibc)*(4*(nnx+ibc)*(nny+ibc)) + (((blockIdx.y * blockDim.y + threadIdx.y )<<1)+1 +ibc)*2*(nnx+ibc) + (((blockIdx.x * blockDim.x + threadIdx.x )<<1)+1 +ibc) )] + fine[( (((blockIdx.z * blockDim.z + threadIdx.z )<<1)+ibc)*(4*(nnx+ibc)*(nny+ibc)) + (((blockIdx.y * blockDim.y + threadIdx.y )<<1)+1 +ibc)*2*(nnx+ibc) + (((blockIdx.x * blockDim.x + threadIdx.x )<<1)+ibc) )] +
                                       fine[( (((blockIdx.z * blockDim.z + threadIdx.z )<<1)+1 +ibc)*(4*(nnx+ibc)*(nny+ibc)) + (((blockIdx.y * blockDim.y + threadIdx.y )<<1)+ibc)*2*(nnx+ibc) + (((blockIdx.x * blockDim.x + threadIdx.x )<<1)+ibc) )] + fine[( (((blockIdx.z * blockDim.z + threadIdx.z )<<1)+1 +ibc)*(4*(nnx+ibc)*(nny+ibc)) + (((blockIdx.y * blockDim.y + threadIdx.y )<<1)+ibc)*2*(nnx+ibc) + (((blockIdx.x * blockDim.x + threadIdx.x )<<1)+1 +ibc) )] + fine[( (((blockIdx.z * blockDim.z + threadIdx.z )<<1)+1 +ibc)*(4*(nnx+ibc)*(nny+ibc)) + (((blockIdx.y * blockDim.y + threadIdx.y )<<1)+1 +ibc)*2*(nnx+ibc) + (((blockIdx.x * blockDim.x + threadIdx.x )<<1)+1 +ibc) )] + fine[( (((blockIdx.z * blockDim.z + threadIdx.z )<<1)+1 +ibc)*(4*(nnx+ibc)*(nny+ibc)) + (((blockIdx.y * blockDim.y + threadIdx.y )<<1)+1 +ibc)*2*(nnx+ibc) + (((blockIdx.x * blockDim.x + threadIdx.x )<<1)+ibc) )] ) ;
}
 __global__ void Refine_Cuda_3D(__CUFLOAT *fine,
 __CUFLOAT *coarse,
 const __CINT i3f, const __CINT i3c, const __CINT nnx, const __CINT nny, const __CINT nnz) {
 __shared__ float clocal[ 2 * ( ( _BRFNE_X + 2 ) * ( _BRFNE_Y + 2 ) ) ];
 float *cfront = clocal ;
 float *cback = cfront + ( ( _BRFNE_X + 2 ) * ( _BRFNE_Y + 2 ) ) ;
 float *cswap ;
 float cwsf, cesf, cwnf, cenf, cwsb, cesb, cwnb, cenb ;
 unsigned int ixc = (blockIdx.x * _BRFNE_X) + threadIdx.x;
 unsigned int iyc = (blockIdx.y * _BRFNE_Y) + threadIdx.y;
 bool qedgeblock;
 for (int k=-1 ; k<=nnz ; k++) {
  cback[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x+1))]=0.0156250000000*coarse[( i3c - 1 + (k+1)*(nnx+2)*(nny+2) + (iyc+1)*(nnx+2) + (ixc+1) )];
  if (threadIdx.y < 2) {
   cback[(((_BRFNE_Y+1)*threadIdx.y)*(_BRFNE_X+2) + (threadIdx.x+1))] = 0.0156250000000*coarse[( i3c - 1 + (k+1)*(nnx+2)*(nny+2) + (iyc - threadIdx.y + (_BRFNE_Y+1)*threadIdx.y)*(nnx+2) + (ixc+1) )];
  }
  if (threadIdx.x < 2) {
   cback[((threadIdx.y+1)*(_BRFNE_X+2) + ((_BRFNE_X+1)*(threadIdx.x%2)))] = 0.0156250000000*coarse[( i3c - 1 + (k+1)*(nnx+2)*(nny+2) + (iyc+1)*(nnx+2) + (ixc - threadIdx.x + (_BRFNE_X+1)*(threadIdx.x%2)) )];
  }
  qedgeblock = ( (blockIdx.x==0 || blockIdx.x == gridDim.x-1) && (blockIdx.y==0 || blockIdx.y == gridDim.y-1) ) ;
  if (qedgeblock) __syncthreads();
  if ( ((threadIdx.x==0) || (threadIdx.x==_BRFNE_X-1) || (ixc==nnx-1) ) && ((threadIdx.y==0) || (threadIdx.y==_BRFNE_Y-1) || (iyc==nny-1))) {
   cback[((threadIdx.y+1-(1-2*(bool)(threadIdx.y)))*(_BRFNE_X+2) + (threadIdx.x+1-(1-2*(bool)(threadIdx.x))))]=0.0156250000000*coarse[( i3c - 1 + (k+1)*(nnx+2)*(nny+2) + (iyc+1-(1-2*(bool)(threadIdx.y)))*(nnx+2) + (ixc+1-(1-2*(bool)(threadIdx.x))) )] ;
   if (qedgeblock) {
    cback[((threadIdx.y+1-(1-2*(bool)(iyc)))*(_BRFNE_X+2) + (threadIdx.x+1-(1-2*(bool)(ixc))))] = cback[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x+1-(1-2*(bool)(ixc))))] + cback[((threadIdx.y+1-(1-2*(bool)(iyc)))*(_BRFNE_X+2) + (threadIdx.x+1))] - cback[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x+1))] ;
   }
  }
  if (k==0) {
   __syncthreads();
   if ( ((ixc==0) || (ixc==nnx-1) ) && ((iyc==0) || (iyc==nny-1)) ) {
    cfront[((threadIdx.y+1-(1-2*(bool)(iyc)))*(_BRFNE_X+2) + (threadIdx.x+1-(1-2*(bool)(ixc))))] = 0.6666666666666 * ( cfront[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x+1-(1-2*(bool)(ixc))))] + cfront[((threadIdx.y+1-(1-2*(bool)(iyc)))*(_BRFNE_X+2) + (threadIdx.x+1))] + cback[((threadIdx.y+1-(1-2*(bool)(iyc)))*(_BRFNE_X+2) + (threadIdx.x+1-(1-2*(bool)(ixc))))] ) - cback[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x+1))] ;
   }
  } else if (k==nnz) {
   __syncthreads();
   if ( ((ixc==0) || (ixc==nnx-1) ) && ((iyc==0) || (iyc==nny-1)) ) {
    cback[((threadIdx.y+1-(1-2*(bool)(iyc)))*(_BRFNE_X+2) + (threadIdx.x+1-(1-2*(bool)(ixc))))] = 0.6666666666666 * ( cback[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x+1-(1-2*(bool)(ixc))))] + cback[((threadIdx.y+1-(1-2*(bool)(iyc)))*(_BRFNE_X+2) + (threadIdx.x+1))] + cfront[((threadIdx.y+1-(1-2*(bool)(iyc)))*(_BRFNE_X+2) + (threadIdx.x+1-(1-2*(bool)(ixc))))] ) - cfront[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x+1))] ;
   }
  }
  __syncthreads() ;
  cwsf=cfront[((threadIdx.y)*(_BRFNE_X+2) + (threadIdx.x))] ;
  cesf=cfront[((threadIdx.y)*(_BRFNE_X+2) + (threadIdx.x+1))] ;
  cenf=cfront[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x+1))] ;
  cwnf=cfront[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x))] ;
  cwsb=cback[((threadIdx.y)*(_BRFNE_X+2) + (threadIdx.x))] ;
  cesb=cback[((threadIdx.y)*(_BRFNE_X+2) + (threadIdx.x+1))] ;
  cenb=cback[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x+1))] ;
  cwnb=cback[((threadIdx.y+1)*(_BRFNE_X+2) + (threadIdx.x))] ;
  cswap =cfront;
  cfront=cback;
  cback =cswap;
  if (k<0) continue;
  fine[( i3f - 1 + (2*k)*4*(nnx+1)*(nny+1) + ((iyc<<1)+1)*2*(nnx+1) + ((ixc<<1)+1) )] += (cwsb + (cesb + cwnb + cwsf)*3.0 + (cenb + cesf + cwnf)*9.0 + cenf*27.0);
  fine[( i3f - 1 + (2*k)*4*(nnx+1)*(nny+1) + ((iyc<<1)+1)*2*(nnx+1) + ((ixc<<1)) )] += (cesb + (cwsb + cenb + cesf)*3.0 + (cwnb + cwsf + cenf)*9.0 + cwnf*27.0);
  fine[( i3f - 1 + (2*k)*4*(nnx+1)*(nny+1) + ((iyc<<1))*2*(nnx+1) + ((ixc<<1)+1) )] += (cwnb + (cenb + cwsb + cwnf)*3.0 + (cesb + cenf + cwsf)*9.0 + cesf*27.0);
  fine[( i3f - 1 + (2*k)*4*(nnx+1)*(nny+1) + ((iyc<<1))*2*(nnx+1) + ((ixc<<1)) )] += (cenb + (cwnb + cesb + cenf)*3.0 + (cwsb + cwnf + cesf)*9.0 + cwsf*27.0);
  fine[( i3f - 1 + (2*k+1)*4*(nnx+1)*(nny+1) + ((iyc<<1)+1)*2*(nnx+1) + ((ixc<<1)+1) )] += (cwsf + (cesf + cwnf + cwsb)*3.0 + (cenf + cesb + cwnb)*9.0 + cenb*27.0);
  fine[( i3f - 1 + (2*k+1)*4*(nnx+1)*(nny+1) + ((iyc<<1)+1)*2*(nnx+1) + ((ixc<<1)) )] += (cesf + (cwsf + cenf + cesb)*3.0 + (cwnf + cwsb + cenb)*9.0 + cwnb*27.0);
  fine[( i3f - 1 + (2*k+1)*4*(nnx+1)*(nny+1) + ((iyc<<1))*2*(nnx+1) + ((ixc<<1)+1) )] += (cwnf + (cenf + cwsf + cwnb)*3.0 + (cesf + cenb + cwsb)*9.0 + cesb*27.0);
  fine[( i3f - 1 + (2*k+1)*4*(nnx+1)*(nny+1) + ((iyc<<1))*2*(nnx+1) + ((ixc<<1)) )] += (cenf + (cwnf + cesf + cenb)*3.0 + (cwsf + cwnb + cesb)*9.0 + cwsb*27.0);
 }
}
extern "C" void AllocDevMem(__CUFLOAT **p, __CINT n) {
 checkCudaErrors(cudaMalloc(p, n*sizeof(__CUFLOAT)));
}
extern "C" void FreeDevMem(__CUFLOAT **p) {
 checkCudaErrors(cudaFree(*p));
}
extern "C" void CopyHostToDevice(__CFLOAT *hostp, __CUFLOAT *devp, __CINT n){
 checkCudaErrors(cudaMemcpy(devp, hostp, sizeof(__CUFLOAT)*n, cudaMemcpyHostToDevice));
}
extern "C" void CopyDeviceToHost(__CFLOAT *hostp, __CUFLOAT *devp, __CINT n){
 checkCudaErrors(cudaMemcpy(hostp, devp, sizeof(__CUFLOAT)*n, cudaMemcpyDeviceToHost));
}
extern "C" void InitDevMem(__CUFLOAT **devp, __CINT i3b, __CINT v, __CINT n) {
 checkCudaErrors(cudaMemset( *devp + i3b - 1, v, n*sizeof(__CUFLOAT)));
}
extern "C" void BindTextures( __CFLOAT *devrhs, __CFLOAT *devkappa, __CFLOAT *deveps,
                              __CFLOAT *devbcw, __CFLOAT *devbce, __CFLOAT *devbcs, __CFLOAT *devbcn, __CFLOAT *devbcf, __CFLOAT *devbcb,
                              const __CINT len3D, const __CINT len3Dbc,
                              const __CINT len2Dyz, const __CINT len2Dxz, const __CINT len2Dxy ) {
}
extern "C" void UnbindTextures() {
}
extern "C" void Residual_Cuda(__CUFLOAT *devres, __CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, __CUFLOAT *devoodx, __CUFLOAT *devoody, __CUFLOAT *devoodz,
                              const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz,
                              const __CINT1 i2d, const __CINT qmaxres, const __CINT qresnorm, __CFLOAT *maxres, __CINT *imax) {
 int nnx=nx-2;
 int nny=ny-2;
 int nblk ;
 float datasize ;
 float *restile = NULL, *drestile = NULL;
 __CINT *imaxtile = NULL, *dimaxtile = NULL;
 if (!i2d) {
  dim3 block(_BGSMAX_X, _BGSMAX_Y);
  dim3 grid( (nnx)/(_SX*_BGSMAX_X) + ( (nnx) % (_SX*_BGSMAX_X) > 0 ) , (nny)/(_SY*_BGSMAX_Y) + ( (nny) % (_SY*_BGSMAX_Y) > 0 ));
  nblk = grid.x * grid.y ;
  datasize = nblk * ( sizeof(__CFLOAT) + sizeof(__CINT) ) ;
  if (qmaxres) {
   checkCudaErrors(cudaMalloc(&drestile, datasize));
   dimaxtile = (__CINT *) drestile + nblk;
   restile = (__CFLOAT *) malloc(datasize);
   imaxtile = (__CINT *) restile + nblk;
  }
  Residual_Cuda_3D<<<grid,block,0>>>(devres, devp, devrhs, deveps, devkappa, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, qmaxres, qresnorm, drestile, dimaxtile);
  if (qmaxres) {
   checkCudaErrors(cudaMemcpy(restile, drestile, datasize, cudaMemcpyDeviceToHost));
   maxres[0]=restile[0];
   imax[0]=imaxtile[0];
   for (unsigned int i=1 ; i < nblk ; i++) {
    if (maxres[0]<restile[i]) {
     maxres[0]=restile[i];
     imax[0]=imaxtile[i];
    }
   }
   free(restile);
   cudaFree(drestile);
  }
 }
}
extern "C" void Coarsen_Cuda(__CUFLOAT *fine, __CUFLOAT *coarse, const __CINT i3, const __CINT nx, const __CINT ny, const __CINT nz, const __CINT1 i2d, const __CINT ibc) {
 int nnx=nx/2;
 int nny=ny/2;
 int nnz=nz/(2-i2d);
 if (!i2d) {
  dim3 block(_BCRSE_X, _BCRSE_Y, _BCRSE_Z);
  dim3 grid( (nnx)/(_BCRSE_X) + ( (nnx) % (_BCRSE_X) > 0 ), (nny)/(_BCRSE_Y) + ( (nny) % (_BCRSE_Y) > 0 ), (nnz)/(_BCRSE_Z) + ( (nnz) % (_BCRSE_Z) > 0 ));
  Coarsen_Cuda_3D<<<grid,block>>>(fine, coarse, i3, nnx, nny, nnz, ibc);
 }
}
extern "C" void Refine_Cuda(__CUFLOAT *fine, __CUFLOAT *coarse, const __CINT i3f, const __CINT i3c, const __CINT nx, const __CINT ny, const __CINT nz, const __CINT1 i2d) {
 int nnx=nx/2;
 int nny=ny/2;
 int nnz=nz/(2-i2d);
 if (!i2d) {
  dim3 block(_BRFNE_X, _BRFNE_Y);
  dim3 grid((nnx)/(_BRFNE_X) + ( (nnx) % (_BRFNE_X) > 0 ), (nny)/(_BRFNE_Y) + ( (nny) % (_BRFNE_Y) > 0 ));
  Refine_Cuda_3D<<<grid,block>>>(fine, coarse, i3f, i3c, nnx, nny, nnz);
 }
}
extern "C" void GaussSeidel_Cuda(__CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, __CUFLOAT *devoodx, __CUFLOAT *devoody, __CUFLOAT *devoodz,
                                 const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz,
                                 const __CFLOAT dt, const __CINT1 i2d, __CINT1 *qpinitzero) {
 int nnx=nx-2;
 int nny=ny-2;
 if (!i2d) {
  dim3 block(_BGSMAX_X, _BGSMAX_Y);
  dim3 grid( (nnx)/(_SX*_BGSMAX_X) + ( (nnx) % (_SX*_BGSMAX_X) > 0 ) , (nny)/(_SY*_BGSMAX_Y) + ( (nny) % (_SY*_BGSMAX_Y) > 0 ));
  Gauss_Seidel_Cuda_3D<<<grid, block, 0>>>(devp, devrhs, deveps, devkappa, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, dt, _REDBLACK, *qpinitzero);
 }
 *qpinitzero=0;
}
extern "C" const int bcwest, bceast, bcnorth, bcsouth, bcback, bcfront ;
extern "C" void ApplyBC_Cuda(__CUFLOAT *devp, __CUFLOAT *devbcw, __CUFLOAT *devbce, __CUFLOAT *devbcn, __CUFLOAT *devbcs, __CUFLOAT *devbcf, __CUFLOAT *devbcb,
                             const __CINT i3b, const __CINT i2, const __CINT j2, const __CINT k2, const __CINT nx, const __CINT ny, const __CINT nz,
                             __CINT *bc_type, __CFLOAT *bc_wgt, const __CINT1 i2d, const __CINT1 qpinitzero_) {
 int nnx=nx-2;
 int nny=ny-2;
 int nnz=nz-2;
 bool qpinitzero =0 ;
 dim3 block, grid ;
 if (!i2d) {
  block.x=1;
  block.y=_BCDIM_Y;
  block.z=_BCDIM_Z;
  grid.x =1;
  grid.y =(nny)/(_BCDIM_Y) + ( (nny) % (_BCDIM_Y) > 0 );
  grid.z =(nnz)/(_BCDIM_Z) + ( (nnz) % (_BCDIM_Z) > 0 );
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, devbcw, devbce, i3b, i2, nx, ny, nz, bcwest, bceast, bc_type[bcwest-1], bc_type[bceast-1], bc_wgt[bcwest-1], bc_wgt[bceast-1], qpinitzero);
  block.x=_BCDIM_X;
  block.y=1;
  grid.x =(nnx)/(_BCDIM_X) + ( (nnx) % (_BCDIM_X) > 0 );
  grid.y =1;
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, devbcn, devbcs, i3b, j2, nx, ny, nz, bcnorth, bcsouth, bc_type[bcnorth-1], bc_type[bcsouth-1], bc_wgt[bcnorth-1], bc_wgt[bcsouth-1], qpinitzero);
  block.y=_BCDIM_Y;
  block.z=1;
  grid.y =(nny)/(_BCDIM_Y) + ( (nny) % (_BCDIM_Y) > 0 );
  grid.z =1;
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, devbcf, devbcb, i3b, k2, nx, ny, nz, bcfront, bcback, bc_type[bcfront-1], bc_type[bcback-1], bc_wgt[bcfront-1], bc_wgt[bcback-1], qpinitzero);
 }
}
