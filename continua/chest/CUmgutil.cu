#define _MAX(B,C) (((B)>(C)) ? (B) : (C))
#define _MIN(B,C) (((C)>(B)) ? (B) : (C))

#ifdef __BCTEX
// declare textures
 texture<float> texbcw;
 texture<float> texbce;
 texture<float> texbcs;
 texture<float> texbcn;
 texture<float> texbcf;
 texture<float> texbcb;
#endif
#ifdef __MGTEX
 texture<float> texp;
 texture<float> texrhs;
 texture<float> texeps;
 texture<float> texkappa;
 texture<float> texfine;
 texture<float> texcoarse;
//
#endif

#include <helper_cuda.h>
#include <stdio.h>
#include "mgkernels.cu"
#include "bckernels.cu"
#include "residual.cu"
#include "coarsen.cu"
#include "refine.cu"


extern "C" void AllocDevMem(__CUFLOAT **p, __CINT n) {
// note that p is a double pointer, so that we do not have to take its address
// not sure everything is done correctly
// cudaMalloc(&p, n*sizeof(__CUFLOAT));
 checkCudaErrors(cudaMalloc(p, n*sizeof(__CUFLOAT)));
// checkCudaErrors(cudaMemset(*p, 0, sizeof(float)));
//
}

extern "C" void FreeDevMem(__CUFLOAT **p) {
// cudaFree(p);
 checkCudaErrors(cudaFree(*p));
}

extern "C" void CopyHostToDevice(__CFLOAT *hostp, __CUFLOAT *devp, __CINT n){
 checkCudaErrors(cudaMemcpy(devp, hostp, sizeof(__CUFLOAT)*n, cudaMemcpyHostToDevice));
}

extern "C" void CopyDeviceToHost(__CFLOAT *hostp, __CUFLOAT *devp, __CINT n){
 checkCudaErrors(cudaMemcpy(hostp, devp, sizeof(__CUFLOAT)*n, cudaMemcpyDeviceToHost));
}

extern "C" void InitDevMem(__CUFLOAT **devp, __CINT i3b, __CINT v, __CINT n) { // note that every byte of devp will be set to v
 checkCudaErrors(cudaMemset( *devp + i3b - 1, v, n*sizeof(__CUFLOAT)));
}

extern "C" void BindTextures( __CFLOAT *devrhs, __CFLOAT *devkappa, __CFLOAT *deveps,
                              __CFLOAT *devbcw, __CFLOAT *devbce, __CFLOAT *devbcs, __CFLOAT *devbcn, __CFLOAT *devbcf, __CFLOAT *devbcb, 
                              const __CINT len3D, const __CINT len3Dbc,
                              const __CINT len2Dyz, const __CINT len2Dxz, const __CINT len2Dxy ) {
/*bind to textures
#ifdef __BCTEX
 checkCudaErrors(cudaBindTexture(NULL, texbcw, devbcw, len2Dyz*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbce, devbce, len2Dyz*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbcs, devbcs, len2Dxz*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbcn, devbcn, len2Dxz*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbcf, devbcf, len2Dxy*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbcb, devbcb, len2Dxy*sizeof(float)));
#endif
//
#ifdef __MGTEX
 checkCudaErrors(cudaBindTexture(NULL, texrhs, devrhs, len3D*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texkappa, devkappa, len3D*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texeps, deveps, len3Dbc*sizeof(float)));
#endif
*/
}

extern "C" void UnbindTextures() {
/*unbind textures
#ifdef __BCTEX
 checkCudaErrors(cudaUnbindTexture(texbcw));
 checkCudaErrors(cudaUnbindTexture(texbce));
 checkCudaErrors(cudaUnbindTexture(texbcs));
 checkCudaErrors(cudaUnbindTexture(texbcn));
 checkCudaErrors(cudaUnbindTexture(texbcf));
 checkCudaErrors(cudaUnbindTexture(texbcb));
#endif
#ifdef __MGTEX
 checkCudaErrors(cudaUnbindTexture(texrhs));
 checkCudaErrors(cudaUnbindTexture(texkappa));
 checkCudaErrors(cudaUnbindTexture(texeps));
#endif
*/
}
//=========================================================================================================================================================================//
#define _NBLK(N,n) (N)/(n) + ( (N) % (n) > 0 )
extern "C" void Residual_Cuda(__CUFLOAT *devres, __CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, __CUFLOAT *devoodx, __CUFLOAT *devoody, __CUFLOAT *devoodz,
                              const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz, 
                              const int8_t i2d, const __CINT qmaxres, const __CINT qresnorm, __CFLOAT *maxres, __CINT *imax) {
 int nnx=nx-2; // inner points
 int nny=ny-2;
#ifdef __MGTEX
 int nnz=nz-2;
#endif
 int nblk ;
 float datasize ;
 float  *restile = NULL, *drestile = NULL; // host and device, respectively
 __CINT *imaxtile = NULL, *dimaxtile = NULL;
//
 if (!i2d) {
#ifdef __DSHMEM
// compute block size based on the array sizes
// same code as for GS kernel below
  unsigned int bx=4 ; while ( ( (_SX*bx) <= nnx) && (bx <= _BGSMAX_X) && !(nnx % (_SX*bx))) { bx<<=1; }; bx>>=1;
  unsigned int by=4 ; while ( ( (_SY*by) <= nny) && (by <= _BGSMAX_Y) && !(nny % (_SY*by))) { by<<=1; }; by>>=1;
#define _BX bx
#define _BY by

#define ntx block.x
#define nty block.y
#define sizep      (( _SX * ntx + 2 ) * ( _SY * nty + 2 ))
#define sizeeps    sizep
#define sizeres    (( _SX * ntx     ) * ( _SY * nty     ))
#define sizedxcen  (_SX * ntx + 1)
#define sizedxcor  (_SX * ntx)
#define sizedycen  (_SY * nty + 1)
#define sizedycor  (_SY * nty)
#define _SHMEMSIZE sizeof(float)*(sizep + sizeeps + sizeres + sizeres + sizedxcen + sizedxcor + sizedycen + sizedycor)
#else
#define _BX _BGSMAX_X
#define _BY _BGSMAX_Y
#define _SHMEMSIZE 0
#endif
  dim3 block(_BX, _BY);
  dim3 grid( _NBLK(nnx,_SX*_BX) , _NBLK(nny,_SY*_BY));
  nblk = grid.x * grid.y ;
//
  datasize = nblk * ( sizeof(__CFLOAT) + sizeof(__CINT) ) ;
//
  if (qmaxres) { // allocate memory --  combine allocations for value and index
   checkCudaErrors(cudaMalloc(&drestile, datasize));
   dimaxtile = (__CINT *) drestile + nblk;
   restile  = (__CFLOAT *) malloc(datasize);
   imaxtile = (__CINT *)  restile + nblk;
  }
//
#ifndef __MGTEX
  Residual_Cuda_3D<<<grid,block,_SHMEMSIZE>>>(devres, devp, devrhs, deveps, devkappa, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, qmaxres, qresnorm, drestile, dimaxtile);
#else
  checkCudaErrors(cudaBindTexture(NULL, texrhs, devrhs, (i3-1+nnx*nny*nnz)*sizeof(float)));
  checkCudaErrors(cudaBindTexture(NULL, texkappa, devkappa, (i3-1+nnx*nny*nnz)*sizeof(float)));
  checkCudaErrors(cudaBindTexture(NULL, texp, devp, (i3b-1+nx*ny*nz)*sizeof(float)));
  checkCudaErrors(cudaBindTexture(NULL, texeps, deveps, (i3b-1+nx*ny*nz)*sizeof(float)));
//
  Residual_Cuda_3D<<<grid,block,_SHMEMSIZE>>>(devres, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, qmaxres, qresnorm, drestile, dimaxtile);
//
  checkCudaErrors(cudaUnbindTexture(texrhs));
  checkCudaErrors(cudaUnbindTexture(texkappa));
  checkCudaErrors(cudaUnbindTexture(texp));
  checkCudaErrors(cudaUnbindTexture(texeps));
#endif
  if (qmaxres) { // copy tile residuals to main memory and find the maximum
   checkCudaErrors(cudaMemcpy(restile, drestile, datasize, cudaMemcpyDeviceToHost));
//
   maxres[0]=restile[0];
   imax[0]=imaxtile[0];
//
   for (unsigned int i=1 ; i < nblk ; i++) {
    if (maxres[0]<restile[i]) {
     maxres[0]=restile[i];
     imax[0]=imaxtile[i];
    } //if
   } //for
//
   free(restile);
   cudaFree(drestile);
  } // qmaxres
 }
#undef _SHMEMSIZE
#undef _BX
#undef _BY
}

//=========================================================================================================================================================================//
extern "C" void Coarsen_Cuda(__CUFLOAT *fine, __CUFLOAT *coarse, const __CINT i3, const __CINT nx, const __CINT ny, const __CINT nz, const int8_t i2d, const __CINT ibc) {
// NOTE that the dimensions passed in correspond to the fine grid
 int nnx=nx/2; // coarse grid points (INNER GRID)
 int nny=ny/2;
 int nnz=nz/(2-i2d);
//
 if (!i2d) {
  dim3 block(_BCRSE_X, _BCRSE_Y, _BCRSE_Z);
  dim3 grid( _NBLK(nnx,_BCRSE_X), _NBLK(nny,_BCRSE_Y), _NBLK(nnz,_BCRSE_Z));
#ifndef __MGTEX
  Coarsen_Cuda_3D<<<grid,block>>>(fine, coarse, i3, nnx, nny, nnz, ibc);
#else
  checkCudaErrors(cudaBindTexture(NULL, texfine, fine, (nx+2*ibc)*(ny+2*ibc)*(nz+2*ibc)*sizeof(float)));
  Coarsen_Cuda_3D<<<grid,block>>>(coarse, i3, nnx, nny, nnz, ibc);
  checkCudaErrors(cudaUnbindTexture(texfine));
#endif
 } //i2d
}
//=========================================================================================================================================================================//
extern "C" void Refine_Cuda(__CUFLOAT *fine, __CUFLOAT *coarse, const __CINT i3f, const __CINT i3c, const __CINT nx, const __CINT ny, const __CINT nz, const int8_t i2d) {
// NOTE that the dimensions passed in correspond to the fine grid
//
 int nnx=nx/2; // coarse grid points (INNER GRID)
 int nny=ny/2;
 int nnz=nz/(2-i2d);
//
 if (!i2d) {
  dim3 block(_BRFNE_X, _BRFNE_Y);
  dim3 grid(_NBLK(nnx,_BRFNE_X), _NBLK(nny,_BRFNE_Y));
//#ifndef __MGTEX
  Refine_Cuda_3D<<<grid,block>>>(fine, coarse, i3f, i3c, nnx, nny, nnz);
//#else
//  size_t coarse_byte_offset = sizeof(__CUFLOAT)*(i3c-1) ;
// note : texture memory needs to be aligned, which cannot be guaranteed here
//  checkCudaErrors(cudaBindTexture(NULL, texcoarse, coarse+i3c-1, (nnx+2)*(nny+2)*(nnz+2)*sizeof(float))); // this does not work
//  checkCudaErrors(cudaBindTexture(&coarse_byte_offset, texcoarse, coarse + i3c-1, (nnx+2)*(nny+2)*(nnz+2)*sizeof(float))); // this does not work
//  Refine_Cuda_3D<<<grid,block>>>(fine, i3f, i3c, nnx, nny, nnz);
//  checkCudaErrors(cudaUnbindTexture(texcoarse));
//#endif
 } //i2d
}

//=========================================================================================================================================================================//
extern "C" void GaussSeidel_Cuda(__CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, __CUFLOAT *devoodx, __CUFLOAT *devoody, __CUFLOAT *devoodz,
                                 const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz, 
                                 const __CFLOAT dt, const int8_t i2d, int8_t *qpinitzero) {

 int nnx=nx-2; // inner points
 int nny=ny-2;
#ifdef __MGTEX
 int nnz=nz-2;
#endif
//printf("2D calculation ? %1u\n",i2d);
//
 if (!i2d) {
//  dim3 grid( nnx / (_SX*_BSIZE_X) + ( nnx % (_SX*_BSIZE_X) > 0 ), nny / (_SY*_BSIZE_Y) + ( nnx % (_SY*_BSIZE_Y) > 0 ));
//
#ifdef __DSHMEM
// compute block size based on the array sizes
  unsigned int bx=4 ; while ( ( (_SX*bx) <= nnx) && (bx <= _BGSMAX_X) && !(nnx % (_SX*bx))) { bx<<=1; }; bx>>=1;
  unsigned int by=4 ; while ( ( (_SY*by) <= nny) && (by <= _BGSMAX_Y) && !(nny % (_SY*by))) { by<<=1; }; by>>=1;
#define _BX bx
#define _BY by

#define ntx block.x
#define nty block.y
#define sizep      (( _SX * ntx + 2 ) * ( _SY * nty + 2 ))
#define sizeeps    sizep
#define sizerhs    (( _SX * ntx     ) * ( _SY * nty     ))
#define sizekappa  sizerhs
#define sizedxcen  (_SX * ntx + 1)
#define sizedxcor  (_SX * ntx)
#define sizedycen  (_SY * nty + 1)
#define sizedycor  (_SY * nty)
#define _SHMEMSIZE sizeof(float)*(sizep + sizeeps + sizerhs + sizekappa + sizedxcen + sizedxcor + sizedycen + sizedycor)
#else
#define _BX _BGSMAX_X
#define _BY _BGSMAX_Y
#define _SHMEMSIZE 0
#endif
  dim3 block(_BX, _BY); // static block size
  dim3 grid( _NBLK(nnx,_SX*_BX) , _NBLK(nny,_SY*_BY));

//printf("calling GS-CUDA-3D with the %5d x %5d block grid\n",_BX,_BY);
#ifndef __MGTEX
  Gauss_Seidel_Cuda_3D<<<grid, block, _SHMEMSIZE>>>(devp, devrhs, deveps, devkappa, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, dt, _REDBLACK, *qpinitzero);
#else
// note that the code below binds the array fron the beginning ; might want to only bind the current level data
  checkCudaErrors(cudaBindTexture(NULL, texrhs, devrhs, (i3-1+nnx*nny*nnz)*sizeof(float)));
  checkCudaErrors(cudaBindTexture(NULL, texkappa, devkappa, (i3-1+nnx*nny*nnz)*sizeof(float)));
  checkCudaErrors(cudaBindTexture(NULL, texeps, deveps, (i3b-1+nx*ny*nz)*sizeof(float)));
//
  Gauss_Seidel_Cuda_3D<<<grid, block, _SHMEMSIZE>>>(devp, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, dt, _REDBLACK, *qpinitzero);
//
  checkCudaErrors(cudaUnbindTexture(texrhs));
  checkCudaErrors(cudaUnbindTexture(texkappa));
  checkCudaErrors(cudaUnbindTexture(texeps));
#endif
//  Gauss_Seidel_Cuda_3D<<<grid, block>>>(devp, devrhs, deveps, devkappa, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, dt, _RED);
//  Gauss_Seidel_Cuda_3D<<<grid, block>>>(devp, devrhs, deveps, devkappa, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, dt, _BLACK);
 }
#undef _SHMEMSIZE
 *qpinitzero=0; // this is essential -- needed so that we do not restart from zero every time
}

extern "C" const int bcwest, bceast, bcnorth, bcsouth, bcback, bcfront ; // from FORTRAN

extern "C" void ApplyBC_Cuda(__CUFLOAT *devp, __CUFLOAT *devbcw, __CUFLOAT *devbce, __CUFLOAT *devbcn, __CUFLOAT *devbcs, __CUFLOAT *devbcf, __CUFLOAT *devbcb,
                             const __CINT i3b, const __CINT i2, const __CINT j2, const __CINT k2, const __CINT nx, const __CINT ny, const __CINT nz, 
                             __CINT *bc_type, __CFLOAT *bc_wgt, const int8_t i2d, const int8_t qpinitzero) {
 int nnx=nx-2; // inner points
 int nny=ny-2;
 int nnz=nz-2;

// bool qpinitzero =0 ; // DEBUG

#ifdef __BCTEX
//bind to textures
 checkCudaErrors(cudaBindTexture(NULL, texbcw, devbcw, (i2-1+nny*nnz)*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbce, devbce, (i2-1+nny*nnz)*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbcs, devbcs, (j2-1+nnx*nnz)*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbcn, devbcn, (j2-1+nnx*nnz)*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbcf, devbcf, (k2-1+nnx*nny)*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbcb, devbcb, (k2-1+nnx*nny)*sizeof(float)));
//
#endif

 dim3 block, grid ;
 if (!i2d) {
// x
  block.x=1;
  block.y=_BCDIM_Y;
  block.z=_BCDIM_Z;
  grid.x =1;
  grid.y =_NBLK(nny,_BCDIM_Y);
  grid.z =_NBLK(nnz,_BCDIM_Z);
//
//  Apply_BC_Cuda_3D<<<grid, block>>>(devp, devbcw, i3b, i2, nx, ny, nz, bcwest, bc_type[bcwest-1], bc_wgt[bcwest-1]);
//  Apply_BC_Cuda_3D<<<grid, block>>>(devp, devbce, i3b, i2, nx, ny, nz, bceast, bc_type[bceast-1], bc_wgt[bceast-1]);
// printf(" pinitzero %1d\n",qpinitzero);
#ifndef __BCTEX
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, devbcw, devbce, i3b, i2, nx, ny, nz, bcwest, bceast, bc_type[bcwest-1], bc_type[bceast-1], bc_wgt[bcwest-1], bc_wgt[bceast-1], qpinitzero);
#else
//  checkCudaErrors(cudaBindTexture(NULL, texbc1, devbcw, (i2-1+nny*nnz)*sizeof(float)));
//  checkCudaErrors(cudaBindTexture(NULL, texbc2, devbce, (i2-1+nny*nnz)*sizeof(float)));
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, i3b, i2, nx, ny, nz, bcwest, bceast, bc_type[bcwest-1], bc_type[bceast-1], bc_wgt[bcwest-1], bc_wgt[bceast-1], qpinitzero);
//  checkCudaErrors(cudaUnbindTexture(texbc1));
//  checkCudaErrors(cudaUnbindTexture(texbc2));
#endif
// y
  block.x=_BCDIM_X;
  block.y=1;
//  block.z=_BCDIM_Z;
  grid.x =_NBLK(nnx,_BCDIM_X);
  grid.y =1;
//  grid.z =_NBLK(nnz,_BCDIM_Z);
//
//  Apply_BC_Cuda_3D<<<grid, block>>>(devp, devbcn, i3b, j2, nx, ny, nz, bcnorth, bc_type[bcnorth-1], bc_wgt[bcnorth-1]);
//  Apply_BC_Cuda_3D<<<grid, block>>>(devp, devbcs, i3b, j2, nx, ny, nz, bcsouth, bc_type[bcsouth-1], bc_wgt[bcsouth-1]);
#ifndef __BCTEX
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, devbcn, devbcs, i3b, j2, nx, ny, nz, bcnorth, bcsouth, bc_type[bcnorth-1], bc_type[bcsouth-1], bc_wgt[bcnorth-1], bc_wgt[bcsouth-1], qpinitzero);
#else
//  checkCudaErrors(cudaBindTexture(NULL, texbc1, devbcs, (j2-1+nnx*nnz)*sizeof(float)));
//  checkCudaErrors(cudaBindTexture(NULL, texbc2, devbcn, (j2-1+nnx*nnx)*sizeof(float)));
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, i3b, j2, nx, ny, nz, bcnorth, bcsouth, bc_type[bcnorth-1], bc_type[bcsouth-1], bc_wgt[bcnorth-1], bc_wgt[bcsouth-1], qpinitzero);
//  checkCudaErrors(cudaUnbindTexture(texbc1));
//  checkCudaErrors(cudaUnbindTexture(texbc2));
#endif
// z
//  block.x=_BCDIM_X;
  block.y=_BCDIM_Y;
  block.z=1;
//  grid.x =_NBLK(nnx,_BCDIM_X);
  grid.y =_NBLK(nny,_BCDIM_Y);
  grid.z =1;
//
//  Apply_BC_Cuda_3D<<<grid, block>>>(devp, devbcf, i3b, k2, nx, ny, nz, bcfront, bc_type[bcfront-1], bc_wgt[bcfront-1]);
//  Apply_BC_Cuda_3D<<<grid, block>>>(devp, devbcb, i3b, k2, nx, ny, nz, bcback,  bc_type[bcback-1],  bc_wgt[bcback-1]);
#ifndef __BCTEX
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, devbcf, devbcb, i3b, k2, nx, ny, nz, bcfront, bcback, bc_type[bcfront-1], bc_type[bcback-1], bc_wgt[bcfront-1], bc_wgt[bcback-1], qpinitzero);
#else
//  checkCudaErrors(cudaBindTexture(NULL, texbc1, devbcf, (k2-1+nnx*nny)*sizeof(float)));
//  checkCudaErrors(cudaBindTexture(NULL, texbc2, devbcb, (k2-1+nnx*nny)*sizeof(float)));
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, i3b, k2, nx, ny, nz, bcfront, bcback, bc_type[bcfront-1], bc_type[bcback-1], bc_wgt[bcfront-1], bc_wgt[bcback-1], qpinitzero);
//  checkCudaErrors(cudaUnbindTexture(texbc1));
//  checkCudaErrors(cudaUnbindTexture(texbc2));
#endif

#ifdef __BCTEX
//unbind textures
 checkCudaErrors(cudaUnbindTexture(texbcw));
 checkCudaErrors(cudaUnbindTexture(texbce));
 checkCudaErrors(cudaUnbindTexture(texbcs));
 checkCudaErrors(cudaUnbindTexture(texbcn));
 checkCudaErrors(cudaUnbindTexture(texbcf));
 checkCudaErrors(cudaUnbindTexture(texbcb));
//
#endif
 }
/* translated from fortran  code
  call apply_bc_dnp(allp(i3b),allwest (i2),nnx, nny, nnz, west,  bc_type(west),  bc_wgt(west,level),  q2d)
  call apply_bc_dnp(allp(i3b),alleast (i2),nnx, nny, nnz, east,  bc_type(east),  bc_wgt(east,level),  q2d)
  call apply_bc_dnp(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level), q2d) 
  call apply_bc_dnp(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level), q2d)
  if (.not.q2d) then 
   call apply_bc_dnp(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level), q2d)
   call apply_bc_dnp(allp(i3b),allback (k2),nnx, nny, nnz, back,  bc_type(back),  bc_wgt(back,level),  q2d)
  endif
*/
}
