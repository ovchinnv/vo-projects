
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
//
#endif

#include <helper_cuda.h>
#include "mgkernels.cu"
#include "bckernels.cu"
#include "residual.cu"


extern "C" void AllocDevMem(__CUFLOAT **p, __CINT n) {
// note that p is a double pointer, so that we do not have to take its address
// not sure everything is done correctly
// cudaMalloc(&p, n*sizeof(__CUFLOAT));
 checkCudaErrors(cudaMalloc(p, n*sizeof(__CUFLOAT)));
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
  float  *restile = NULL, *drestile = NULL; // host and device, respectively
  __CINT *imaxtile = NULL, *dimaxtile = NULL;
//
 if (!i2d) {
  dim3 block(_BSIZE_X, _BSIZE_Y);
  dim3 grid( _NBLK(nnx,_BSIZE_X) , _NBLK(nny,_BSIZE_Y));
//
  if (qmaxres) { // allocate memory
   checkCudaErrors(cudaMalloc(&drestile, grid.x * grid.y *sizeof(__CUFLOAT)));
   checkCudaErrors(cudaMalloc(&dimaxtile, grid.x * grid.y *sizeof(__CINT)));
   restile  = (__CUFLOAT *) malloc(grid.x * grid.y *sizeof(__CUFLOAT));
   imaxtile = (__CINT *)    malloc(grid.x * grid.y *sizeof(__CINT));
  }
//
#ifndef __MGTEX
  Residual_Cuda_3D<<<grid, block>>>(devres, devp, devrhs, deveps, devkappa, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, qmaxres, qresnorm, drestile, dimaxtile);
#else
  checkCudaErrors(cudaBindTexture(NULL, texrhs, devrhs, (i3-1+nnx*nny*nnz)*sizeof(float)));
  checkCudaErrors(cudaBindTexture(NULL, texkappa, devkappa, (i3-1+nnx*nny*nnz)*sizeof(float)));
  checkCudaErrors(cudaBindTexture(NULL, texp, devp, (i3b-1+nx*ny*nz)*sizeof(float)));
  checkCudaErrors(cudaBindTexture(NULL, texeps, deveps, (i3b-1+nx*ny*nz)*sizeof(float)));
//
  Residual_Cuda_3D<<<grid, block>>>(devres, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, qmaxres, qresnorm, drestile, dimaxtile);
//
  checkCudaErrors(cudaUnbindTexture(texrhs));
  checkCudaErrors(cudaUnbindTexture(texkappa));
  checkCudaErrors(cudaUnbindTexture(texp));
  checkCudaErrors(cudaUnbindTexture(texeps));
#endif
  if (qmaxres) { // copy tile residuals to main memory and find the maximum
   checkCudaErrors(cudaMemcpy(restile, drestile, grid.x * grid.y * sizeof(__CUFLOAT), cudaMemcpyDeviceToHost));
   checkCudaErrors(cudaMemcpy(imaxtile, dimaxtile, grid.x * grid.y *sizeof(__CINT), cudaMemcpyDeviceToHost)); // maybe combine into restile ?
//
   maxres[0]=restile[0];
   imax[0]=imaxtile[0];
//
   for (unsigned int i=1 ; i < grid.x * grid.y ; i++) {
    if (maxres[0]<restile[i]) {
     maxres[0]=restile[i];
     imax[0]=imaxtile[i];
    } //if
   } //for
//
   free(restile);
   free(imaxtile);
   cudaFree(drestile);
   cudaFree(dimaxtile);
  }
 }
}


extern "C" void GaussSeidel_Cuda(__CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, __CUFLOAT *devoodx, __CUFLOAT *devoody, __CUFLOAT *devoodz,
                                 const __CINT i3b, const __CINT i3, const __CINT i1, const __CINT j1, const __CINT k1, const __CINT nx, const __CINT ny, const __CINT nz, 
                                 const __CFLOAT dt, const int8_t i2d) {

 int nnx=nx-2; // inner points
 int nny=ny-2;
#ifdef __MGTEX
 int nnz=nz-2;
#endif
//
 if (!i2d) {
  dim3 block(_BSIZE_X, _BSIZE_Y);
//  dim3 grid( nnx / (_SX*_BSIZE_X) + ( nnx % (_SX*_BSIZE_X) > 0 ), nny / (_SY*_BSIZE_Y) + ( nnx % (_SY*_BSIZE_Y) > 0 ));
  dim3 grid( _NBLK(nnx,_SX*_BSIZE_X) , _NBLK(nny,_SY*_BSIZE_Y));
//
#ifndef __MGTEX
  Gauss_Seidel_Cuda_3D<<<grid, block>>>(devp, devrhs, deveps, devkappa, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, dt, _REDBLACK);
#else
  checkCudaErrors(cudaBindTexture(NULL, texrhs, devrhs, (i3-1+nnx*nny*nnz)*sizeof(float)));
  checkCudaErrors(cudaBindTexture(NULL, texkappa, devkappa, (i3-1+nnx*nny*nnz)*sizeof(float)));
  checkCudaErrors(cudaBindTexture(NULL, texeps, deveps, (i3b-1+nx*ny*nz)*sizeof(float)));
//
  Gauss_Seidel_Cuda_3D<<<grid, block>>>(devp, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, dt, _REDBLACK);
//
  checkCudaErrors(cudaUnbindTexture(texrhs));
  checkCudaErrors(cudaUnbindTexture(texkappa));
  checkCudaErrors(cudaUnbindTexture(texeps));
#endif
//  Gauss_Seidel_Cuda_3D<<<grid, block>>>(devp, devrhs, deveps, devkappa, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, dt, _RED);
//  Gauss_Seidel_Cuda_3D<<<grid, block>>>(devp, devrhs, deveps, devkappa, devoodx, devoody, devoodz, i3b, i3, i1, j1, k1, nx, ny, nz, dt, _BLACK);
 }
}

extern "C" const int bcwest, bceast, bcnorth, bcsouth, bcback, bcfront ; // from FORTRAN

extern "C" void ApplyBC_Cuda(__CUFLOAT *devp, __CUFLOAT *devbcw, __CUFLOAT *devbce, __CUFLOAT *devbcn, __CUFLOAT *devbcs, __CUFLOAT *devbcf, __CUFLOAT *devbcb,
                             const __CINT i3b, const __CINT i2, const __CINT j2, const __CINT k2, const __CINT nx, const __CINT ny, const __CINT nz, 
                             __CINT *bc_type, __CFLOAT *bc_wgt, const int8_t i2d) {
 int nnx=nx-2; // inner points
 int nny=ny-2;
 int nnz=nz-2;

#ifdef __BCTEX
//bind to textures
 checkCudaErrors(cudaBindTexture(NULL, texbcw, devbcw, (i2-1+nny*nnz)*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbce, devbce, (i2-1+nny*nnz)*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbcs, devbcs, (j2-1+nnx*nnz)*sizeof(float)));
 checkCudaErrors(cudaBindTexture(NULL, texbcn, devbcn, (j2-1+nnx*nnx)*sizeof(float)));
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
#ifndef __BCTEX
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, devbcw, devbce, i3b, i2, nx, ny, nz, bcwest, bceast, bc_type[bcwest-1], bc_type[bceast-1], bc_wgt[bcwest-1], bc_wgt[bceast-1]);
#else
//  checkCudaErrors(cudaBindTexture(NULL, texbc1, devbcw, (i2-1+nny*nnz)*sizeof(float)));
//  checkCudaErrors(cudaBindTexture(NULL, texbc2, devbce, (i2-1+nny*nnz)*sizeof(float)));
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, i3b, i2, nx, ny, nz, bcwest, bceast, bc_type[bcwest-1], bc_type[bceast-1], bc_wgt[bcwest-1], bc_wgt[bceast-1]);
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
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, devbcn, devbcs, i3b, j2, nx, ny, nz, bcnorth, bcsouth, bc_type[bcnorth-1], bc_type[bcsouth-1], bc_wgt[bcnorth-1], bc_wgt[bcsouth-1]);
#else
//  checkCudaErrors(cudaBindTexture(NULL, texbc1, devbcs, (j2-1+nnx*nnz)*sizeof(float)));
//  checkCudaErrors(cudaBindTexture(NULL, texbc2, devbcn, (j2-1+nnx*nnx)*sizeof(float)));
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, i3b, j2, nx, ny, nz, bcnorth, bcsouth, bc_type[bcnorth-1], bc_type[bcsouth-1], bc_wgt[bcnorth-1], bc_wgt[bcsouth-1]);
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
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, devbcf, devbcb, i3b, k2, nx, ny, nz, bcfront, bcback, bc_type[bcfront-1], bc_type[bcback-1], bc_wgt[bcfront-1], bc_wgt[bcback-1]);
#else
//  checkCudaErrors(cudaBindTexture(NULL, texbc1, devbcf, (k2-1+nnx*nny)*sizeof(float)));
//  checkCudaErrors(cudaBindTexture(NULL, texbc2, devbcb, (k2-1+nnx*nny)*sizeof(float)));
  Apply_BC_Cuda_3Dx2<<<grid, block>>>(devp, i3b, k2, nx, ny, nz, bcfront, bcback, bc_type[bcfront-1], bc_type[bcback-1], bc_wgt[bcfront-1], bc_wgt[bcback-1]);
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
