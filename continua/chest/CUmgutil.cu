
#include <helper_cuda.h>
#include "mgkernels.cu"

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

extern "C" void GaussSeidel_Cuda(__CUFLOAT *devp, __CUFLOAT *devrhs, __CUFLOAT *deveps, __CUFLOAT *devkappa, __CUFLOAT *devdx, __CUFLOAT *devdy, __CUFLOAT *devdz,
                                 __CINT i3b, __CINT i3, __CINT i1, __CINT j1, __CINT k1, __CINT nx, __CINT ny, __CINT nz, __CFLOAT dt, int8_t i2d) {

 int nnx=nx-2; // inner points
 int nny=ny-2;
//
 if (!i2d) {
  dim3 block(_BSIZE_X, _BSIZE_Y);
  dim3 grid( nnx / (_SX*_BSIZE_X) + ( nnx % (_SX*_BSIZE_X) > 0 ), nny / (_SY*_BSIZE_Y) + ( nnx % (_SY*_BSIZE_Y) > 0 ));
//
  Gauss_Seidel_Cuda_3D<<<grid, block>>>(devp, devrhs, deveps, devkappa, devdx, devdy, devdz, i3b, i3, i1, j1, k1, nx, ny, nz, dt);
 }
}

extern "C" void ApplyBC_Cuda(__CUFLOAT *devp, __CUFLOAT *devbcw, __CUFLOAT *devbce, __CUFLOAT *devbcn, __CUFLOAT *devbcs, __CUFLOAT *devbcf, __CUFLOAT *devbcb,
                                 __CINT i3b, __CINT i2, __CINT j2, __CINT k2, __CINT nx, __CINT ny, __CINT nz, __CINT *bc_type, __CFLOAT *bc_wgt, int8_t i2d) {

}
