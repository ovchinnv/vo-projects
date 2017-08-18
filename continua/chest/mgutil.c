// C functions for CHEST
// to bridge fortran/CPU & Cuda functionality
// some of the routines are for testing or development only
//
#include <math.h>
#include <stdint.h>
//#include <stdio.h>

void GaussSeidel_C(__CFLOAT *p, __CFLOAT *rhs, __CFLOAT *w, __CFLOAT *e, __CFLOAT *s, __CFLOAT *n, __CFLOAT *f, __CFLOAT *b,
                    __CINT nx, __CINT ny, __CINT nz, __CFLOAT dt, int8_t i2d) {
#define IO(i,j,k) ((((k)-1)*(ny) + ((j)-1))*(nx) + i-1)
#define II(i,j,k) ((((k)-1)*(nny) + ((j)-1))*(nnx) + i-1)
// translated from f90
// note that nx,ny,nz include ghost points
// note that rhs passed into this routine is normalized by o (see below)
// float, dimension(nx,ny,nz) :: p 
// float, dimension(nx-2,ny-2,nz-2) :: rhs,w,e,s,n,f,b // note that the metrics e--b are variable (because they include eps & kappa; oo=1/o)
 __CINT nnx, nny, nnz;
 __CINT nxp, nyp, nzp;
 __CINT i, j, k, m, im, jm, km, ii, io;
 __CFLOAT pw, pe, po;

 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
//
 if (i2d) { //2D
//
//##################################################################################
// NOTE : index multiplications do not make any difference in speed;
// they are useless, so I favor full index notation as in the (untouched) 3D case below
  io=IO(2,2,2);
  ii=II(1,1,1);
  k=2; km=k-1;
  for (j=2;j<=nyp;j++){
   jm=j-1;
   po=p[io]; pw=p[io-1];
   for (i=2;i<=nxp;i++){
    im=i-1;
    pe=p[io+1];
    po=po - dt * (             w[ii]*pw       + e[ii]*pe +\
                               s[ii]*p[io-nx] + n[ii]*p[io+nx] +\
                               po - rhs[ii] );
    p[io]=po;
    pw=po; po=pe;

    io+=1;
    ii+=1;
   };
   io+=2; // skip over ghost points
  };
//
 } else { // 3D
//##################################################################################
 for (k=2;k<=nzp;k++){    km=k-1;
  for (j=2;j<=nyp;j++){   jm=j-1;
//
  po=p[IO(2,j,k)]; pw=p[IO(1,j,k)];
   for (i=2;i<=nxp;i++){    im=i-1;
    pe=p[IO(i+1,j,k)];
    po=po - dt * (             w[II(im,jm,km)]*pw + e[II(im,jm,km)]*pe +\
                               s[II(im,jm,km)]*p[IO(i,jm,k)] + n[II(im,jm,km)]*p[IO(i,j+1,k)] +\
                               f[II(im,jm,km)]*p[IO(i,j,km)] + b[II(im,jm,km)]*p[IO(i,j,k+1)] +\
                               po - rhs[II(im,jm,km)] );
    p[IO(i,j,k)]=po;
    pw=po; po=pe;
   };
  };
 }; // k
//
 }; // q2d
#undef IO
#undef II
} // GaussSeidel
//*****************************************************************************************//
