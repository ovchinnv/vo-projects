// NOTE : constants are not defined on device, so I will "hack" the code and use preprocessor macros to match fortran
//
#define bcleft	1
#define bcright	2
#define bcbottom 3
#define bctop	4
#define bcfront	5
#define bcback	6

#define bcdirichlet	1
#define bcneumann	2
#define bcperiodic	3
#define bcdirichletg	4
#define bcneumanng	5

/*
      int, parameter, public :: left=1, right=2, bottom=3, top=4, front=5, back=6
!      int, parameter, public :: west=1, east=2, south=3,  north=4
      int, parameter, public :: west=1, east=2, south=3,  north=4
      int, parameter, public :: bcx0=1, bcx1=2, bcy0=3, bcy1=4, bcz0=5, bcz1=6
!     corresponding boundary names:
      character(len=4), parameter :: bc_names(numbc) = (/'BCX0','BCX1','BCY0','BCY1','BCZ0','BCZ1'/);
      int, public, parameter :: dirichlet=1,  &
     &                          neumann=2,    &
     &                          periodic=3,   &
     &                          dirichletg=4, &
     &                          neumanng=5,   &      ! boundary condition codes
*/

//extern "C" const int bcleft, bcright, bctop, bcbottom, bcback, bcfront ; // boundary definitions from FORTRAN
//extern "C" const int bcdirichlet, bcdirichletg, bcneumann, bcneumanng, bcperiodic ; // boundary condition definitions

 __global__ void  Apply_BC_Cuda_3D(__CUFLOAT *devu, __CUFLOAT *devg, const __CINT i3b, const __CINT ind2, const __CINT nx, const __CINT ny, const __CINT nz,
                                   const __CINT boundary, const __CINT bctype, const __CUFLOAT wgt) {
// translated from FORTRAN apply_bc_dnp
// g has inner points only ; u of course has ghost points

      int ix = blockIdx.x*_BCDIM_X+threadIdx.x ;
      int iy = blockIdx.y*_BCDIM_Y+threadIdx.y ;
      int iz = blockIdx.z*_BCDIM_Z+threadIdx.z ;

// define index functions
#define IOG(i,j,k)  ( i3b - 1 + (k)*(nx*ny) + (j)*(nx) + (i) )
#define I2D(i,j,k,nnx,nny,nnz) ( ind2 - 1 + (k)*((nnx)*(nny)) + (j)*(nnx) + (i) )

      int ib, ie, jb, je, kb, ke ;
      int di, dj, dk ;
      int dir ;
//      float uval ;
//      float gval ;
//
//    preset indices: zero-based rather than one-based as in FORTRAN
//
      ib=1; ie=nx-2; // limiting indices of inner points 
      jb=1; je=ny-2;
      kb=1; ke=nz-2;
//
      if (boundary==bcleft) {

         ib--; ie=ib; di=+1; dj=0; dk=0; dir=di;

        } else if (boundary==bcright) {

         ie++; ib=ie; di=-1; dj=0; dk=0; dir=di;

        } else if (boundary==bcbottom) {

         jb--; je=jb; dj=+1; di=0; dk=0; dir=dj;

        } else if (boundary==bctop) {

         je++; jb=je; dj=-1; di=0; dk=0; dir=dj;

        } else if (boundary==bcfront) {

         kb--; ke=kb; dk=+1; di=0; dj=0; dir=dk;

        } else if (boundary==bcback) {

         ke++; kb=ke; dk=-1; di=0; dj=0; dir=dk;

        } else {

         ke=0; // set invalid indices so that no loops are executed; better than early return
      }
//************************************
#define         gval devg[I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1)]
#define         uval devu[IOG(ix+ib+di,iy+jb+dj,iz+kb+dk)]
      if (bctype==bcdirichlet || bctype==bcdirichletg) {
//         gval=devg[I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1)];
//         uval=devu[IOG(ix+ib+di,iy+jb+dj,iz+kb+dk)];
         devu[IOG(ix+ib,iy+jb,iz+kb)] = uval + wgt * (gval - uval);
//***************************
        } else if (bctype==bcneumann || bctype==bcneumanng) {
//         gval=devg[I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1)];
//         uval=devu[IOG(ix+ib+di,iy+jb+dj,iz+kb+dk)];
         devu[IOG(ix+ib,iy+jb,iz+kb)] = uval - dir * wgt * gval;
//***************************
        } else if (bctype==bcperiodic ) {
         dk*=(nz-2);     dj*=(ny-2);    di*=(nx-2); // wrap around
//         uval=devu[IOG(ix+ib+di,iy+jb+dj,iz+kb+dk)];
         devu[IOG(ix+ib,iy+jb,iz+kb)] = uval ;
//***************************
      } // bctype
//***********************************************************************************
//      __syncthreads(); // should not be needed
#undef gval
#undef uval
 }


//combined 2-bc kernel
#ifndef __BCTEX
 __global__ void  Apply_BC_Cuda_3Dx2(__CUFLOAT *devu, __CUFLOAT *devg, __CUFLOAT *devg2, const __CINT i3b, const __CINT ind2, const __CINT nx, const __CINT ny, const __CINT nz,
#else
 __global__ void  Apply_BC_Cuda_3Dx2(__CUFLOAT *devu, const __CINT i3b, const __CINT ind2, const __CINT nx, const __CINT ny, const __CINT nz,
#endif
                                   const __CINT boundary, const __CINT boundary2, const __CINT bctype, const __CINT bctype2, const __CFLOAT wgt, const __CUFLOAT wgt2) {
// translated from FORTRAN apply_bc_dnp
// g has inner points only ; u of course has ghost points
// this version applies BC to two OPPOSING sides simultaneously ; user is responsible for passing in the correct boundaries

      int ix = blockIdx.x*_BCDIM_X+threadIdx.x ;
      int iy = blockIdx.y*_BCDIM_Y+threadIdx.y ;
      int iz = blockIdx.z*_BCDIM_Z+threadIdx.z ;

// define index functions
#define IOG(i,j,k)  ( i3b - 1 + (k)*(nx*ny) + (j)*(nx) + (i) )
#define I2D(i,j,k,nnx,nny,nnz) ( ind2 - 1 + (k)*((nnx)*(nny)) + (j)*(nnx) + (i) )

      int ib, ie, jb, je, kb, ke ;
      int di, dj, dk ;
      int dir ;
#ifdef __BCTEX
      float gval = 0.f;
#endif
//
//    preset indices: zero-based rather than one-based as in FORTRAN
//
      ib=1; ie=nx-2; // limiting indices of inner points 
      jb=1; je=ny-2;
      kb=1; ke=nz-2;
//
      if (boundary==bcleft) {

         ib--; ie=ib; di=+1; dj=0; dk=0; dir=di;
#ifdef __BCTEX
         gval=tex1Dfetch(texbcw, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else if (boundary==bcright) {

         ie++; ib=ie; di=-1; dj=0; dk=0; dir=di;
#ifdef __BCTEX
         gval=tex1Dfetch(texbce, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else if (boundary==bcbottom) {

         jb--; je=jb; dj=+1; di=0; dk=0; dir=dj;
#ifdef __BCTEX
         gval=tex1Dfetch(texbcs, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else if (boundary==bctop) {

         je++; jb=je; dj=-1; di=0; dk=0; dir=dj;
#ifdef __BCTEX
         gval=tex1Dfetch(texbcn, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else if (boundary==bcfront) {

         kb--; ke=kb; dk=+1; di=0; dj=0; dir=dk;
#ifdef __BCTEX
         gval=tex1Dfetch(texbcf, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else if (boundary==bcback) {

         ke++; kb=ke; dk=-1; di=0; dj=0; dir=dk;
#ifdef __BCTEX
         gval=tex1Dfetch(texbcb, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
      }
//************************************
#ifndef __BCTEX
#define         gval devg[I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1)]
#endif
#define         uval devu[IOG(ix+ib+di,iy+jb+dj,iz+kb+dk)]
      if (bctype==bcdirichlet || bctype==bcdirichletg) {
         devu[IOG(ix+ib,iy+jb,iz+kb)] = uval + wgt * (gval - uval);
//***************************
        } else if (bctype==bcneumann || bctype==bcneumanng) {
         devu[IOG(ix+ib,iy+jb,iz+kb)] = uval - dir * wgt * gval;
//***************************
        } else if (bctype==bcperiodic ) {
         dk*=(nz-2);     dj*=(ny-2);    di*=(nx-2); // wrap around
         devu[IOG(ix+ib,iy+jb,iz+kb)] = uval;
//***************************
      } // bctype
//***********************************************************************************
//    reset indices:
      ib=1; ie=nx-2;
      jb=1; je=ny-2;
      kb=1; ke=nz-2;
//
      if (boundary2==bcleft) {

         ib--; ie=ib; di=+1; dj=0; dk=0; dir=di;
#ifdef __BCTEX
         gval=tex1Dfetch(texbcw, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else if (boundary2==bcright) {

         ie++; ib=ie; di=-1; dj=0; dk=0; dir=di;
#ifdef __BCTEX
         gval=tex1Dfetch(texbce, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else if (boundary2==bcbottom) {

         jb--; je=jb; dj=+1; di=0; dk=0; dir=dj;
#ifdef __BCTEX
         gval=tex1Dfetch(texbcs, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else if (boundary2==bctop) {

         je++; jb=je; dj=-1; di=0; dk=0; dir=dj;
#ifdef __BCTEX
         gval=tex1Dfetch(texbcn, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else if (boundary2==bcfront) {

         kb--; ke=kb; dk=+1; di=0; dj=0; dir=dk;
#ifdef __BCTEX
         gval=tex1Dfetch(texbcf, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else if (boundary2==bcback) {

         ke++; kb=ke; dk=-1; di=0; dj=0; dir=dk;
#ifdef __BCTEX
         gval=tex1Dfetch(texbcb, I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1));
#endif
        } else {

         gval=0.f;
      }
//************************************
#ifndef __BCTEX
#undef gval
#define         gval devg2[I2D(ix,iy,iz,ie-ib+1,je-jb+1,ke-kb+1)]
#endif
      if (bctype2==bcdirichlet || bctype2==bcdirichletg) {
         devu[IOG(ix+ib,iy+jb,iz+kb)] = uval + wgt2 * (gval - uval);
//***************************
        } else if (bctype2==bcneumann || bctype2==bcneumanng) {
         devu[IOG(ix+ib,iy+jb,iz+kb)] = uval - dir * wgt2 * gval;
//***************************
        } else if (bctype2==bcperiodic ) {
         dk*=(nz-2);     dj*=(ny-2);    di*=(nx-2); // wrap around
         devu[IOG(ix+ib,iy+jb,iz+kb)] = uval;
//***************************
      } // bctype2
//***********************************************************************************
 }


#undef IOG
#undef I2D
#undef bcleft
#undef bcright
#undef bcbottom
#undef bctop
#undef bcfront
#undef bcback

#undef bcdirichlet
#undef bcneumann
#undef bcperiodic
#undef bcdirichletg
#undef bcneumanng
