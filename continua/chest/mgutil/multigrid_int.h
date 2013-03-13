# 1 "multigrid_int.H"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "multigrid_int.H"
interface
 subroutine copy2d(a,b,nx,ny)
 implicit none
 integer :: nx, ny
 real*8, dimension(nx,ny) :: a, b
 end subroutine copy2d
!*****************************************************************************************!
 subroutine copy3d(a,b,nx,ny,nz)
 implicit none
 integer :: nx, ny, nz
 real*8, dimension(nx,ny,nz) :: a, b
 end subroutine copy3d
!*****************************************************************************************!
 subroutine coarsen2d(fine,coarse,nx,ny)
 implicit none
 integer :: nx, ny
 real*8 :: fine(nx,ny), coarse(nx/2,ny/2)
 end subroutine coarsen2d
!*****************************************************************************************!
 subroutine coarsen3d(fine,coarse,nx,ny,nz)
 implicit none
 integer :: nx, ny, nz
 real*8 :: fine(nx,ny,nz), coarse(nx/2,ny/2,nz/2)
 end subroutine coarsen3d
!*****************************************************************************************!
 subroutine refine2d(fine,coarse,nx,ny)
 implicit none
 integer :: nx, ny
 real*8 :: fine(nx+2,ny+2), coarse(nx/2+2,ny/2+2)
 end subroutine refine2d
!*****************************************************************************************!
 subroutine refine3d(fine,coarse,nx,ny,nz)
 implicit none
 integer :: nx, ny, nz
 real*8 :: fine(nx+2,ny+2,nz+2), coarse(nx/2+2,ny/2+2,nz/2+2)
 end subroutine refine3d
!*****************************************************************************************!
 subroutine residual2d(r,p,rhs,e,w,s,n,o,nx,ny)
 implicit none
 integer :: nx, ny
 real*8, dimension(nx,ny) :: p
 real*8, dimension(nx-2,ny-2) :: r,rhs,e,w,s,n,o
 end subroutine residual2d
!*****************************************************************************************!
 subroutine residual3d(r,p,rhs,e,w,s,n,f,b,o,nx,ny,nz)
 implicit none
 integer :: nx, ny, nz
 real*8, dimension(nx,ny,nz) :: p
 real*8, dimension(nx-2,ny-2,nz-2) :: r,rhs,e,w,s,n,f,b,o
 end subroutine residual3d
!*****************************************************************************************!
 subroutine GaussSeidel3d(p,rhs,e,w,s,n,f,b,nx,ny,nz,dt,iter)
 implicit none
 integer :: nx, ny, nz, iter
 real*8 :: dt
 real*8, dimension(nx,ny,nz) :: p
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs,e,w,s,n,f,b ! note that the metrics e--b are variable (because they include eps & kappa; oo=1/o)
 end subroutine GaussSeidel3d
!*****************************************************************************************!
 subroutine GaussSeidel2d(p,rhs,e,w,s,n,nx,ny,dt,iter)
 implicit none
 integer :: nx, ny, iter
 real*8 :: dt
 real*8, dimension(nx,ny) :: p
 real*8, dimension(nx-2,ny-2) :: rhs,e,w,s,n
 end subroutine GaussSeidel2d
!*****************************************************************************************!
 subroutine compute_fd_coef3d(w,e,s,n,f,b,o,eps,kappa,odxcen,odxcor,odycen,odycor,odzcen,odzcor,nx,ny,nz)
 implicit none
 integer :: nx, ny, nz
 real*8, dimension(nx,ny,nz) :: w,e,s,n,f,b,o,eps,kappa
 real*8 :: odxcen(nx+1), odxcor(nx)
 real*8 :: odycen(ny+1), odycor(ny)
 real*8 :: odzcen(nz+1), odzcor(nz)
 end subroutine compute_fd_coef3d
!*****************************************************************************************!
 subroutine compute_fd_coef2d(w,e,s,n,o,eps,kappa,odxcen,odxcor,odycen,odycor,nx,ny)
 implicit none
 integer :: nx, ny
 real*8, dimension(nx,ny) :: w,e,s,n,o,eps,kappa
 real*8 :: odxcen(nx+1), odxcor(nx)
 real*8 :: odycen(ny+1), odycor(ny)
 end subroutine compute_fd_coef2d
end interface
