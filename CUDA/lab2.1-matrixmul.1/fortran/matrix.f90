 program matrix
  use, intrinsic :: iso_c_binding
  implicit none
  interface
   subroutine matmul_cuda(M,N,P,my,mx,nx) bind(c,name='matmul_cuda')
   use, intrinsic :: iso_c_binding
   implicit none
   integer, VALUE :: my,mx,nx
   type(C_PTR), VALUE :: M, N, P
   end subroutine matmul_cuda
  end interface

  integer, parameter :: n=3072
  integer, parameter :: mx=n, my=n, nx=n
  real*4, target :: A(my,mx);
  real*4, target :: B(mx,nx);
  real*4, target :: C(my,nx);
  real*4 :: CF(my,nx);
  real*8 :: error(1)

  integer :: i,j

  ! populate matrices:
  do i=1,mx
   do j=1,my
    A(j,i)=i+i-1
   enddo
  enddo

  do i=1,nx
   do j=1,nx
    B(j,i)=i+i-1
   enddo
  enddo

!  C=0d0

  ! call fortran
  write(0,*) 'Calling CUDA matrix multiply ...'
  call matmul_cuda(C_LOC(A), C_LOC(B), C_LOC(C), my, mx, nx)
  error=maxval(abs(C-(matmul(B,A))))
  write(0,*) 'Maximum unsigned difference between CPU/GPU:',error
  ! compare matrices
  !
!  write(0,*) C
!  write(0,*)
!  write(0,*) (matmul(B,A))

  end