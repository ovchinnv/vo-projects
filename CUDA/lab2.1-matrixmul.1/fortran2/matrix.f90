 program matrix
  use, intrinsic :: iso_c_binding
  implicit none

  include 'interface.h'

  integer, parameter :: n=2048
  integer, parameter :: mx=n, my=n, nx=n
  real*4, target :: A(my,mx);
  real*4, target :: B(mx,nx);
  real*4, target :: C(my,nx);
  real*4 :: CF(my,nx);

  type(c_ptr) :: deviceA, deviceB, deviceC

  real*8 :: error(1)

  integer :: i,j

  ! populate matrices:
  do i=1,mx
   do j=1,my
    A(j,i)=i+j-1
   enddo
  enddo

  do i=1,nx
   do j=1,nx
    B(j,i)=i+j-1
   enddo
  enddo

!  C=0d0

  ! allocate cuda arrays
  call AllocDevMem_C(deviceA,mx*my);
  call AllocDevMem_C(deviceB,mx*my);
  call AllocDevMem_C(deviceC,mx*my);
  ! copy to device
  call CopyHostToDevice_C(C_LOC(A), deviceA, mx*my);
  call CopyHostToDevice_C(C_LOC(B), deviceB, mx*my);
  ! call computation
  call matmul_cuda(deviceA, deviceB, deviceC, my, mx, nx)
  ! copy result from device
  call CopyDeviceToHost_C(C_LOC(C), deviceC, mx*my)

!  call matmul_cuda((A), (B), (C), deviceA, deviceB, deviceC, my, mx, nx) ! does not work
  if (n<=16) then 
   write(0,'(8F12.5)') C
   write(0,*)
   write(0,'(8F12.5)') matmul(B,A)
  endif
  error=maxval(abs(C-(matmul(B,A))))
  write(0,*) 'Maximum unsigned difference between CPU/GPU:',error
  ! compare matrices

  ! allocate cuda arrays
  call FreeDevMem_C(deviceA);
  call FreeDevMem_C(deviceB);
  call FreeDevMem_C(deviceC);


  end



