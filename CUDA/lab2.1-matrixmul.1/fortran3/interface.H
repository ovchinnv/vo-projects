interface

 subroutine matmul_cuda(devM,devN,devP,my,mx,nx) bind(c,name='matmul_cuda')
 use, intrinsic :: iso_c_binding
 implicit none
 integer, VALUE :: my,mx,nx
 type(C_PTR), VALUE :: devM, devN, devP
 end subroutine matmul_cuda

 subroutine AllocDevMem_C(devpp, len) bind(c,name='AllocDevMem')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR) :: devpp ! double pointer, passed by reference
 __FCINT, VALUE :: len
 end subroutine AllocDevMem_C

 subroutine CopyHostToDevice_C(hostp, devp, len) bind(c,name='CopyHostToDevice')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 __FCFLOAT :: hostp(len)
 type(C_PTR), VALUE :: devp
 __FCINT, VALUE :: len
 end subroutine CopyHostToDevice_C

 subroutine CopyDeviceToHost_C(hostp, devp, len) bind(c,name='CopyDeviceToHost')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 __FCFLOAT, dimension(len) :: hostp ! this type of addressing does not work
! type(C_PTR), VALUE :: hostp
 type(C_PTR), VALUE :: devp
 __FCINT, VALUE :: len
 end subroutine CopyDeviceToHost_C

 subroutine FreeDevMem_C(devp) bind(c,name='FreeDevMem')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR) :: devp
 end subroutine FreeDevMem_C

end interface
