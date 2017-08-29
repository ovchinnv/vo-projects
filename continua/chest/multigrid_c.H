interface
 subroutine GaussSeidel_C(p,rhs,w,e,s,n,f,b,nx,ny,nz,dt,i2d) bind(c,name='GaussSeidel_C')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 __FCINT, VALUE :: nx, ny, nz
 __FCFLOAT, VALUE :: dt
 TYPE (C_PTR), VALUE :: p
 TYPE (C_PTR), VALUE :: rhs,w,e,s,n,f,b
 int(C_INT8_T), VALUE :: i2d
 end subroutine GaussSeidel_C

#ifdef __CUDA

 subroutine AllocDevMem_C(dxall, dxlen, dyall, dylen, dzall, dzlen, allp, alleps, len3Dbc, allkappa, allrhs, len3D, &
&                         alleast, allwest, len2Dyz, allnorth, allsouth, len2Dxz, allback, allfront, len2Dxy) bind(c,name='AllocDevMem_C')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR), value :: dxall, dyall, dzall, allp, alleps, allkappa, allrhs, alleast, allwest, allnorth, allsouth, allback, allfront
 __FCINT, value :: dxlen, dylen, dzlen, len3Dbc, len3D, len2Dyz, len2Dxz, len2Dxy
 end subroutine AllocDevMem_C

 subroutine CopyHostToDevice(hostp, devp, len)
 use, intrinsic :: iso_c_binding
 __IMPNONE
 __FCFLOAT, dimension(len) :: hostp
 type(C_PTR) :: devp
 __FCINT, value :: len
 end subroutine CopyHostToDevice

 subroutine FreeDevMem_C(dxall, dyall, dzall, allp, alleps, allkappa, allrhs, &
&                         alleast, allwest, allnorth, allsouth, allback, allfront) bind(c,name='FreeDevMem_C')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR), value :: dxall, dyall, dzall, allp, alleps, allkappa, allrhs, alleast, allwest, allnorth, allsouth, allback, allfront
 end subroutine FreeDevMem_C

#endif //CUDA

end interface
