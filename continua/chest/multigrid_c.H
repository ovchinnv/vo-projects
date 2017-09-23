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

 subroutine AllocDevMem_C(devpp, len) bind(c,name='AllocDevMem')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR) :: devpp ! double pointer, passed by reference
 __FCINT, VALUE :: len
 end subroutine AllocDevMem_C

 subroutine CopyHostToDevice_C(hostp, devp, len) bind(c,name='CopyHostToDevice')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 __FCFLOAT, dimension(len) :: hostp
 type(C_PTR), VALUE :: devp
 __FCINT, VALUE :: len
 end subroutine CopyHostToDevice_C

 subroutine CopyDeviceToHost_C(hostp, devp, len) bind(c,name='CopyDeviceToHost')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 __FCFLOAT, dimension(len) :: hostp
 type(C_PTR), VALUE :: devp
 __FCINT, VALUE :: len
 end subroutine CopyDeviceToHost_C

 subroutine FreeDevMem_C(devp) bind(c,name='FreeDevMem')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR) :: devp
 end subroutine FreeDevMem_C

 subroutine InitDevMem_C(devp, i3b, v, len) bind(c,name='InitDevMem')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR) :: devp
 __FCINT, VALUE :: i3b
 __FCINT, VALUE :: v
 __FCINT, VALUE :: len
 end subroutine InitDevMem_C

 subroutine BindTextures_C(rhs, kappa, eps, bcw, bce, bcs, bcn, bcb, bcf, len3d, len3dbc, sx, sy, sz) bind(c,name='BindTextures')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR), VALUE :: rhs, kappa, eps, bcw, bce, bcs, bcn, bcb, bcf
 __FCINT, VALUE :: len3d, len3dbc, sx, sy, sz
 end subroutine BindTextures_C

 subroutine UnbindTextures_C() bind(c,name='UnbindTextures')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 end subroutine UnbindTextures_C

 subroutine GaussSeidel_Cuda(CUallp,CUallrhs,CUalleps,CUallkappa,CUdxall,CUdyall,CUdzall,i3b,i3,i1,j1,k1,nx,ny,nz,omega,i2d,qpinitzero) bind(c,name='GaussSeidel_Cuda')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR), VALUE :: CUallp, CUallrhs, CUalleps, CUallkappa, CUdxall, CUdyall, CUdzall
 __FCINT, VALUE :: i3b, i3 ! indices into data arrays (with and without ghostpoint layer, resp.)
 __FCINT, VALUE :: i1, j1, k1 ! indices into metric arrays
 __FCINT, VALUE :: nx, ny, nz ! array dimensions, incl. ghost points
 __FCFLOAT, VALUE :: omega ! overrelaxation parameter
 int(C_INT8_T), VALUE :: i2d ! whether this is a 2D calculation
 int(C_INT8_T) :: qpinitzero ! whether to initialize solution to zero
 end subroutine GaussSeidel_Cuda

 subroutine Residual_Cuda(CUres,CUallp,CUallrhs,CUalleps,CUallkappa,CUdxall,CUdyall,CUdzall,i3b,i3,i1,j1,k1,nx,ny,nz,i2d,qmaxres,qresnorm,curr_residual,imax) bind(c,name='Residual_Cuda')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR), VALUE :: CUres, CUallp, CUallrhs, CUalleps, CUallkappa, CUdxall, CUdyall, CUdzall
 __FCINT, VALUE :: i3b, i3 ! indices into data arrays (with and without ghostpoint layer, resp.)
 __FCINT, VALUE :: i1, j1, k1 ! indices into metric arrays
 __FCINT, VALUE :: nx, ny, nz ! array dimensions, incl. ghost points
 __FCFLOAT, intent(in) :: curr_residual
 __FCINT, intent(in) :: imax
 __FCINT, VALUE :: qmaxres ! whether to compute maximum residual
 __FCINT, VALUE :: qresnorm ! whether to normalize the residual by central coefficient (o)
 int(C_INT8_T), VALUE :: i2d
 end subroutine Residual_Cuda

 subroutine Coarsen_Cuda(CUfine,CUcoarse,i3,nx,ny,nz,i2d,ibc) bind(c,name='Coarsen_Cuda')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR), VALUE :: CUfine, CUcoarse
 __FCINT, VALUE :: i3 ! index into rhs array
 __FCINT, VALUE :: nx, ny, nz ! fine array dimensions, _not_ incl. ghost points
 int(C_INT8_T), VALUE :: i2d
 __FCINT, VALUE :: ibc ! whether arrays include ghostpoints (notwithstanding the inner dimensions nx/ny/nz)
 end subroutine Coarsen_Cuda

 subroutine Refine_Cuda(CUfine,CUcoarse,i3f,i3c,nx,ny,nz,i2d) bind(c,name='Refine_Cuda')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR), VALUE :: CUfine, CUcoarse
 __FCINT, VALUE :: i3f, i3c ! indices into fine and coarse arrays, respectively
 __FCINT, VALUE :: nx, ny, nz ! fine array dimensions, _not_ incl. ghost points
 int(C_INT8_T), VALUE :: i2d
 end subroutine Refine_Cuda

 subroutine ApplyBC_Cuda(CUallp,CUallwest,CUalleast,CUallnorth,CUallsouth,CUallfront,CUallback,i3b,i2,j2,k2,nx,ny,nz,bc_type,bc_wgt,i2d,qpinitzero) bind(c,name='ApplyBC_Cuda')
 use, intrinsic :: iso_c_binding
 __IMPNONE
 type(C_PTR), VALUE :: CUallp,CUallwest,CUalleast,CUallnorth,CUallsouth,CUallfront,CUallback
 __FCINT, VALUE :: i3b ! index into data arrays (with ghostpoint layer)
 __FCINT, VALUE :: i2, j2, k2 ! indices into boundary condition
 __FCINT, VALUE :: nx, ny, nz ! array dimensions, incl. ghost points
 __FCINT, dimension(6) :: bc_type ! boundary condition codes
 __FCFLOAT, dimension(6) :: bc_wgt ! boundary condition values at current level only
 int(C_INT8_T), VALUE :: i2d
 int(C_INT8_T), VALUE :: qpinitzero
 end subroutine ApplyBC_Cuda


#endif //CUDA

end interface
