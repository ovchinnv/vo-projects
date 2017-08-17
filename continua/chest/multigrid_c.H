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
end interface
