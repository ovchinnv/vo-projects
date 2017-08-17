interface
 subroutine GaussSeidel_C(p,rhs,w,e,s,n,f,b,nx,ny,nz,dt,i2d) bind(c,name='GaussSeidel_C')
 use, intrinsic :: iso_c_binding
 implicit none
 integer(c_int), VALUE :: nx, ny, nz
 real(c_float), VALUE :: dt
 TYPE (C_PTR), VALUE :: p
 TYPE (C_PTR), VALUE :: rhs,w,e,s,n,f,b
 integer(C_INT8_T), VALUE :: i2d
 end subroutine GaussSeidel_C
end interface
