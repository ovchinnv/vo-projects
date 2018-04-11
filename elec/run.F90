 program elec
 
 use vars
 
 float totcharge
 
 call setup()
 
! filter charge density
 __OUT('Computing smoothed charge density');
 call filt3p()
!
 totcharge=sum(rho(1:nx-1, 1:ny-1, 1:nz-1))
  __OUT('Total smoothed charge is ', totcharge) 
!
! compute laplacian for FD Poisson equation
 lap=-(rho-totcharge/((nx-1)*(ny-1)*(nz-1)))/eps ;

! evaluate long-range energy
 call direct3p()
! call direct3p_nosym()
!
 __OUT('Short-range electrostatic potential energy is: ',elsr)
 __OUT('Long-range electrostatic potential energy is: ',el)
 __OUT('Self-energy part of long-range energy is: ',el_self)
 __OUT('Total electrostatic potential energy is: ', elsr+el-el_self);
!
! output some quantities
 call write_field(rho,'rho.txt',ascii,size(rho));
 call write_field(grad_el,'grad_el.txt', ascii, size(grad_el))
 call write_field(grad_elsr,'grad_elsr.txt',ascii,size(grad_el))
! call write_field(el_dx,'el_dx.txt', ascii, size(el_dx))
 end program elec
 