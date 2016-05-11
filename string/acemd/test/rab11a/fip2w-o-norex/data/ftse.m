%
 function [epar, eprp]=ftse(i,ibeg,iend,kpar_,kprp_,dpar,dperp,dperp0)

  if (i==ibeg)
   kpar=kpar_;
   dpar0=0
   dperp=dperp / 2 ; % scale to 2D
  elseif (i==iend)
   kpar=kpar_;
   dpar0=1;
   dperp=dperp / 2 ; % scale to 2D
  else
   kpar=kpar_ * 4 ;
   dpar0=0.5;
  end
  epar=0.5*kpar*(dpar-dpar0)^2;
% perpenducular energy
  kprp=kprp_ * 4; % restraint energy is always expressed in units of 2D
  eprp=0.5*kprp*max( dperp - dperp0 / 2, 0 )^2;
