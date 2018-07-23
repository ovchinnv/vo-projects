% compute average value from _steady-state_ PDF given by Fokker-Planck equation for Brownian motion
% in 1D with reflective boundary conditions, assuming a linear potential (below)
% The FPE is given by
% .
% P = D [γP + P']' ; γ = α / kBT ; assuming the potential has the form U = αx + C
% the domain is I=(0..1) and the boundary conditions are P'(0,t) = P'(1,t) = 0
% for the initial condition, take δ(x-x0) where x0 is in I.
%
function [x,dx]=xave_s(g)

% g    : drift velocity (gamma above)
tol=1e-8 ; % tolerance at which to switch to taylor expansion of steady-state contribution to average

if ( abs(g) > tol )
 og=1/g; eg=exp(g);
 x = og - 1/(eg-1);
 if (nargout>1)
  dx=-og^2 + eg/(eg-1)^2 ;
 end
else
 x = 0.5 * (1-g/6);
 if (nargout>1)
  dx=-1/12;
 end
end
