% compute free energy from forces acting on the CV
% assuming 1D
%
close all;

cv     =load('deca_tamd_cv.dat');
fcv    =load('deca_tamd_cv_force.dat');
minlen =min(length(cv), length(fcv));
data   = [ cv(1:minlen), fcv(1:minlen) ] ;


x0=12 ;
x1=32 ;
nxcor=50 ;

xcor=linspace(x0, x1, nxcor);
xcen=xcor(2:end)-xcor(1:end-1);
nxcen=nxcor-1;
dx=xcor(2)-xcor(1);

% bin cvs and gradients
% determine bin number
data(:,3)= floor ( ( data(:,1)-x0 )/ dx ) + 1;
data=sortrows(data,3);
% compute average force in each bin
fx =zeros(1,nxcen);
fxe=zeros(1,nxcen);
for ibin=1:nxcen
 inds=(find(data(:,3)==ibin));
 dd=data(inds,2);
 fx(ibin)=mean(dd);
 if (~isnan(fx(ibin)))
  fxe(ibin)=std(dd);
 end
end

% integrate derivative :
fe=cumsum(fx)*dx ; fe=fe-min(fe);
plot(xcor(1:nxcen), fe,'k*-');
%errorbar(xcor(1:nxcen), fe, fxe);

