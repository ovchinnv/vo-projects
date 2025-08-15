% load spline file:
addpath('~/scripts/matlab');
splinefile='../sbf-md-rcut8spl2dkz.dat';
spl2d=readspl2d(splinefile);
%
qcrvsmooth=1; % smooth by Gaussian convolution
qcurvgrad=1; % whether to include forces from the curvature gradient
%
% show a simple force plot at 0 curvature :
zzgrid=spl2d.breaks{2};
kkgrid=spl2d.breaks{1};
if (qplot)
 zz=zzgrid;
 kk=0;
 f0=fspl2d(spl2d,kk,zz); 
 figure(4); clf ; hold on ;
 plot(zz,f0,'k*-') ; % note that the plot is (-) outside (z<0) ; (+) inside wshell boundary ; need to flip sign for fc, bc normal is inward [recall, +grad of density]
end
%
% first, evaluate the 0-curvature, surface force ; decide whether to use derivative of drho_inverse (F code parameter), which can be very large
ddsurf_solv(ifr) = sig * sqrt(2*pi) * exp ( (dsurf(ifr) + padding).^2 / (2*sig^2)  ) ;# note sign convention
iforced=find ( (dsurf(ifr)-surface_distance) <= zzgrid(end) ) ; iforced=ifr(iforced);

fsurf0=zeros(1,nsolv);
fsurf0(iforced) = fspl2d(spl2d,0,dsurf(iforced)-surface_distance) ;

frcx0 = -fsurf0(ifr) .* ddsurf_solv(ifr) .* drhox(ifr);
frcy0 = -fsurf0(ifr) .* ddsurf_solv(ifr) .* drhoy(ifr);
frcz0 = -fsurf0(ifr) .* ddsurf_solv(ifr) .* drhoz(ifr);

% compare with F
fprintf(':checking 0-curvature surface force:\n');
forces_solv=load('../forces_solv.dat');
forces_solv=reshape(forces_solv,3,[]);
fprintf(':fx:');
max( abs ( frcx0 - forces_solv(1,ifr) ) )
fprintf(':fy:');
max( abs ( frcy0 - forces_solv(2,ifr) ) )
fprintf(':fz:');
max( abs ( frcz0 - forces_solv(3,ifr) ) )
corr(frcx0',forces_solv(1,ifr)')

if (qplot)
 figure(1) ; clf ; hold on ; 
% plot(drhoxxx(ifr),d3rho_solv(1,ifr),'k.')
 plot(frcx0,'k.-')
 plot(forces_solv(1,ifr),'r.');
 fprintf(':c_P:');
end

% compute surface force (but without curvature gradient) :
fsurf=zeros(1,nsolv);
fpot=zeros(1,nsolv);
fcurv=zeros(1,nsolv);
kmin=kkgrid(1);
kmax=kkgrid(end);
%kmin=-0.1;
%kmax= 0.1;
ksig=0.03; % NOTE different sigma from the one for rho !

if qcurvgrad % define new splines for evaluating derivative wrt curvature
 spl2dintd=splint2dy(spl2d);
 spl2did_dk=spldiff2dx(spl2dintd);
end

for i=iforced
% need to limit curvature to grid support (or to manual spec from wshell input)
 if (qcrvsmooth)
  k1 = (kmin-Curv(i))/(ksig*sqrt(2));
  k2 = (kmax-Curv(i))/(ksig*sqrt(2));
  ka = 0.5 * ( kmin+kmax + ksig*sqrt(2) * ( k1*erf(k1) - k2*erf(k2) ) ) + ksig*oosq2pi*(exp(-k1^2)-exp(-k2^2));
  kadiff=0.5*(erf(k2)-erf(k1)); % curvature correction
 else % cap
  ka=min(kmax,max(kmin,Curv(i)));
  kadiff=1;
 end
 fsurf(i)=fspl2d(spl2d,ka,dsurf(i)-surface_distance); % custom function
% fsurf(i)=ppval(spl2d,[ka,dsurf(i)-surface_distance]'); % native -- does not work
 if (qcurvgrad)
  fpot(i) =             fspl2d(spl2dintd,ka,dsurf(i)-surface_distance)-fspl2d(spl2dintd,ka,zzgrid(end)); % 2D potential of mean force
  fcurv(i) = kadiff * ( fspl2d(spl2did_dk,ka,dsurf(i)-surface_distance)-fspl2d(spl2did_dk,ka,zzgrid(end)) ) ; % derivative of potential wrt curvature
 end
end
% compare splined prefactors:
fprintf(':checking spline prefactors:\n');
sbfe_solv=load('../sbfe_solv.dat');
sbfe_solv=reshape(sbfe_solv,[],nsolv);
fprintf(':distance force:');
max( abs (sbfe_solv(1,ifr) + sbfe_solv(2,ifr) - fsurf(ifr) ) )
corr(sbfe_solv(1,ifr)' + sbfe_solv(2,ifr)',fsurf(ifr)')
%
if qcurvgrad
 fprintf(':2D PMF:');
 max( abs (sbfe_solv(3,ifr) + sbfe_solv(4,ifr) - fpot(ifr) ) )
 corr(sbfe_solv(3,ifr)' + sbfe_solv(4,ifr)',fpot(ifr)') % note : for this to match have to set outputfreq=1 in wshell, otherwise, will get 0 for potential
% sbfe_pot=sbfe_solv(3,:)+sbfe_solv(4,:) ;
% plot(sbfe_pot)
%
 fprintf(':curvature force:');
 max( abs (sbfe_solv(5,ifr) - fcurv(ifr) ) )
 corr(sbfe_solv(5,ifr)',fcurv(ifr)')
end

% note the composition of the surface force :
% first part is from the spline, which has the "sign" with which (-) z-positions give (-) force [ so we reverse it ]; the grad is inward pointing, and the
% distance derivative does not change the sign because it increases monotonically with rho in view of the (-) sign convention 
frcx = -fsurf(ifr) .* ddsurf_solv(ifr) .* drhox(ifr);
frcy = -fsurf(ifr) .* ddsurf_solv(ifr) .* drhoy(ifr);
frcz = -fsurf(ifr) .* ddsurf_solv(ifr) .* drhoz(ifr);

% curvature force :
if qcurvgrad
 fcvx =  fcurv(ifr) .* dCurvx(ifr);
 fcvy =  fcurv(ifr) .* dCurvy(ifr);
 fcvz =  fcurv(ifr) .* dCurvz(ifr);
%
 frcx=frcx-fcvx;
 frcy=frcy-fcvy;
 frcz=frcz-fcvz;
end

if (qcurvgrad)
 fprintf(':checking surface force (with curv. gradient):\n');
else
 fprintf(':checking surface force (no curv. gradient):\n');
end

forces_solv=load('../forces_solv.dat');
forces_solv=reshape(forces_solv,3,[]);
fprintf(':fx:');
max( abs ( frcx - forces_solv(1,ifr) ) )
fprintf(':fy:');
max( abs ( frcy - forces_solv(2,ifr) ) )
fprintf(':fz:');
max( abs ( frcz - forces_solv(3,ifr) ) )
corr(frcx',forces_solv(1,ifr)')

