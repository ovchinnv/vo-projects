% load spline file:
addpath('~/scripts/matlab');
splinefile='../sbf-md-rcut8spl2dkz.dat';
spl2d=readspl2d(splinefile);
%
qcrvsmooth=0;
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
 plot(frcx,'k.-')
 plot(forces_solv(1,ifr),'r.');
 fprintf(':c_P:');
end

% compute surface force (but without curvature gradient) :
fsurf=zeros(1,nsolv);
for i=iforced
% need to limit curvature to grid support (or to manual spec from wshell input)
 if (qcrvsmooth)
 else % cap
  ka=min(kkgrid(end),max(kkgrid(1),Curv(i)));
 end
 fsurf(i)=fspl2d(spl2d,ka,dsurf(i)-surface_distance);
end
frcx = -fsurf(ifr) .* ddsurf_solv(ifr) .* drhox(ifr);
frcy = -fsurf(ifr) .* ddsurf_solv(ifr) .* drhoy(ifr);
frcz = -fsurf(ifr) .* ddsurf_solv(ifr) .* drhoz(ifr);

fprintf(':checking surface force (no curv. gradient):\n');
forces_solv=load('../forces_solv.dat');
forces_solv=reshape(forces_solv,3,[]);
fprintf(':fx:');
max( abs ( frcx - forces_solv(1,ifr) ) )
fprintf(':fy:');
max( abs ( frcy - forces_solv(2,ifr) ) )
fprintf(':fz:');
max( abs ( frcz - forces_solv(3,ifr) ) )
corr(frcx',forces_solv(1,ifr)')




