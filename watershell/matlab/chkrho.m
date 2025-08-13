% adapted from ../shape.m with modifications
% adapted from shape3 with mods
if exist('OCTAVE_VERSION'); graphics_toolkit('gnuplot'); end

% type of density :
if ~exist('qsimple'); qsimple=1; end
if ~exist('qmed'); qmed=0; end
if ~exist('qfull'); qfull=0; end
qsimple=1;qmed=0;qfull=0; %
%
qcurv=1 ;% whether to compute gaussian curvature (not needed for cvapx)
qdcurv=1 ;% whether to compute curvature derivatives
%
qdcurv=qdcurv&qcurv ;% otherwise, nonsense
%
% compute various quantities below
rho=zeros(1,nsolv);
%rhop=rho;
%rhopp=rho;
%rhoppp=rho;
% grad (3)
drhox=rho ; % init
drhoy=rho ;
drhoz=rho ; 
% hess (6 by symmetry)
drhoxx=rho ; % init
drhoxy=rho ;
drhoyy=rho ;
drhoxz=rho ;
drhoyz=rho ;
drhozz=rho ;
% third (10 by symmetry);
if (qdcurv)
drhoxxx=rho;
drhoxxy=rho;
drhoxxz=rho;
drhoxyy=rho;
drhoxyz=rho;
drhoxzz=rho;
drhoyyy=rho;
drhoyyz=rho;
drhoyzz=rho;
drhozzz=rho;
end

arad=load('../radii.dat');
padding=3 ;% extra padding of radii (from F code)
arad=arad+padding ;
sig=2 ; % from wshell F output

for i=1:nsurf ;
 if (mod(i,10)==1) ; fprintf(' Computing density for atom %d \n',i) ; end
% NOTE : below we have dr = rwat - rprot ; so grad wrt rwat (needed for forces on the water) does not need (-) sign ; unlike F code ; will double check via FD
 dX=r_solv(1,:)-r_surf(1,i) ; dY=r_solv(2,:)-r_surf(2,i) ; dZ=r_solv(3,:)-r_surf(3,i) ;
 dR2 = dX.^2 + dY.^2 + dZ.^2 ;
 dR = sqrt(dR2) ;
%
 A = (dR+arad(i))/sig/sqrt(2);
 B = (dR-arad(i))/sig/sqrt(2);
 expA=exp(-A.^2);
 expB=exp(-B.^2);
 oosq2pi=1/sqrt(2*pi) ;
 oosig=1/sig;
%
 oodR=1./dR;
% normalize displacement vector
 dX=dX.*oodR;
 dY=dY.*oodR;
 dZ=dZ.*oodR;
%
 if (qsimple)
  rho    = rho + 0.5*( 1 - erf( B )) ;
% NOTE : only rho above is accumulated here, i.e. _not_ rhop, rhopp, rhoppp
% for gradients :
  rhop   = -expB * oosig * oosq2pi ;
  rhopp  =  B.*expB*oosig^2/sqrt(pi);
  if (qdcurv)
   rhoppp =  oosig^3 * oosq2pi * expB.*(1-2*B.^2);
  end
 elseif(qmed)
  rho    = rho + 0.5*(erf( A ) - erf( B )) ;
  rhop   = (expA - expB) * oosig * oosq2pi;
  rhopp  = (B.*expB - A.*expA)*oosig^2/sqrt(pi);
  if(qdcurv)
   rhoppp = oosig^3 * oosq2pi * ( expB.*(1-2*B.^2) - expA.*(1-2*A.^2));
  end
 elseif(qfull)
  rho    = rho + 0.5*(erf( A ) - erf( B )) + (expA - expB)*sig.*oodR*oosq2pi;
  rhop   = -oosq2pi * oodR .* (sig.*oodR .* (expA-expB) + arad(i)*oosig*(expA+expB));
  rhopp  =  oosig*oosq2pi*(oodR.*(1+2*sig^2*oodR.^2 + (arad(i)*oosig)^2) .* (expA-expB) + arad(i)*(2*oodR.^2 + oosig^2).*(expA+expB));
  if(qdcurv)
   rhoppp = -oosig*oosq2pi*( ( 3*oodR.^2 .* (1 + (arad(i)/sig)^2 + 2*(sig.*oodR).^2) + (1+2*(arad(i)/sig)^2)/sig^2 ) .* (expA-expB) + ...
                                      arad(i) * ( ( 6*oodR.^2 + (3+(arad(i)/sig)^2+(dR/sig).^2)/sig^2).*oodR ).* (expA+expB) );
  end
 end
% gradient : NOTE: the derivatives written here are singular at r=0 (could use Taylor expansion) but we will only use them far from zero
 drhox = drhox + rhop .* dX ;
 drhoy = drhoy + rhop .* dY ;
 drhoz = drhoz + rhop .* dZ ;

% Hessian
% diagonal
 d1=( rhopp - rhop .*oodR );
% drhoxx = drhoxx + rhopp .* dX.*dX + rhop .* oodR .* (1-dX.*dX) ;
% drhoyy = drhoyy + rhopp .* dY.*dY + rhop .* oodR .* (1-dY.*dY) ;
% drhozz = drhozz + rhopp .* dZ.*dZ + rhop .* oodR .* (1-dZ.*dZ) ;
% same as above, but fewer ops:
 drhoxx = drhoxx + d1 .* dX.*dX + rhopp - d1 ;
 drhoyy = drhoyy + d1 .* dY.*dY + rhopp - d1 ;
 drhozz = drhozz + d1 .* dZ.*dZ + rhopp - d1 ;
% off-diagonal
 drhoxy = drhoxy + d1 .* dX.*dY ;
 drhoxz = drhoxz + d1 .* dX.*dZ ;
 drhoyz = drhoyz + d1 .* dY.*dZ ;
% Third derivatives
 if(qdcurv)
  d2=d1.*oodR ;
%  drhoxxx = drhoxxx +  dX.^3 .* ( rhoppp - 3*d2 ) + 3 * dX .* d2 ;
%  drhoyyy = drhoyyy +  dY.^3 .* ( rhoppp - 3*d2 ) + 3 * dY .* d2 ;
%  drhozzz = drhozzz +  dZ.^3 .* ( rhoppp - 3*d2 ) + 3 * dZ .* d2 ;
  drhoxxx = drhoxxx +  dX.* (dX.^2 .* ( rhoppp - 3*d2 ) + 3*d2 ) ;
  drhoyyy = drhoyyy +  dY.* (dY.^2 .* ( rhoppp - 3*d2 ) + 3*d2 ) ;
  drhozzz = drhozzz +  dZ.* (dZ.^2 .* ( rhoppp - 3*d2 ) + 3*d2 ) ;
%
  drhoxxy = drhoxxy +  dY.* (dX.^2 .* ( rhoppp - 3*d2 ) + d2 ) ;
  drhoxxz = drhoxxz +  dZ.* (dX.^2 .* ( rhoppp - 3*d2 ) + d2 ) ;
  drhoxyy = drhoxyy +  dX.* (dY.^2 .* ( rhoppp - 3*d2 ) + d2 ) ;
  drhoxyz = drhoxyz +  dX.* dY.*dZ .* ( rhoppp - 3*d2 ) ;
  drhoxzz = drhoxzz +  dX.* (dZ.^2 .* ( rhoppp - 3*d2 ) + d2 ) ;
  drhoyyz = drhoyyz +  dZ.* (dY.^2 .* ( rhoppp - 3*d2 ) + d2 ) ;
  drhoyzz = drhoyzz +  dY.* (dZ.^2 .* ( rhoppp - 3*d2 ) + d2 ) ;
 end
end % natoms


% compare to F output :
% NOTE that the F code will only compute the rho for molecules whhich will experience a force, NOT ALL !
% the max is written in max_rho ; it is 
fprintf(':checking rho against fortran wshell output:') ;
rho_solv=load('../rho_solv.dat');
max_rho=load('../max_rho.dat');
ifr=find(rho<=max_rho) ;
max ( abs ( rho(rho<max_rho)-rho_solv(rho_solv<max_rho) ) ) % this also checks that the sizes match
% also, check computed surface distance :
fprintf(':checking surface distance:') ;
dsurf_solv=load('../dsurf_solv.dat');
dsurf = - ( sqrt(2)*sig*erfinv(1-2*rho) + padding ) ;% this corresponds to "unsafe" inverse in the code b/c low quality appx could put |rho|>1 ; note sign convention
% Need to subtract surface_distance to complete dsurf (i.e. distance to surface, in the negative convention)
surface_distance=load('../surface_distance.dat');
max(abs( (dsurf(ifr)-surface_distance) -dsurf_solv(ifr))) % this will generally be to single prec b/c fortran erfinv code is single prec
corr(dsurf(ifr)',dsurf_solv(ifr)') ; % the correlation is usually very good (O(single) but the absolute error above could be O(-2) for appx

% gradients :
fprintf(':checking ∇rho fortran wshell output:\n') ;
drho_solv=load('../drho_solv.dat');
drho_solv=reshape(drho_solv,3,[]);
fprintf(':∇x:');
max( abs ( drhox(ifr) - drho_solv(1,ifr) ) )
fprintf(':∇y:');
max( abs ( drhoy(ifr) - drho_solv(2,ifr) ) )
fprintf(':∇z:');
max( abs ( drhoz(ifr) - drho_solv(3,ifr) ) )

% hessian :
fprintf(':checking ∇∇rho fortran wshell output:\n') ;
d2rho_solv=load('../d2rho_solv.dat');
d2rho_solv=reshape(d2rho_solv,6,[]);
fprintf(':∇∇xx:');
max( abs ( drhoxx(ifr) - d2rho_solv(1,ifr) ) )
fprintf(':∇∇xy:');
max( abs ( drhoxy(ifr) - d2rho_solv(2,ifr) ) )
fprintf(':∇∇xz:');
max( abs ( drhoxz(ifr) - d2rho_solv(3,ifr) ) )
fprintf(':∇∇yy:');
max( abs ( drhoyy(ifr) - d2rho_solv(4,ifr) ) )
fprintf(':∇∇yz:');
max( abs ( drhoyz(ifr) - d2rho_solv(5,ifr) ) )
fprintf(':∇∇zz:');
max( abs ( drhozz(ifr) - d2rho_solv(6,ifr) ) )
corr(drhoxx(ifr)',d2rho_solv(1,ifr)')

if (qplot)
 figure(1) ; clf ; hold on ; 
% plot(drhoxxx(ifr),d3rho_solv(1,ifr),'k.')
 plot(drhoxx(ifr),'k.-')
 plot(d2rho_solv(1,ifr),'r.-');
 fprintf(':c_P:');
 corr(drhoxx(ifr)',d2rho_solv(1,ifr)')
end

% turd :
fprintf(':checking ∇∇∇rho fortran wshell output:\n') ;
d3rho_solv = load('../d3rho_solv.dat');
d3rho_solv=reshape(d3rho_solv,10,[]);
% NOTE the sign difference here!
fprintf(':∇∇∇xxx:');
max( abs ( drhoxxx(ifr) - d3rho_solv(1,ifr) ) )
fprintf(':∇∇∇xxy:');
max( abs ( drhoxxy(ifr) - d3rho_solv(2,ifr) ) )
fprintf(':∇∇∇xxz:');
max( abs ( drhoxxz(ifr) - d3rho_solv(3,ifr) ) )
fprintf(':∇∇∇xyy:');
max( abs ( drhoxyy(ifr) - d3rho_solv(4,ifr) ) )
fprintf(':∇∇∇xyz:');
max( abs ( drhoxyz(ifr) - d3rho_solv(5,ifr) ) )
fprintf(':∇∇∇xzz:');
max( abs ( drhoxzz(ifr) - d3rho_solv(6,ifr) ) )
fprintf(':∇∇∇yyy:');
max( abs ( drhoyyy(ifr) - d3rho_solv(7,ifr) ) )
fprintf(':∇∇∇yyz:');
max( abs ( drhoyyz(ifr) - d3rho_solv(8,ifr) ) )
fprintf(':∇∇∇yzz:');
max( abs ( drhoyzz(ifr) - d3rho_solv(9,ifr) ) )
fprintf(':∇∇∇zzz:');
max( abs ( drhozzz(ifr) - d3rho_solv(10,ifr) ) )


if (qplot)
% figure(1) ; clf ; hold on ; 
% plot(drhoxxx(ifr),d3rho_solv(1,ifr),'k.')
% plot(drhoxxx(ifr),'k.')
% plot(d3rho_solv(1,ifr),'r.');
 fprintf(':c_P:');
 corr(drhoxxx(ifr)',d3rho_solv(1,ifr)')
end

%
% curvature
if (qcurv)
% below is a single curvature number per point ("total"), but we actualy have two curvature components, both of which are computed in rhosurf
% 12/24 : note that this is the "total" curvature -- not the average ;  to get the average -- halve ; NOTE which curvature is used in MD parametrization
 drhon=sqrt(drhox.^2 + drhoy.^2 + drhoz.^2) ;
 odrhon=1./drhon ;
% NOTE that the formula in rhosurf has an extra 's' at the end of var names; I think that is to differentiate "surface" points from all grid points
 Curv = odrhon .* ( drhoxx + drhoyy + drhozz - odrhon.^2 .* ( drhoxx.*drhox.^2 + drhoyy.*drhoy.^2 + ... 
  2* ( drhoxy.*drhox.*drhoy + ( 0.5 * drhozz.*drhoz + drhoxz.*drhox + drhoyz.*drhoy).*drhoz ) ...
  ) ) / 2; % average curvature
%
 fprintf(':checking curvature :\n') ;
 acurv_solv=load('../acurv_solv.dat');
 fprintf(':C(r):');
 max( abs ( Curv(ifr) - acurv_solv(ifr) ) )
 fprintf(':c_P:');
 corr(Curv(ifr)',acurv_solv(ifr)')
% also check optional curvature correction here:
 iconvex=find(Curv(ifr)<=0); iconvex=ifr(iconvex);
 cCurv(ifr)=Curv(ifr); cCurv(iconvex)=Curv(iconvex)./(1+abs(Curv(iconvex)).* ( dsurf(iconvex)-surface_distance) );
 fprintf(':checking corrected curvature :\n') ;
 max( abs ( cCurv(ifr) - acurv_solv(ifr) ) ) % note that this is limited by the accuracy of the inverse error function in F
 fprintf(':c_P:');
 corr(cCurv(ifr)',acurv_solv(ifr)')

 if (qplot)
  figure(2) ; clf ; hold on ;
  plot(cCurv(ifr),acurv_solv(ifr),'k.')
% plot(drhoxxx(ifr),'k.')
% plot(d3rho_solv(1,ifr),'r.');
 end
%
 if (qdcurv)
% define some aux variables, which maybe (slightly) optimizable later
  L = drhoxx+drhoyy+drhozz ; % laplacian
  HoG1 = drhoxx.*drhox + drhoxy.*drhoy + drhoxz.*drhoz ;  % H . G
  HoG2 = drhoxy.*drhox + drhoyy.*drhoy + drhoyz.*drhoz ;
  HoG3 = drhoxz.*drhox + drhoyz.*drhoy + drhozz.*drhoz ;
% alternative curvature calculation (check)
% note: could also normalize everything by grad magnitude, unless want them available for checking
  GoHoGn = ( drhox.*HoG1 + drhoy.*HoG2 + drhoz.*HoG3 ) .*odrhon.^2  ;
  Curv2 = odrhon.*(L - GoHoGn) /2;
  assert(max(abs(Curv(:) - Curv2(:))<1e-12)) % make sure they agree befire proceeeding
% first term (divergence of Hessian):
  divH1 = drhoxxx+drhoxyy+drhoxzz; % to be scaled by 1/g
  divH2 = drhoxxy+drhoyyy+drhoyzz;
  divH3 = drhoxxz+drhoyyz+drhozzz;
% second term : (could do L <= 0.5 * (L - 3*GoHoGn) ; then take 2 outside
% T21 = ( L-3*GoHoGn + 2 * drhoxx ).* HoG1 + ( L-3*GoHoGn + 2 * drhoxy ).* HoG2 + ( L-3*GoHoGn + 2 * drhoxz ).* HoG3 ; % to be scaled by 1/g^3
% T22 = ( L-3*GoHoGn + 2 * drhoxy ).* HoG1 + ( L-3*GoHoGn + 2 * drhoyy ).* HoG2 + ( L-3*GoHoGn + 2 * drhoyz ).* HoG3 ;
% T23 = ( L-3*GoHoGn + 2 * drhoxz ).* HoG1 + ( L-3*GoHoGn + 2 * drhoyz ).* HoG2 + ( L-3*GoHoGn + 2 * drhozz ).* HoG3 ;
  T21 = ( L-3*GoHoGn + 2 * drhoxx ).* HoG1 + (              2 * drhoxy ).* HoG2 + (              2 * drhoxz ).* HoG3 ; % to be scaled by 1/g^3
  T22 = (              2 * drhoxy ).* HoG1 + ( L-3*GoHoGn + 2 * drhoyy ).* HoG2 + (              2 * drhoyz ).* HoG3 ;
  T23 = (              2 * drhoxz ).* HoG1 + (              2 * drhoyz ).* HoG2 + ( L-3*GoHoGn + 2 * drhozz ).* HoG3 ;
% third term (reuse HoG)
% hess grad wrt x:
  HoG1 = drhoxxx.*drhox + drhoxxy.*drhoy + drhoxxz.*drhoz ;
  HoG2 = drhoxxy.*drhox + drhoxyy.*drhoy + drhoxyz.*drhoz ;
  HoG3 = drhoxxz.*drhox + drhoxyz.*drhoy + drhoxzz.*drhoz ;
  T31 = drhox.*HoG1 + drhoy.*HoG2 + drhoz.*HoG3;
% grad wrt y:
  HoG1 = drhoxxy.*drhox + drhoxyy.*drhoy + drhoxyz.*drhoz ;
  HoG2 = drhoxyy.*drhox + drhoyyy.*drhoy + drhoyyz.*drhoz ;
  HoG3 = drhoxyz.*drhox + drhoyyz.*drhoy + drhoyzz.*drhoz ;
  T32 = drhox.*HoG1 + drhoy.*HoG2 + drhoz.*HoG3;
%
  HoG1 = drhoxxz.*drhox + drhoxyz.*drhoy + drhoxzz.*drhoz ;
  HoG2 = drhoxyz.*drhox + drhoyyz.*drhoy + drhoyzz.*drhoz ;
  HoG3 = drhoxzz.*drhox + drhoyzz.*drhoy + drhozzz.*drhoz ;
  T33 = drhox.*HoG1 + drhoy.*HoG2 + drhoz.*HoG3;
%
  dCurvx=odrhon .* ( divH1 - (T21+T31).*odrhon.^2 )/2;
  dCurvy=odrhon .* ( divH2 - (T22+T32).*odrhon.^2 )/2;
  dCurvz=odrhon .* ( divH3 - (T23+T33).*odrhon.^2 )/2;
%
% debugging hacks :
% Curv=L; %pass
% dCurvx=divH1;
% dCurvy=divH2;
% dCurvz=divH3;
% ============
% Curv=L.*odrhon ; % fails w/ rmse ; ok with median error -- unclear why, in light of what's below
% dCurvx = (divH1 - L.*(drhoxx.*drhox + drhoxy.*drhoy + drhoxz.*drhoz).*odrhon.^2).*odrhon ;
% dCurvy = (divH2 - L.*(drhoxy.*drhox + drhoyy.*drhoy + drhoyz.*drhoz).*odrhon.^2).*odrhon ;
% dCurvz = (divH3 - L.*(drhoxz.*drhox + drhoyz.*drhoy + drhozz.*drhoz).*odrhon.^2).*odrhon ;

% Curv=odrhon ; % ?  the FD rmse error is huge but is does scale as O2 ...
% dCurvx = -(drhoxx.*drhox + drhoxy.*drhoy + drhoxz.*drhoz).*odrhon.^3 ;
% dCurvy = -(drhoxy.*drhox + drhoyy.*drhoy + drhoyz.*drhoz).*odrhon.^3 ;
% dCurvz = -(drhoxz.*drhox + drhoyz.*drhoy + drhozz.*drhoz).*odrhon.^3 ;

% Curv=drhon ; % ok w/ rmse
% dCurvx = (drhoxx.*drhox + drhoxy.*drhoy + drhoxz.*drhoz).*odrhon ;
% dCurvy = (drhoxy.*drhox + drhoyy.*drhoy + drhoyz.*drhoz).*odrhon ;
% dCurvz = (drhoxz.*drhox + drhoyz.*drhoy + drhozz.*drhoz).*odrhon ;

% Curv=sqrt(drhox.^2 + drhoy.^2 + drhoz.^2) ; % fails (where curv is near zero)
% dCurvx = (drhoxx.*drhox + drhoxy.*drhoy + drhoxz.*drhoz)./Curv ;
% dCurvy = (drhoxy.*drhox + drhoyy.*drhoy + drhoyz.*drhoz)./Curv ;
% dCurvz = (drhoxz.*drhox + drhoyz.*drhoy + drhozz.*drhoz)./Curv ;

% Curv=0.5*(drhox.^2 + drhoy.^2 + drhoz.^2) ; % works
% Curv=0.5*(drhon.^2) ; % works
% dCurvx = drhoxx.*drhox + drhoxy.*drhoy + drhoxz.*drhoz;
% dCurvy = drhoxy.*drhox + drhoyy.*drhoy + drhoyz.*drhoz;
% dCurvz = drhoxz.*drhox + drhoyz.*drhoy + drhozz.*drhoz;

% Curv=1./(drhon.^2) ; % fails ; looks like the gradient is simply too small in some regions
% dCurvx = -2*(drhoxx.*drhox + drhoxy.*drhoy + drhoxz.*drhoz).*Curv.^2;
% dCurvy = -2*(drhoxy.*drhox + drhoyy.*drhoy + drhoyz.*drhoz).*Curv.^2;
% dCurvz = -2*(drhoxz.*drhox + drhoyz.*drhoy + drhozz.*drhoz).*Curv.^2;
%
  fprintf(':checking curvature gradient:\n') ; % to output in code !
  dacurv_solv=load('../dacurv_solv.dat');
  dacurv_solv=reshape(dacurv_solv,3,[]);
  fprintf(':∇x:');
  max( abs ( dCurvx(ifr) - dacurv_solv(1,ifr) ) )
  fprintf(':∇y:');
  max( abs ( dCurvy(ifr) - dacurv_solv(2,ifr) ) )
  fprintf(':∇z:');
  max( abs ( dCurvz(ifr) - dacurv_solv(3,ifr) ) )
  corr(dCurvx(ifr)',dacurv_solv(1,ifr)')
  if (qplot)
   figure(3) ; clf ; hold on ;
   plot(dCurvy(ifr),dacurv_solv(2,ifr),'k.')
% plot(drhoxxx(ifr),'k.')
% plot(d3rho_solv(1,ifr),'r.');
  end
%
 end
end

