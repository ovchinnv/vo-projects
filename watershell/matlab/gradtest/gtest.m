% compute filtered bump and its gradient in 2D to test (see notes ; consistent w/ code)
% NOTE : this is being done in 3D, because the filtering formula given here is incorrect for 2D
% however, the derivative test should still be the same !
%
% define grid

xc=0;
yc=0;

dx=1;
dy=1;

x0=-dx ;
x1= dx ;
y0=-dy ;
y1= dy ;

npx=200;
npy=202;

s=0.11 ; %filter width
radpad=0.5 ; % radius + padding

xx=linspace(x0,x1,npx);
yy=linspace(y0,y1,npy);

[X,Y]=meshgrid(xx,yy) ;

dX = X-xc ;
dY = Y-yc ;

dR = sqrt(dX.^2 + dY.^2) ;

c=dR/s;

a=1/sqrt(2) * ( c + radpad/s ) ;
b=1/sqrt(2) * ( c - radpad/s ) ; 
expa=exp(-a.^2);
expb=exp(-b.^2);

G=1/2 * (erf(a) - erf(b)) + 1/sqrt(2*pi) * ( expa - expb )./c ;
%pcolor(X,Y,G); colorbar ; shading interp ; 
%meshc(X,Y,G); colorbar ; shading interp ; 


dGdX_code = - 1/sqrt(2*pi) * ( (expa + expb) * radpad/s + 0.5 * (expa - expb)./c)./(c*s).^2 .*dX ;% from the watershell code
dGdX_fix =  - 1/sqrt(2*pi) * ( (expa + expb) * radpad/s +       (expa - expb)./c)./(c*s).^2 .*dX ;% corrected analytical expression

dGdX_fd = zeros(size(G));
dGdX_fd(:,2:end-1) = ( G(:,3:end) - G(:,1:end-2) )/(xx(3)-xx(1));

dGdXerr = dGdX_fd - dGdX_fix; % CORRECT !!!
%dGdXerr = dGdX_fd - dGdX_code; % WRONG !!!!

pcolor( X(:,2:end-1) ,Y(:,2:end-1), dGdXerr(:,2:end-1) ); colorbar ; shading interp ; 

%pcolor(X,Y,dGdX_code); colorbar ; shading interp ; 
%pcolor(X,Y,dGdX_fix); colorbar ; shading interp ; 
%pcolor(X,Y,dGdX_fd); colorbar ; shading interp ; 
%pcolor(X,Y,dGdX_fix-dGdX_code); colorbar ; shading interp ; 



