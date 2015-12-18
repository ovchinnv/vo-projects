% diagonalize covariance matrix
%

delf='int32';
intf='int32';
ff='double';

fname='rab11a_cov.dat';

fid=fopen(fname,'r');

d=fread(fid,1,delf);

n=fread(fid,2,intf);
nx=n(1); ny=n(2);

d=fread(fid,2,delf);

a=fread(fid,[nx ny],ff); pcolor(a); shading interp; colorbar;

% 
%mass=12.0110;
%a=a/mass/mass;

%save 'a.dat' -ascii a

%return

%diagonalize
[ev,e]=eig(a);
e=diag(e)