% visualize exact and computed solutions
clf;

rlfmt='float32';

gridfile='xyz_test.xyz';

f=fopen(gridfile,'r');

rl=fread(f,1,'int32')
nx=fread(f,1,'int32')
ny=fread(f,1,'int32')
nz=fread(f,1,'int32')
rl=fread(f,2,'int32')

x=fread(f,nx,rlfmt);
y=fread(f,ny,rlfmt);
z=fread(f,nz,rlfmt);

fclose(f);

uexactfile='uexact_test.dat';
ufile='solution.dat';
%ufile='rhs_test.dat';

f=fopen(uexactfile,'r');
rl=fread(f,1,'int32')
dim=fread(f,3,'int32')
rl=fread(f,2,'int32')
uex=fread(f,nx*ny*nz,rlfmt);
uex=reshape(uex,nx,ny,nz);
fclose(f);
f=fopen(ufile,'r');
rl=fread(f,1,'int32')
dim=fread(f,3,'int32')
rl=fread(f,2,'int32')
uch=fread(f,nx*ny*nz,rlfmt);
uch=reshape(uch,nx,ny,nz);

d=uch-uex;
% show slice
kz=floor(nz/2)-1
%kz=20;
off=1;
%
%pcolor(x,y,uex(2:nx-1,2:ny-1,kz)); shading flat; colorbar; box on
%mesh(x(1+off:nx-off),y(1+off:ny-off),squeeze(d(1+off:nx-off,1+off:ny-off,kz))); shading flat; colorbar; box on
pcolor(x(1+off:nx-off),y(1+off:ny-off),squeeze(d(1+off:nx-off,1+off:ny-off,kz))); shading flat; colorbar; box on
%pcolor(y(1+off:nx-off),z(1+off:ny-off),squeeze(d(kz,1+off:nx-off,1+off:ny-off))); shading flat; colorbar; box on
%pcolor(y(1+off:nx-off),z(1+off:ny-off),squeeze(d(1+off:nx-off,kz,1+off:ny-off))); shading flat; colorbar; box on

%mesh(x(1+off:nx-off),y(1+off:ny-off),uch(1+off:nx-off,1+off:ny-off,kz)); shading flat; colorbar; box on
%mesh(x(2:nx-1),y(2:ny-1),d(2:nx-1,2:ny-1,kz)); shading flat; colorbar; box on
%hold on ;
%caxis([-50 50])


