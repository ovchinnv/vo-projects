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

ufile='eps.dat';

if (0)
f=fopen(uexactfile,'r');
rl=fread(f,1,'int32')
dim=fread(f,3,'int32')
rl=fread(f,2,'int32')
uex=fread(f,nx*ny*nz,rlfmt);
uex=reshape(uex,nx,ny,nz);
fclose(f);
end

f=fopen(ufile,'r');
rl=fread(f,1,'int32')
dim=fread(f,3,'int32')
rl=fread(f,2,'int32')
uch=fread(f,nx*ny*nz,rlfmt);
uch=reshape(uch,nx,ny,nz);

%d=uch-uex;
d=uch;
% show slice
kz=floor(nz/2)-1
%kz=20;
off=1;
%
%pcolor(x(1+off:nx-off),y(1+off:ny-off),squeeze(d(1+off:nx-off,1+off:ny-off,kz))); shading interp; colorbar; box on
contourslice(x,y,z,d,0,0,0)

colormap(autumn)


