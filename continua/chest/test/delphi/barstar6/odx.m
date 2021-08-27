% load a hacked openDX file to see what the eps values are

e=load('eps.dat');

nx=98;
ny=98;
nz=98;

size(e);
e=reshape(e,nz,ny,nx);

kx=floor(nx/2)-1
off=0;
pcolor(squeeze(e(1+off:nz-off,1+off:ny-off,kx))); shading interp; colorbar; box on
