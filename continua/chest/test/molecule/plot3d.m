% read data in plot3D ascii format
% close all


q3d=true;
ascii=true;

gridname='xyz.xyz';
%dataname='eps.dat';
dataname='eps.dat';

if ~exist('read')
 read=1;
end

if (read)
 if (ascii)
  grd=textread(gridname, '%f','bufsize', 1000000);
  dat=textread(dataname, '%f','bufsize', 1000000);
 end
 read=0;
end 
%
nx=floor(grd(1));
ny=floor(grd(2));
%
nx2=floor(dat(1));
ny2=floor(dat(2));
% other operations later
if (q3d)
 nz=floor(grd(3));
 nz2=floor(dat(3));
 i=4; i2=8;
else
 nz=1;
 nz2=1; 
 i=4; i2=3;
end
%
if ( (nx~=nx2) || (ny~=ny2) || (nz~=nz2) )
 error ( 'DATA SIZE DOES NOT MATCH GRID SIZE. ABORT.' );
 return
end
%
%
if (q3d)
% grid
% in case the grid is full as in plot3d
 x=grd(i:nx+i-1);
 y=grd(nx*ny*nz+i:nx:nx*ny*nz+nx*ny+i-1);
 z=grd(2*nx*ny*nz+i:ny*ny:end);
%
 d=reshape(dat (i2:nx*ny*nz+i2-1), nx,ny,nz);
else
% grid
% x=reshape(grd(i:nx*ny+i-1),  nx,ny); xx=x(:,1); % in case the grid is full
% y=reshape(grd(nx*ny+i:end),nx,ny); yy=y(1,:);
 x=grd(i:nx+i-1);
 y=grd(nx+i:nx+i+ny-1);
%
 d=reshape(dat(i2:nx*ny+i2-1), nx,ny);
end

% show slice
kz=round(nz/2);

%
%mesh(x,y,d(:,:,kz)); shading flat; colorbar; box on
%mesh(x(2:nx-1),y(2:ny-1),d(2:end-1,2:end-1,kz)); shading flat; colorbar; box on
%hold on ;
%pcolor(x,y,d(:,:,kz)); shading interp; colorbar; box on
%meshc(x,y,d); shading flat; colorbar; box on

% plot isosurface :

isosurface(x,y,z,d,5)
axis equal;

return

jy=10;
mesh(x,z,squeeze(d(:,jy,:))); shading flat; colorbar; box on
%mesh(x,z,squeeze(g(:,jy,:))); shading flat; colorbar; box on
%mesh(x,z,squeeze(dg(:,jy,:))); shading flat; colorbar; box on
mesh(x(2:end-1),z(2:end-1),squeeze(d(2:end-1,jy,2:end-1))); shading flat; colorbar; box on

