% read data in plot3D ascii format
% close all


q3d=false;
ascii=true;

gridname='xy_test.xyz';
dataname='solution.dat';

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
 i=4; i2=4;
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
else
% grid
% x=reshape(grd(i:nx*ny+i-1),  nx,ny); xx=x(:,1); % in case the grid is full
% y=reshape(grd(nx*ny+i:end),nx,ny); yy=y(1,:);
 x=grd(i:nx+i-1);
 y=grd(nx+i:nx+i+ny-1);
% 
 d=reshape(dat(i2:nx*ny+i2-1), nx,ny);
end

%pcolor(x,y,d); shading flat; colorbar; box on
meshc(x,y,d); shading flat; colorbar; box on


