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

ufile='eps.dat';

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
clf
slice(x,y,z,d,[],[],0) ; shading interp ; colorbar % not a clear plot ;
if (1)
p=patch(isosurface(x,y,z,d,5));
p.FaceLighting='gouraud'
p.FaceColor=[1 1 1]*0.75
p.EdgeColor='none'
end
%axis manual
%p=patch(isosurface(x,y,z,d,5)) ; hold on
%p.FaceColor=''
%p.EdgeColor='none'

lighting gouraud
camlight

colormap(jet)
axis equal;
view(3)

camlight

set(gca, 'fontsize',12)
xlabel('\it X(\AA)','interpreter','latex');
ylabel('\it Y(\AA)','interpreter','latex');
zlabel('\it Z(\AA)','interpreter','latex');

set(gcf, 'paperpositionmode','auto')
%print(gcf, '-depsc2', '-painters', 'eps_map.eps'); % errors in file
%print(gcf, '-depsc2', '-opengl', 'eps_map.eps'); % low quality
print(gcf, '-djpeg100', '-opengl', 'eps_map.jpg'); % acceptable quality

p.Visible='off'
view(2)
box on;
set(gca, 'TickDir','out')
print(gcf, '-djpeg100', '-opengl', 'eps_map-noprot.jpg'); % acceptable quality
