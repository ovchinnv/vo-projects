
if exist('OCTAVE_VERSION','builtin')
 clear graphics_toolkit;
 graphics_toolkit('gnuplot')
% graphics_toolkit('fltk')
end

close;

map=load('adia.dat');

phi=map(:,1);
psi=map(:,2);
ene=map(:,3);

n=sqrt(length(phi));
phi=reshape(phi,n,n);
psi=reshape(psi,n,n);
ene=reshape(ene,n,n);

pcolor(phi, psi, ene);shading interp;%colorbar;
hold on;colormap pink

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now can overlay paths
name='dihe_dcd16.dat';

%
d=load(name);
%
skip=500; % number of data-points in a single trajectory (THIS HAS TO BE CORRECT)
%
% reorganize data into a matrix
id=d(:,1);
phi=d(:,2);
psi=d(:,3);
%
phi=reshape(phi,skip,[]);
psi=reshape(psi,skip,[]);
id=reshape(id,skip,[]);
%
styles={'r','g','b','m','y','w'};
%
%plot(phi,psi,char(styles(1)))
%
[m,nrep]=size(phi)

reps=[1:6, 10:16]
reps=[2:2:16]
reps=[1:8]
reps=[1:16]
%
step=1;
j=0;
for k=1:length(reps)
 i=reps(k);
 j=j+1;
 plot(phi(1:step:end,i), psi(1:step:end,i),[char(styles( mod(j-1,length(styles))+1 )),'.'],'MarkerSize',5);
end

end % traj

if (0)
% overlay (initial) string
name='../data/diala22_zts_.dihe'
d=load(name);
id=d(:,1);
phi=d(:,2);
psi=d(:,3);

% average 
%phi=reshape(phi,nrep,[]); phi=mean(phi,2);
%psi=reshape(psi,nrep,[]); psi=mean(psi,2);

plot(phi, psi,'k.','MarkerSize',25); hold on

end

set(gca,'FontSize',16);

xlabel('\phi', 'FontSize',16);
ylabel('\psi', 'FontSize',16);

set(gcf, 'paperpositionmode', 'auto');
%print(gcf, '-dpdf', 'surf1.pdf');
%print(gcf, '-depsc2', 'surf1.eps');
print(gcf, '-dtiff', 'surf1.tif');

