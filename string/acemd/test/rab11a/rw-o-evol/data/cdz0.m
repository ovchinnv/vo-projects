% plot FE computed in CHARMM
% it is the same!
close all;

rmsa=load('rmsd.dat');
[niter, nrep]=size(rmsa);

figure; hold on;box on;

shift=0.1;
subplot(1,2,1);hold on

for i=2:nrep
 plot(rmsa(:,i)+i*shift,'-');
end
box on;

ylabel('Replica index');
xlabel('iteration');
%%%%%%%%%%%%%%%%%%%%%%% average
subplot(1,2,2);hold on
plot( sqrt(mean(rmsa(:,2:end).^2,2))/(nrep-1) );
box on;
ylabel('Replica index');
xlabel('iteration');

set(gcf, 'paperpositionmode', 'auto');

%print(gcf, '-depsc2', 'dz0.eps');
%print(gcf, '-djpeg100', 'dz0.jpg');
