%
close all;

if (~exist('read'))
 read=1;
end

read=1

if (read==1)
%%%%%%%%%% load multiple files %%%%%%%%%%%%
fnames={'../data/rmsd.dat'};
%
clear rmsa;
for i=1:length(fnames)
 fname=char(fnames(i));
 if (i==1)
  rmsa=load(fname);
 else
  rmsa=[rmsa; load(fname)];
 end 
end
read=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[niter, nrep]=size(rmsa);

figure('position', [100 100 1000 400]); hold on;box on;

shift=0.2;
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
title('RMSD from initial string');
set(gcf, 'paperpositionmode', 'auto');

print(gcf, '-depsc2', 'fts2-rmsd0.eps');
print(gcf, '-djpeg100', 'fts2-rmsd0.jpg');
