
s=load('arcl.dat');
[nsamples,nrep]=size(s) ;
plot(1:nsamples, sum(s(:,2:end),2),'kx-');

box on;
ylabel('Replica index');
xlabel('iteration');

set(gcf, 'paperpositionmode', 'auto');

%print(gcf, '-depsc2', 'slen.eps');
%print(gcf, '-djpeg100', 'slen.jpg');



