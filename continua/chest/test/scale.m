% 9.25.17 : comparison of multigrid execution times for several grids ; CPU and GPU

nx=[98 128 256 512];
ny=[98 128 256 256];
nz=[98 128 256 256];
gput=[0.042 0.071 0.368 0.706];
cput=[1.015 2.389 19.0 37.67];

%98x98x98	0.042		1.015		4 levels
%128x128x128	0.071		2.389		6 levels
%256x256x256	0.368		19.0		7 levels
%512x256x256	0.706		37.67		7 levels

npt=nx.*ny.*nz;

cog = cput./gput;

close all;
h=figure(1);

hold on;
plot(npt,cput,'k-x') ; hold on;
plot(npt, gput,'r-x') ;

% plot a linear law
plot(npt, 5*npt/npt(end), 'k--') ;

box on ; grid on;

legend('AMD FX-8350','GTX 1080','Linear law',2);

set(gca, 'xscale','log') ;
set(gca, 'yscale','log') ;

xlabel('\it DOFs', 'Fontsize',14)
ylabel('\it Time(s)', 'Fontsize',14);

print(gcf, '-depsc2', 'poisson.eps');


