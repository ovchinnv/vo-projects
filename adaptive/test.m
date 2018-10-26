% plot temperature history
%
graphics_toolkit('gnuplot');

close all;

files={ 'adaptive-test1.dat'};
files={ 'adaptive-test2.dat'};
files={ 'adaptive-test3.dat'};
%files={ 'adaptive-test4.dat'};


for i=1:length(files)
 file=char(files(i));
 if (i==1)
  d=load(file);
 else
  d=[d; load(file)];
 end 
end

time=d(:,1);
param=d(:,2);

plot(time, param,'k-','linewidth',0.5) ; box on; hold on;
xlabel('\it t(ns)', 'Fontsize',14) ;
ylabel('\it k_{fc}(kcal/mol/A^2)', 'Fontsize',14);

set(gca, 'yscale','log')
%ylim([250 500]);

set(gcf,'position',[100 100 800 200]);

set(gcf, 'paperpositionmode','auto') ;
print(gcf, '-depsc2', 'kfc_hist');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid :
npt=30;
[pdf,bins]=hist(param,npt); pdf=pdf/length(time);
figure ; hold on; grid on; box on
plot(bins, pdf,'ko-')
% superpose powerlaw
expo=-1;
plot(bins,(bins.^expo) *pdf(end/2) / bins(end/2)^expo ,'r--', 'linewidth',2)

set(gca, 'xscale','log')
set(gca, 'yscale','log')

xlabel('\it k_{fc}(kcal/mol/A^2)', 'fontsize',14)
ylabel('\it PDF', 'fontsize',14);
legend('pdf','power-law');
set(gcf, 'paperpositionmode','auto') ;
%print(gcf, '-depsc2', 'pdf');
