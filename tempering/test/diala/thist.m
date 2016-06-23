% plot temperature history
%
close all;
kboltz=1.98e-3 ;

file='tempering.series.txt';

d=load(file);

step=d(:,1);
temp=d(:,2);
ener=d(:,3);

dt = 10 * 1 / 1000000 ; % plugin_freq x timestep (fs) x (ns / fs) => to convert to time (ns)

plot(step*dt, temp,'k-','linewidth',0.5) ; box on; hold on;
xlabel('\it t(ns)', 'Fontsize',14) ;
ylabel('\it T(K)', 'Fontsize',14);

ylim([250 550]);

set(gcf,'position',[100 100 800 200]);

set(gcf, 'paperpositionmode','auto') ;
print(gcf, '-depsc2', 'temp_hist');

% compute and plot density (this can also be done from the restart file itself of course)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temperature grid :
npt=100;
[pdft,bint]=hist(temp,npt); pdft=pdft/length(step);
figure ; hold on; grid on; box on
plot(bint, pdft,'k.-')
% superpose powerlaw
plot(bint,4*bint.^-1,'r')
%plot(bint,1500*bint.^-2,'b')

set(gca, 'xscale','log')
set(gca, 'yscale','log')

xlabel('\it T(K)', 'fontsize',14)
ylabel('\it PDF(T) a.u.', 'fontsize',14);
legend('pdf','power-law fit');
set(gcf, 'paperpositionmode','auto') ;
print(gcf, '-depsc2', 'temp_pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse temperature grid :
npb=100;
[pdfb,binb]=hist(1./(kboltz*temp),npb); pdfb=pdfb/length(step);
figure ; hold on; grid on; box on
plot(binb, pdfb,'k.-')
% superpose powerlaw
plot(binb,0.013*binb.^-1,'r')

set(gca, 'xscale','log')
set(gca, 'yscale','log')

xlabel('\it \beta(kcal/mol)^{-1}', 'fontsize',14)
ylabel('\it PDF(\beta) a.u.', 'fontsize',14);
legend('pdf','power-law fit');
set(gcf, 'paperpositionmode','auto') ;
print(gcf, '-depsc2', 'beta_pdf');




