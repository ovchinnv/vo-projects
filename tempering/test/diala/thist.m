% plot temperature history
%
close all;
kboltz=1.98e-3 ;

file='tempering.series.txt';

d=load(file);

step=d(:,1);
temp=d(:,2);
ener=d(:,3);

dt = 10 * 10 * 1 / 1000000 ; %stat_freq x plugin_freq x timestep (fs) x (ns / fs) => to convert to time (ns)

plot(step*dt, temp) ; box on; hold on;


xlabel('\it t(ns)', 'Fontsize',14) ;
ylabel('\it T(K)', 'Fontsize',14);

% compute and plot density (this can also be done from the restart file itself of course)

% temperature grid :
npt=100;
[pdft,bint]=hist(temp,npt); pdft=pdft/length(step);
figure ; hold on; grid on
plot(bint, pdft,'k.-')
% superpose powerlaw
plot(bint,4*bint.^-1,'r')

set(gca, 'xscale','log')
set(gca, 'yscale','log')


% inverse temperature grid :
npb=100;
[pdfb,binb]=hist(1./(kboltz*temp),npb); pdfb=pdfb/length(step);
figure ; hold on; grid on
plot(binb, pdfb,'k.-')
% superpose powerlaw
plot(binb,0.01*binb.^-1,'r')

set(gca, 'xscale','log')
set(gca, 'yscale','log')

