% plot some averages
gdata;

close all;

% interpolate average using the method of Zhang&Ma10:
n=length(bet) ;
db=50; % number of interpolation bins on either side of central bin
for i=1:n
 ib=max(1, i-db);
 ie=min(i+db,n);
 x=bet (ib:ie)';
 y=eavg(ib:ie)';
 w=wgt (ib:ie)';
% determine fit
 m=length(x);
 ws=sum(w);
 wx=w.*x;
 xs=sum(wx) ; 
 ys=w*y' ;
 xy=wx*y' ;
 x2s=wx*x' ;
 d = ws*x2s-xs^2 ;
 a = (x2s*ys-xs*xy)/d;
 b = (ws*xy-xs*ys)/d;
% evaluate :
 eint(i)=a+b*bet(i) ;
end

% average energy
figure ; hold on; box on;
plot(bet, eavg,'k:') % noisy average from the bins
plot(bet, eint,'r','linewidth',2) % smoothed interpolated average
legend('bin average', ['smooth average (',num2str(2*db+1),' bins)']);
ylabel('\it <E>(kcal/mol)', 'fontsize',14);
xlabel('\it \beta(kcal/mol)^{-1}', 'fontsize',14);

%set(gca,'xscale','log')
%set(gca,'yscale','log')
set(gcf, 'paperpositionmode','auto') ;
print(gcf, '-depsc2', 'energy');

% heat capacity
figure ; hold on; box on;

betc=0.5*( bet(1:end-1)+bet(2:end)) ;% midpoints
cvm=-kb*betc.^2 .* (diff(eint')./diff(bet)) ; % from average energy derivative
plot(betc  , cvm, 'r');

cvv=kb*bet.^2 .* (eeavg-eavg.^2) ; % from variance
plot(bet, cvv, 'k-');

legend(['from smoothed average (',num2str(2*db+1),' bins)'],'from variance');
ylabel('\it C_V(kcal/mol/K)', 'fontsize',14);
xlabel('\it \beta(kcal/mol)^{-1}', 'fontsize',14);

set(gcf, 'paperpositionmode','auto') ;
print(gcf, '-depsc2', 'CV');


format long;
mean(cvv)
mean(cvm)

% temperature distribution

figure ; hold on ; box on ; grid on;
plot(bet, nsampl);
plot(bet,13000*bet.^(-1),'k-');

set(gca,'xscale','log')
set(gca,'yscale','log')




