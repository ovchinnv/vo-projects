% process tempering files

file='gdata.dat';

% extract grid
system(['gdata > ',file]);

d=load(file);

bet=d(:,1);
esum=d(:,2);
eesum=d(:,3);
eavg=d(:,4);
eeavg=d(:,5);
wgt=d(:,6);
nsampl=d(:,7);

kb=1.98e-3 ;
temp=1./(kb*bet);

figure; hold on;


%plot(temp, eavg)
%plot(bet, eavg)
%plot(temp, bet.^2 .* ( eeavg-eavg.^2) ) ; % heat capacity
%return



plot(bet, nsampl);

plot(bet,100*bet.^(-1),'ko');

set(gca,'xscale','log')
set(gca,'yscale','log')