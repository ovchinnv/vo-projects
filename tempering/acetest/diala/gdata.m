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
