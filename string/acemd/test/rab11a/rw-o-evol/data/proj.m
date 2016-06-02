% 4.14: plot projection variable profiles for FTSM
%

close all;

styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
%%%%%%%%%% load multiple files %%%%%%%%%%%%
fnames={'proj.dat'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('read'))
 read=1;
end
%read=1;
if (read==1)
%
clear proj;
for i=1:length(fnames)
 fname=char(fnames(i));
 if (i==1)
  d=load(fname);
 else
  d=[d; load(fname)];
 end 
end
%
[niter,m]=size(d);
ncv=2; % two projection variables : parallel and perpendicular
niter=niter/ncv;

dpari=d(1:2:end,:);
dperpi=d(2:2:end,:);

% remove normalization by arclength
% read arclength file
%arcfile = 'arcl3.dat' ;
%ds = load(arcfile) ;
%ds=mean(ds(:, 2:end),1) ;
%dsm=mean(ds) ; % average value
dsm=1 ;

dperpi=dperpi*2*dsm ; 
dperpi(:,1) = dperpi(:,1) *0.5 ; % endpoints were divided by dsm, inner points by 2*dsm
dperpi(:,m) = dperpi(:,m) *0.5 ; 
read=0;
end % read
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ie   =niter; % last sample to take
ib   =1;     % first sample
%ib=round(niter * 0.2); % might want to discard initial part of sample
%ib=100;
%
nbox =2;     % number of boxes
bsize=ceil( (ie-ib+1)/nbox); % # samples in bos (used internally below)
%
lw=1.2;
%
leg={};
nsample=floor((ie-ib+1)/bsize)+sign(mod(ie-ib+1,bsize));
dpar =zeros(nsample,m);
dperp=zeros(nsample,m);
j=1;
for i=ib:bsize:ie
 dpar(j,:) =mean(dpari (i : i+min(bsize-1,ie-i), : ),1);
 dperp(j,:)=mean(dperpi(i : i+min(bsize-1,ie-i), : ),1);
 j=j+1;
 leg=[leg {['iteration ',num2str(i-1)]}];
end 

fpar=figure('position',[50,150,600,250]); hold on; box on;
fprp=figure('position',[650,150,600,250]); hold on; box on;
alpha=[0:m-1]; alpha=alpha/alpha(end);

for i=1:nsample
 figure(fpar)
 plot(alpha,dpar(i,1:end),[char(styles(mod(i-1,length(styles))+1)),'x'], 'linewidth', lw)
 figure(fprp)
 plot(alpha,dperp(i,1:end),[char(styles(mod(i-1,length(styles))+1)),'o'], 'linewidth', lw)
end 

figure(fpar)
legend(leg,2);
box on;
ylabel('\rm d_{||} (\alpha) ', 'fontsize',14);
xlabel('\it \alpha', 'fontsize',14);
xlim([0 1]);
%ylim([-15 15]);
set(gcf, 'paperpositionmode', 'auto');
%print(gcf, '-dpsc', 'dpar.eps');
%
figure(fprp)
legend(leg,3);
box on;
ylabel('\rm d_{\perp} (\alpha) ', 'fontsize',14);
xlabel('\it \alpha', 'fontsize',14);
xlim([0 1]);
%ylim([-15 15]);
set(gcf, 'paperpositionmode', 'auto');
%print(gcf, '-dpsc', 'dperp.eps');

%
return
% plot evolution of FE drop
iw=10;
dfe=zeros(1,nsample);
for i=1:nsample
 dfe(i)=min(fe(i,end-iw+1:end) - fe(i,1:iw));
end
figure(3);
plot(dfe,'k-x');
ylabel('\it Endpoint Free Energy Difference(kcal/mol)', 'fontsize',14);
xlabel('\it Iteration \# ', 'fontsize',14);

