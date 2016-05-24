%
close all;
styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};

if (~exist('read'))
 read=1;
end

if (read==1)
%%%%%%%%%% load multiple files %%%%%%%%%%%%
fnames={'curv.dat'};
%
clear crv;
for i=1:length(fnames)
 fname=char(fnames(i));
 if (i==1)
  crv=load(fname);
 else
  crv=[crv; load(fname)];
 end 
end
crv=crv(:,2:end)'; % leave off inxed column
read=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nstring,niter]=size(crv);
%
%
nbox =3;
ib   =1;
%ib   =300;
%ib   =3000;
%ib   =niter;
ie   =niter;
bsize=ceil( (ie-ib+1)/nbox);

%bsize=90

%
leg={};
nsample=floor((ie-ib+1)/bsize)+sign(mod(ie-ib+1,bsize));
cav=zeros(nsample,nstring);
j=1;
for i=ib:bsize:ie
 cav(j,:)=mean(crv(:, i : i+min(bsize-1,ie-i) ),2);
 j=j+1;
 leg=[leg {['iteration ',num2str(i-1)]}];
end 

%figure; hold on;box on;
figure('position',[200,200,600,250]); hold on; box on;


alpha=([1:nstring]-1)/(nstring-1);

cavs=zeros(nsample,nstring);
for i=1:nsample
% smooth curvature
 dfilter=3;
% cavs(i,:)=smooth2(alpha,cav(i,:),dfilter);
%
 plot(cav(i,1:end),[char(styles(mod(i-1,length(styles))+1)),'x'], 'linewidth', 2)
% plot(alpha,cavs(i,:),['*',char(styles(mod(i-1,length(styles))+1))], 'linewidth', 2)
end 

box on;
ylabel('\kappa (Ang^{-1})', 'fontsize',14);
xlabel('\it \alpha ', 'fontsize',14);

%axis([0 1 0 4.5]);
set(gcf, 'paperpositionmode', 'auto');
%print(gcf, '-dpsc', 'curvp.eps');
%print(gcf, '-djpeg100', 'fe64_p1.jpg');
