% compute FE from double half-harmonic window simulations
%

close all;

if ~exist('styles')
 styles={'r-x','g-v','b-s','m-*','c-','k-','r--','g--','b--','m--','c--','k--'};
end
lw=1;
leg={};

%addpath '~/scripts/matlab'; %for xave_

kboltz=1.987191d-3;        % boltzmann constant
Temp=300;

if ~exist('nofig')
 nofig=0;
end

if ~exist('read')
 read=1;
end

if (read)
%%%%%%%%%% process windows
 fbw=0.5; % only applies to the first position component
 iwin=1; % can be 0 or 1 depending on whether the equilibrium point is included
 nwin=7;
% nsamples=4;
% [status, result]=system('grep "will quit" pmf3.out | tail -n1 | awk ''{print $3}'''); nsamples=str2num(result)-1 ; nsamples=nsamples-31 ; % screwed up counts due to crash
% [status, result]=system('grep "will quit" pmf2.log | tail -n1 | awk ''{print $3}'''); nsamples=str2num(result)-251 ;
 nbox=2; % number of statistical samples
 rcind=1 ; %cv index corresponding to the reaction coordinate

 clear df ;
 rc=zeros(1,nwin+1-iwin); % reaction coordinate
 for j=1:nwin+1-iwin;

   ncv=4;
   fname=['data/fbwin',num2str(j-1+iwin),'.dat'];
   dnew=load(fname) ;
   if exist('nsamples')
    d=dnew(1:2*ncv*nsamples,:); % i.e. 2*ncv lines per samples
   else
    d=dnew;
   end
%
   data=reshape(d,ncv,[])' ;
%
% extract the only relevant CV
   cvs=data(1:2:end,rcind) ;
   nsamp=data(2:2:end,rcind) ;
%
   [niter,ncv]=size(cvs); % redefine ncv for rest of script
%
% sample limits and number of boxes
   ie   = niter
   ib   = 2; % skip per equilibration (10 is the minimum -- old restraints are gradually advanced over 10 iterations)
%   ib = 200 ;
   ib=max(1,round(ie * 0.5));
   niters(j)=niter ; % store number of iterations
%
   if (~exist('df'))
    df=zeros(nwin,ncv,nbox);
   end
% loop over all cvs, and compute PMF derivative
% first, need the reference position -- open cv.dat file :
   cvs0 = load(['cv',num2str(j-1+iwin),'.dat']);
   cvs0 = cvs0(rcind);
   rc(j) = cvs0;
%
   for i=1:ncv
% consider only entries with nonzero samples
    inds=find(nsamp(ib:ie,i)>0)+ib-1; % sample for consideration
    if (isempty(inds)) % simulation too short to produce samples within window
     df(j,i,:)=NaN; % will deal with this later
     warning(['Window ',num2str(j),' has no valid samples. Setting F''(',num2str(i),') to NaN.']);
%     return
     continue
    end
    nn=nsamp(inds,i);
    dd=cvs(inds,i) - cvs0(i); % average displacement from restraint center

    dbox = 0.5*fbw ; % half-width
    dd = (dd+dbox)/(2*dbox) ; % normalize data to the interval [0 1]
%
% check to make sure there are no outliers
    if ( find ( abs(dd - 0.5) > 0.5 , 1) )
     error(['ERROR : one or more averages for replica ', num2str(i),' is outside the [0,1] interval. Aborting.']);
    end
% cumulative number of samples
    nnc=cumsum(nn);
    bsize=ceil( nnc(end)/nbox ); %block size
    nnc=nnc/bsize ; % normalize sample count by block size and look for block indices, 1, 2, 3 ... etc
%
    ibeg=1;
    for ibox=1:nbox
     iend=find(nnc>ibox,1)-1 ; % last index
     if (isempty(iend))
      if ibox < nbox
       error(['End of sample reached in box ',num2str(ibox),' of ',num2str(nbox),' of replica ',num2str(j)]);
       return
      else
      iend=length(nnc) ;%take the last sample
      end
     end % iend
% combined average
     ddcomb=sum(dd(ibeg:iend).*nn(ibeg:iend))/sum(nn(ibeg:iend)) ;
% define ibeg for next block (if any)
     ibeg=iend+1;
% solve for FE slope
     f=@(x) xave_s(x)-ddcomb ;
     g=fzero(f,[-100 100]) ;
% stop if a problem with gamma
     if(isnan(g)) ; error('Cannot find gamma: got NaN'); return ; end
% compute df from gamma
     df(j, i, ibox) = g * kboltz * Temp / fbw;
    end % over ibox boxes
   end % over all cvs
%
 end
 read=0;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we know that only one cv is present
dfc0=reshape(df(:,1,:), nwin+1-iwin, nbox)' ;
rc0=rc;
if (0)
% append 0th point, at which we assume df to be zero (stable eq. simulation)
% however, this may not be true is the initial condition is not perfect, or the CVs are not perfect
% in some cases it is best to omit 0th point
 dfc0=[zeros(nbox,1) dfc0];
 rc0 =[2*rc0(1) - rc0(2), rc0]; % assume uniform interval
end

% trapezoidal rule
% compute center derivative
dfc = 0.5 * ( dfc0(:,1:end-1) + dfc0(:,2:end) );
%integral
fe = [ zeros(nbox,1)  cumsum( dfc.*(ones(nbox,1)*diff(rc0)), 2) ];

% to offset from the end : this should expose that the initial replicas are not converging
feoff=mean(fe(:,end/2:end),2);
feoffs=repmat(feoff,1,size(fe,2));
%fe=fe-feoffs;
%============= PLOT FE =============
if ~nofig
 close all;
 figure('position',[200,200,450,350]); hold on; box on;
end
%
for i=1:nbox
 plot(rc0,fe(i,1:end),[char(styles(mod(i-1,length(styles))+1)),''], 'linewidth', lw)
end
%
fave=mean(fe,1);
fstd=std(fe,1);
%mean
%plot(rc0,fave,'k--','linewidth',lw);
%std
%plot(alpha,fave+fstd,'k--','linewidth',1);
%plot(alpha,fave-fstd,'k--','linewidth',1);
%leg=[leg {['Average']}];

legend(leg,4);
box on;
ylabel('\it F(\alpha) (kcal/mol)', 'fontsize',14);
xlabel('\it x', 'fontsize',14);
%
%xlim([0 1]);
%ylim([0 9]);
%set(gcf, 'paperpositionmode', 'auto');
print(gcf, '-dpsc', 'wfe.eps');
%
feoff=min(fe,[],2);
feoffs=repmat(feoff,1,size(fe,2));
fe-feoffs
%fe-repmat(fe(:,3),1,size(fe,2))
mean(ans)
%pause(10)

save -mat wfe.mat fe feoffs fave fstd rc0 niters

