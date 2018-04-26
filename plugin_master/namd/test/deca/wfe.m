% compute FE from double half-harmonic window simulations
%

close all;

if ~exist('styles')
 styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
end
lw=1;
leg={};

addpath '~/scripts/matlab'; %for xave_

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
 ncv=1;
% set limits
 d0=12 ; d1=32;

 fbw=1;

 nwin=32;

 nbox =2 % number of statistical samples

 df=zeros(nwin, ncv, nbox) ; %derivative matrix

 for j=1:nwin
   dref = d0 + (j-1) / (nwin-1) * (d1-d0);
   % round to 2 decimap pts
   dref=0.01 * fix (dref*100);

   ds(j)=dref ; % save for integration

   drefs=sprintf('%.2f',dref);

   fname=['fbwin_',drefs,'.dat'];
   d=load(fname) ;

   data=reshape(d,ncv,[])' ;

% separate averages from number of samples
   cvs=data(1:2:end,:) ;
   nsamp=data(2:2:end,:) ;
%
   [niter,ncv]=size(cvs);
%
% sample limits and number of boxes
   ie   =niter;
   ib   =2; % skip 1st iter per equilibration
   %ib=round(niter * 0.2);
%
   dfe=zeros(nbox,ncv);
% loop over all cvs, and compute PMF derivative

   for i=1:ncv
% consider only entries with nonzero samples
    inds=find(nsamp(ib:ie,i)>0)+ib-1; % sample for consideration
    if (isempty(inds)) % simulation too short to produce samples within window
     dfe(:,i)=NaN; % will deal with this later
     warning(['Window ',num2str(i),' has no valid samples. Setting F''(',num2str(i),') to NaN.']);
     continue
    end
    nn=nsamp(inds,i);
    dd=cvs(inds,i) - dref; % average displacement from restraint center

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
       error(['End of sample reached in box ',num2str(ibox),' of ',num2str(nbox),' of replica ',num2str(m)]);
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
     g=fzero(f,[-20 20]) ;
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

dfc=reshape(df(:,1,:), nwin, nbox)' ;

% compute center derivative
dfc = 0.5 * ( dfc(:,1:end-1) + dfc(:,2:end) );
fe = [ zeros(nbox,1)  cumsum( dfc.*(ones(nbox,1)*diff(ds)), 2) ];

%============= PLOT FE =============
if ~nofig
 close all;
 figure('position',[200,200,450,350]); hold on; box on;
end
%
for i=1:nbox
 plot(ds,fe(i,1:end),[char(styles(mod(i-1,length(styles))+1)),'x'], 'linewidth', lw)
end
%
fave=mean(fe,1);
fstd=std(fe,1);
%mean
plot(ds,fave,'k--','linewidth',lw);
%std
%plot(alpha,fave+fstd,'k--','linewidth',1);
%plot(alpha,fave-fstd,'k--','linewidth',1);
leg=[leg {['Average']}];

legend(leg,2);
box on;
ylabel('\it F(\alpha) (kcal/mol)', 'fontsize',14);
xlabel('\it x', 'fontsize',14);
%
%xlim([0 1]);
%ylim([0 9]);
set(gcf, 'paperpositionmode', 'auto');
print(gcf, '-dpsc', 'wfe.eps');

