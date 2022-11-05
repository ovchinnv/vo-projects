% compute FE from double half-harmonic window simulations
%
if exist('OCTAVE_VERSION')
 clear graphics_toolkit;
 graphics_toolkit('gnuplot');
end

addpath '~/scripts/matlab'; %for xave_

fefile='wfe.mat' ;

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
 ncv=2;
% set limits
 phi0=-180 ; phi1 = 180 ;
 psi0=-180 ; psi1 = 180 ;

 fbw=[5, 5];
 fbwrad=fbw * pi / 180 ;

% fbwrad=[0.175 0.175];

 nphi=20;
 npsi=20;

 nbox =1; % number of statistical samples

 df=zeros(nphi, npsi, ncv, nbox) ; %derivative matrix

 for j=1:nphi
  phi = phi0 + (j-1) / nphi * (phi1-phi0);
  for k=1:npsi
   psi = psi0 + (k-1) / npsi * (psi1-psi0);
%
   if (phi~=0)
    phis=sprintf('%.2f',phi);
   else
    phis='0';
   end
   if (psi~=0)
    psis=sprintf('%.2f',psi);
   else
    psis='0';
   end

   phirad=phi/180*pi;
   psirad=psi/180*pi;
%
   fname=['fbwin_',phis,'_',psis,'.dat'];
   d=load(fname) ;

   data=reshape(d,ncv,[])' ;

% separate averages from number of samples
   cvs=data(1:2:end,:) ;
   nsamp=data(2:2:end,:) ;
%
   [niter,ncv]=size(cvs);
% restraint centers for each replica
   cvs0=[phirad psirad] ;
% sample limits and number of boxes
   ie   =niter;
   ib   =2; % skip as equilibration
%   ib=round(niter * 0.25);
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
    dd=cvs(inds,i) - cvs0(i); % average displacement from restraint center
% account for periodicity
    dd=mod(dd, 2*pi) ; ind=find(dd>pi); dd(ind)=dd(ind)-pi-pi ;

    dbox = 0.5*fbwrad(i) ; % half-width
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
     g=fzero(f,0) ;
% stop if a problem with gamma
     if(isnan(g)) ; error('Cannot find gamma: got NaN'); return ; end
% compute df from gamma
     df(j, k, i, ibox) = g * kboltz * Temp / fbwrad(i);
    end % over ibox boxes
   end % over all cvs
%
  end % k -psi
 end % j -phi
 read=0;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integrate derivatives to get PMF
% define grid
x = ( phi0 + ( (1:nphi) - 1 ) / nphi * (phi1-phi0) ) * pi / 180 ;
y = ( psi0 + ( (1:npsi) - 1 ) / npsi * (psi1-psi0) ) * pi / 180 ;
dx = diff(x);
dy = diff(y) ;
%
f = zeros( nphi, npsi );
if ( length(size(df))==4 )
% dfm=mean(df,4); % derivatives averaged over boxes
else % otherwise the 4th dimensions may be absent, since it is trivial
 dfm=df ;
end
fx=dfm(:,:,1);
fy=dfm(:,:,2);

% corner
f(1,1)=0 ;
% boundaries
for i = 2 : nphi
 f(i,1) = f(i-1,1) + 0.5 * ( fx(i-1,1) + fx(i,1) ) * dx (i-1) ;
end
%
for j = 2 : npsi
 f(1,j) = f(1,j-1) + 0.5 * ( fy(1,j-1) + fy(1,j) ) * dy (j-1) ;
end
%
for i=2:nphi
 for j=2:npsi
  f(i,j) = 0.5 * ( f(i-1,j) + f(i,j-1) ) + 0.25 * ( (fx(i-1,j) + fx(i,j))*dx(i-1) + (fy(i,j) + fy(i,j-1))*dy(j-1) ) ;
 end
end

% shift to zero mean (or elsewhere):
%f=f-mean(mean(f));

todeg=180/pi;
pcolor(y*todeg,x*todeg,f'); shading faceted ;
colorbar ;

print(gcf, '-depsc2','wfe.eps')
%save('-mat', 'diala36-pbc20x20.mat');

%f4=[ f' f' ; f' f' ] ; % extend periodically
%pcolor(f4) ; shading interp ; colorbar ;
