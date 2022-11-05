% compute FE from double half-harmonic window simulations
if ( exist('OCTAVE_VERSION') )
 qoct=1;
end
if qoct
 graphics_toolkit('gnuplot')
% graphics_toolkit('qt')
endif
%
addpath '~/scripts/matlab'; %for xave_
%
fefile='wfe.mat' ;

kboltz=1.987191d-3;        % boltzmann constant
Temp=300;

if ~exist('nofig')
 nofig=0;
end

if ~exist('read')
 read=1;
end

qper=1 ;% whether to include periodicity implicitly ; this is more correct
qlsq=1 ;% least squares integration -- theoretically, the most accurate

todeg=180/pi;
torad=1./todeg;

if (read)
%%%%%%%%%% process windows
 ncv=2;
% set limits
 phi0=-180 ; phi1 = 180 ;
 psi0=-180 ; psi1 = 180 ;

 fbw=[5, 5]; % flat bottom widths in degrees
 fbwrad=fbw * torad ; % in radians

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
%
   phirad=phi*torad;
   psirad=psi*torad;
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
   ib   =2; % skip 1st iter per equilibration
   %ib=round(niter * 0.2);
%
   dfe=zeros(nbox,ncv);
% loop over all cvs, and compute PMF derivative
%
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
if (nbox>1)
 dfm=mean(df,4); % derivatives averaged over boxes
else
 dfm=df;
end

fx=dfm(:,:,1);
fy=dfm(:,:,2);

if (~qper) % direct computation -- approximate ;  note that we are not accounting for periodicity
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
%
elseif (~qlsq)
% construct equation matrix and account for periodicity
% note that even this seems to be a somewhat adhoc construction
% you can see that something is "funny" because the BC point often creates an outlier
% will need to debug this to make sure it is correct ; the artifacts could be due to the system being overdetermined by periodicity
n=nphi*npsi;
M=eye(n,n);
rhs=zeros(n,1);
for i=1:nphi
 im=mod( (i-1)-1, nphi)+1; % 1 maps to nphi+1, and there is no 0 index ; thus, we subtract 1 to get 0 ~ nphi , here: (i-1)-1 mod nphi  + 1 (add one at the end)
 io=(i-1)*npsi ; % compute offset due to previous points
 for j=1:npsi
  o=io+j ;           % center index in M
  l=(im-1)*npsi+j ; % slowly varying index
  jm=mod( (j-1)-1, npsi)+1;
  b=io+jm; % quickly varying index
  M(o,o)=2 ; % NOTE : already set to 1 via eye() above
  M(o,l)=-1 ;
  M(o,b)=-1 ;
  rhs(o)=0.5 * ( (fx(im,j) + fx(i,j))*dx(1) + (fy(i,j) + fy(i,jm))*dy(1) ); % adding two solutions obtained using two 1D integrations
%if (o==2) ;return;end
 end % nphi
end % npsi
% apply "bc"
% add a row to get use LSQ
M(end+1,:)=1 ; rhs(end+1)=0; % necessary to add another eq ; otherwise a point gets corrupted
% without BC get garbage
%
%f=inv(M)*rhs;
f=M\rhs; % note that in octave can automatically switch to approximate methods ; noting this because the periodic system seems to be overdetermined !
% compute residue to check solution quality:
residual = M*f - rhs;
fprintf('%s%12.5e\n','Maximum residual ',max(abs(residual)))
f=reshape(f,nphi, npsi)'; %f=f-mean(mean(f));

else % least squares - similar matrix construction to above, but the individual gradient eqs are separate for each dim
 n=nphi*npsi;
 M=zeros(2*n+1,n);
 rhs=zeros(2*n+1,1);
 for i=1:nphi
  im=mod( (i-1)-1, nphi)+1; % 1 maps to nphi+1, and there is no 0 index ; thus, we subtract 1 to get 0 ~ nphi , here: (i-1)-1 mod nphi  + 1 (add one at the end)
  io=(i-1)*npsi ; % compute offset due to previous points
  for j=1:npsi
   o=io+j ;           % center index in M
   l=(im-1)*npsi+j ; % slowly varying index
   jm=mod( (j-1)-1, npsi)+1;
   b=io+jm; % quickly varying index
% x-grad
   M(o,o)=1 ;
   M(o,l)=-1 ;
   rhs(o)=0.5 * ( (fx(im,j) + fx(i,j))*dx(1) );
% y-grad
   M(2*o,o)=1 ;
   M(2*o,b)=-1 ;
   rhs(2*o)=0.5 * ( (fy(i,jm) + fy(i,j))*dy(1) );
%if (o==2) ;return;end
  end % nphi
 end % npsi
% M(end+1,:)=1 ; rhs(end+1)=0; % average solution
 f=M\rhs; % note that in octave can automatically switch to approximate methods ; noting this because the periodic system seems to be overdetermined !
% f=pinv(M)*rhs ;
 residual = M*f - rhs;
 fprintf('%s%12.5e\n','Maximum residual ',max(abs(residual)))
 f=reshape(f,nphi, npsi)';
 f=f-mean(f(:));
end

save('-mat',fefile) ;
if (nofig==1)
 return
end
%
%figure(1);clf
figure ;clf
%pcolor(y*todeg,x*todeg,f');
mycolor(y*todeg,x*todeg,f');

%f4=[ f' f' ; f' f' ] ; % extend periodically
%x2=[ x x+2*pi ]; y2=[ y y+2*pi];
%pcolor( y2*todeg, x2*todeg, f4 ) ;

shading interp ;
colorbar ; colormap(jet)

axis tight ; set(gca,'tickdir','out', 'fontsize',16)
ylabel('\psi')
xlabel('\phi')

