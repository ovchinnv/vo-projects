% parametrize stiffness matrix
%
%
addpath ~/scripts/matlab/matdcd-1.0

dcdname='cg_large_md.dcd';

if ~exist('read')
 read=1;
end

if (read==1)

h=read_dcdheader(dcdname);
nframes=h.NSET;
natom=h.N;

xyz=zeros(nframes,3*natom);
xyz=readdcd(dcdname,1:natom);

%
% can superpose structure using Procrustes analysis (but need to manually weight the atoms first)
%
% (1) read datafile produced with 'cg'
%
fname='clpx-3hws-cgla';
load([fname,'.dat'],'-mat');
%
% reorganize DCD data :
% inelegant but correct way
%
rr=zeros(natom,3,nframes);
for i=1:nframes
 rr(:,1,i) = xyz (i,1:3:end);
 rr(:,2,i) = xyz (i,2:3:end);
 rr(:,3,i) = xyz (i,3:3:end);
end
%
% to mass-weight coordinates:
m(:)=1; % uniform weights
% NOTE: with the uniform weights as above, the results of the alignments are the same as those in VMD
% this is expected for uniform weights; (the COM is the same)
% when the weights are variable, the best fit alignment by Procrustes below involves a COM using squares of the weights
w=m/sum(m); w=sqrt(w);

rrw=zeros(natom,3,nframes);
for i=1:natom
 rrw(i,:,:)=w(i)*rr(i,:,:);
end
%
% mass-weight frames
%

close all; styles='rgbmck'
msd=zeros(1,nframes);
rave=rrw(:,:,1);
niter=2; % three is usually sufficient, but two is good enough for government work
msds=zeros(1,niter);
for iter=1:niter % iterations, as in generalized Procrustean
 for i=1:nframes
  rold=rrw(:,:,i);
  [d,rnew,trans]=procrustes(rave,rold,'Scaling',false,'Reflection',false);% not quite right b/c the translation done with sqrt(w) not w
  msd(i)=( sum(sum( (rave-rnew).^2) ) );
  rrw(:,:,i)=rnew;
 end
 rave=mean(rrw,3);
% plot(sqrt(msd),styles(iter)); hold on;;
 msds(iter)=sum(msd);
end
msds;

% since we do not want mass-weighted coordinates, scale back
for i=1:natom
 rr(i,:,:)=rrw(i,:,:)/w(i);
 rave(i,:)=rave(i,:)/w(i);
end
%
%
% write an oriented file, just to see what it looks like -- fine in VMD
%dcdout='cg_large_md_orient.dcd';
%writedcd(dcdout,squeeze(rr(:,1,:)),squeeze(rr(:,2,:)),squeeze(rr(:,3,:)));
%
% build correlation matrix
%
C=zeros(natom, natom);
%compute deviations
for i=1:nframes
 rr(:,:,i) = rr(:,:,i) - rave(:,:);
end

for i=1:natom
 fprintf('%d\n',i)
 for j=i:natom
  C(i,j) = ( squeeze(rr(i,1,:))' * squeeze(rr(j,1,:)) + ...
             squeeze(rr(i,2,:))' * squeeze(rr(j,2,:)) + ...
             squeeze(rr(i,3,:))' * squeeze(rr(j,3,:)) ) / nframes;
  C(j,i)=C(i,j);
 end
end

read=0;
end % if read


rmsf=sqrt(diag(C)); % RMSF -- (almost, to within three points of decimal prec.) same as in VMD (if using uniform weights)
%figure;
%plot( rmsf );
%
% fit the stiffness matrix to the fluctuations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) apply the constraint tha all displacements sum up to zero 
% (or mass-weighted displacements) and eliminate last row & column
%

Y=C(1:end-1,1:end-1); % fluctuation matrix
% initial guess for stiffness matrix:

nbonds=size(bonds,1);
ks=1*ones(1,nbonds); % initial guess for force constants
%
K=sparse(bonds(:,1),bonds(:,2),-ks,natom,natom);
K=K+K'; % symmetrize

for i=1:natom
 K(i,i)=-sum(K(i,:));
end
% derivative of the stiffness matrix w.r.t. force constants
K1=zeros(natom, natom, nbonds);
for k=1:nbonds
 i=bonds(k,1); j=bonds(k,2); 
 K1(i,i,k)= 1; K1(j,j,k)= 1;
 K1(i,j,k)=-1; K1(j,i,k)=-1;
end
ndof=natom-1;
% projection matrix to eliminate last dof:
P=sparse([eye(ndof) -ones(ndof,1)]); % constant

kT=0.59; %normalization energy

% derivative of the reduced stiffness matrix w.r.t. spring constants
%Kp1=zeros(ndof, ndof, nbonds);
%Kp1=cell(1,nbonds); for i=1:nbonds;  Kp1(i)={zeros(ndof, ndof)}; end
%Kpi1=cell(1,nbonds); for i=1:nbonds;  Kpi1(i)={zeros(ndof, ndof)}; end
Kp1=cell(1,nbonds); %for i=1:nbonds;  Kp1(i)={sparse(ndof, ndof)}; end
Kpi1=cell(1,nbonds);% for i=1:nbonds;  Kpi1(i)={sparse(ndof, ndof)}; end

err1=zeros(1,nbonds);
%
for k=1:nbonds
 Kp1(k)={sparse(P*K1(:,:,k)*P')};
end

% derivative matrices
H1=zeros(nbonds,nbonds);
H2=zeros(nbonds,nbonds);

% begin NR-iterations
Kp=P*K*P';   % reduced stiffness matrix
Kpi=full(inv(Kp)); % inverse is essentially full

def=(Y-3*kT*Kpi);
err=(sum(sum(def.*def)))/ndof/ndof;
return
% derivatives of the r. stiffness inverse
for k=1:nbonds
 k
 L=cell2mat(Kp1(k));
 M =  - Kpi * L * Kpi; 
 Kpi1(k)={M};
% also: derivative of err w.r.t. bond constants
 err1(k) = sum(sum(def.*M));
% also: compute the derivative matrices (basically, Hessians)
 for l=1:k-1
  N=cell2mat(Kpi1(l));
  H1(k,l)=sum(sum(M.*N));
  H1(l,k)=H1(k,l);
%
  Q=cell2mat(Kp1(l));
  R = - sparse((M*Q + N*L)) * Kpi;

  H2(k,l)=sum(sum(def.*R));
  H2(l,k)=H2(k,l);
  
 end
 
 H1(k,k)=sum(sum(M.*M));
 H2(k,k)=-2*sum(sum(def.*( sparse(M*L) * Kpi)));
%
end

H = 3*kT*H1 - H2;
ks = ks + H\err1;





