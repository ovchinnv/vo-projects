% simple test case of stiffness fitting
%
% three particles connected by two springs :  o---/\/\/---o---/\/\/---o
%                                                 k1=k12      k2=k23


% variables
syms k1 k2 k3
syms V k11 k12 k13 k21 k22 k23 k31 k32 k33
syms x1 x2 x3

k12=-k1  ; k21=k12 ; k13=-k3
k23=-k2  ; k32=k23 ; k31=k13

% partial stiff matrix
K=[0 k12 k13 ; 0 0 k23 ; 0 0 0 ];
%symmetrize:
K=K+transpose(K); 
np=3; for i=1:np; K(i,i)=-sum(K(i,:)); end 
K

% define potential energy
X=[x1 ; x2 ; x3];

V = 1/2 *  transpose(X) * K * X


% put in actual numbers

k1=1.;
k2=1.5;
k3=.5; % spring between 1 and 3

%k1=.5;
%k2=1.;
%k3=1.5;


%%%

K=eval(K)
V=eval(V)

%now apply a constraint to keep the center of geometry fixed: x1+x1+x3=0 => x3=-x1-x2

%projection matrix:
ndof=np-1 ; %remove one dof
P=[eye(ndof)  -ones(ndof,1) ];

Kp=P*K*P'


Kpi=inv(Kp) ; % compliance matrix (also unnormalized fluctuation matrix)

%full fluctuation matrix (call this K inverse):
Ki = P' * Kpi * P

% compute this pseudo inverse by eigenvalue decomposition

[v,e]=eig(K); % first evalue is zero -- change it
e=1./diag(e); e=diag([0; e(2:end)]);
Ki2 = v * e * v' ; %this is equal to Ki above

% the above suggests that we know what we are doing
% now let us attempt to fit the coefficients of the stiffness matrix using the NR algorithm I developed
%%%%%%%%% basically just pasting from stiff.m %%%%%%%%%%%%%%%

Ko=K; % original stiffness

kT=0.59; %normalization energy
alpha=1; % * 3 * kT
C = Ki * alpha ; 

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

bonds = [ 1 2 ; 2 3 ; ]%3 1 ] ;
bonds = [ 1 2 ; 2 3 ;3 1 ] ;
%bonds = [ 2 3 ; 3 1 ] ;

natom=np;
niter=30;
%niter=1

errs=zeros(1,niter);

nbonds=size(bonds,1);
ks=0.1*ones(nbonds,1); % initial guess for force constants

ks(3)=0.;
%ks(1)=0.1;
%

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

% derivative of the (reduced) stiffness matrix w.r.t. spring constants
Kp1=cell(1,nbonds);
Kpi1=cell(1,nbonds);
Kpj1=cell(1,nbonds);
err1=zeros(nbonds,1);

% derivative matrices
H1=zeros(nbonds,nbonds);
H2=zeros(nbonds,nbonds);


% iterate
for iter=1:niter

K=sparse(bonds(:,1),bonds(:,2),-ks,natom,natom);
K=K+K'; % symmetrize

for i=1:natom
 K(i,i)=-sum(K(i,:));
end

%
for k=1:nbonds
 Kp1(k)={sparse(P*K1(:,:,k)*P')};
end

% NR-iteration
Kp=P*K*P';   % reduced stiffness matrix
Kpi=inv(Kp); % inverse is essentially full
Ki = P' * Kpi * P; % full fluctuation matrix
%

normC = sum(sum(C.*C)) ;
def=(C-alpha*Ki);
err=(sum(sum(def.*def))) / normC ;
%
% derivatives of the r. stiffness inverse
for k=1:nbonds
 k
 L=cell2mat(Kp1(k));
 M =  - Kpi * L * Kpi; 
 M1 = P' * M * P;
 Kpi1(k)={M };
 Kpj1(k)={M1};

% also: derivative of err w.r.t. bond constants
 err1(k) = sum(sum(def.*M1 ))/normC;
% also: compute the derivative matrices (basically, Hessians)
 for l=1:k-1
  N1=cell2mat(Kpj1(l));
  H1(k,l)=sum(sum(M1.*N1));
  H1(l,k)=H1(k,l);
%
  Q=cell2mat(Kp1(l));
  N=cell2mat(Kpi1(l));
  R = - P' * sparse((M*Q + N*L)) * Kpi * P;

  H2(k,l)=sum(sum(def.*R));
  H2(l,k)=H2(k,l);
  
 end
 
 H1(k,k)=sum(sum(M1.*M1));
 H2(k,k)=-2*sum(sum(def.* ( P' * sparse(M*L) * Kpi * P)));
%
end

H = ( -alpha*H1 + H2 ) / normC; % this is the derivative of teh objective function w.r.t dof (fc. constants)

% compute H by finite diffs -- pass !
if (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=0.000001;
Hfd=zeros(nbonds);
for l=1:nbonds

 ks(l)=ks(l)+h;
 K=sparse(bonds(:,1),bonds(:,2),-ks,natom,natom);
 ks(l)=ks(l)-h; % revert

 K=K+K'; for i=1:natom ; K(i,i)=-sum(K(i,:)); end
 Kp=P*K*P';   % reduced stiffness matrix
 Kpi=inv(Kp); % inverse is essentially full
 Ki = P' * Kpi * P; % full fluctuation matrix
 def=(C-alpha*Ki);
 
 for k=1:nbonds
  L = cell2mat(Kp1(k));
  M =  - Kpi * L * Kpi; 
  M1 = P' * M * P;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Hfd(k,l) = ( sum(sum(def.*M1 ))/normC - err1(k) ) / h ;
 end
%%% for 2nd order diffs : 
 ks(l)=ks(l)-h;
 K=sparse(bonds(:,1),bonds(:,2),-ks,natom,natom);
 ks(l)=ks(l)+h; % revert
%
 K=K+K'; for i=1:natom ; K(i,i)=-sum(K(i,:)); end
 Kp=P*K*P';   % reduced stiffness matrix
 Kpi=inv(Kp); % inverse is essentially full
 Ki = P' * Kpi * P; % full fluctuation matrix
 def=(C-alpha*Ki);
 
 for k=1:nbonds
  L = cell2mat(Kp1(k));
  M =  - Kpi * L * Kpi; 
  M1 = P' * M * P;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Hfd(k,l) = Hfd(k,l) - ( sum(sum(def.*M1 ))/normC - err1(k) ) / h ;
 end

end
Hfd=Hfd/2;
( H-Hfd ) ./ H 
return
end % fd


dks=H\err1; 
%dks=dks*0.1 ;
ks = ks - dks ; %ks=max(ks,0.0); ks=min(ks,2.0); 

%full fluctuation matrix (call this K inverse):
Ki = P' * Kpi * P;
full(Ki)
ks
err
errs(iter)=err;

end % iter

close
plot(errs,'k*'); ylim([0 1]);
