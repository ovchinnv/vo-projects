%
kb=1.3806503e-23; % J/K
Na=6.0221e23;     %
hbar=1.05457148e-34; % kg m^2 / s   (J . s)
jle=1/4.184/1000;      % (kcal)
T=300;            % K
c=2.998e10;        % cm/s in vacuum
kb2=kb*Na*jle;    % kb2 * T = 0.6
kbt=kb2*T;

h=hbar*2*pi*Na*jle; %   so that h * pf2 gives kcal/mol 

% compute entropy from covariance matrix eigenvalues
%
fname='rab11a_eval.dat';
%fname='carma.PCA.eigenvalues.dat';

[ev]=textread(fname,'','emptyvalue',NaN); % puts missing values where there is a space followed by newline

ev=reshape(ev',[],1); %make into a single column ; by default takes each successive row, but we want to take each successive column
ev=ev(find(~isnan(ev)));

% note : checked these evals against those computed by matlab from covariance matrix -- same to within machine precision

% remove eigenvalues that are nearly zero
%tol=1e-10; ev=ev(find(abs(ev)>tol));
tol=1e-10; ev=ev(find(ev>tol)) % only throw away negatives;

% now have raw eigenvalues from the diagonalization of the covariance matrix (units Ang^2 amu)
% convert frequencies to s^-1 units : 
%
fact  = 1./sqrt(jle * 1e-23) / (2*pi) ; % 2.0455e13 ! this essentially converts kilocalories to Joules in kT
nu_s  = sqrt(kbt./ev) .* fact
nu_ps = nu_s*1e-12
nu_icm= (nu_s/c)'

% Andricioaei formula
pe=nu_s*h/kbt;           % energy of the quasiharmonic modes divided by kBT
%directly

const = 1./sqrt(kbt*4.184) * 6.62606 * 6.0221 * 0.1 / 2 / pi;

pe=1./sqrt(kbt*ev*4.184) * 6.62606 * 6.0221 * 0.1 / 2 / pi;

squant= kb2* sum (pe./(exp(pe) - 1) - log (1-exp(-pe)))  %entropy of PPS ; this is also the same as the quasiharmonic entropy in numata & knapp
sclas = -kb2 * sum(log(pe)+1)

tsquant=T*squant
tsclas=T*sclas

% to joules : 
%ps*4184



