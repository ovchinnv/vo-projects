% compute free energy from forces acting on the CV
% assuming 1D
%
if exist('OCTAVE_VERSION')
 graphics_toolkit('gnuplot')
end
%
%close all;

nrun=5;
cv     =load('deca_tamd_cv1.dat');
fcv    =load('deca_tamd_cv_force1.dat');

for i=2:nrun
 cv=[cv ; load(['deca_tamd_cv',num2str(i),'.dat'])];
 fcv=[fcv ; load(['deca_tamd_cv_force',num2str(i),'.dat'])];
end

minlen =min(length(cv), length(fcv));
data   = [ cv(1:minlen), fcv(1:minlen) ] ;

x0=12 ;
x1=32 ;
nxcor=20 ;

xcor=linspace(x0, x1, nxcor);
xcen=0.5*(xcor(2:end)+xcor(1:end-1));
nxcen=nxcor-1;
dx=xcor(2)-xcor(1);

% bin cvs and gradients
% determine bin number
data(:,3)= floor ( ( data(:,1)-x0 )/ dx ) + 1;
data=sortrows(data,3);
% compute average force in each bin
fx =zeros(1,nxcen);
fxe=zeros(1,nxcen);
for ibin=1:nxcen
 inds=(find(data(:,3)==ibin));
 dd=data(inds,2);
 count(ibin)=numel(dd);
 fx(ibin)=mean(dd);
 if (~isnan(fx(ibin)))
  fxe(ibin)=std(dd)/sqrt(count(ibin));
 end
end

% integrate derivative :
fe=cumsum(fx)*dx ; 
fe=fe-fe(1);
%fe=fe-min(fe);
plot(xcen(1:nxcen), fe,'k*-'); hold on ;
errorbar(xcen(1:nxcen), fe, fxe,'k.-');
xlim([x0 x1])
ylim([-5 26])
hold on ;

%%%%%%%%%%%% plot US data :
fnames={'usfb/fe8.mat','usfb/fe10.mat','usfb/fe16.mat','usfb/fe32.mat'};
styles={'ro-','gs-','bv-','k.-'};
i=0;
inds=[4];
for fname=fnames(inds)
 load(char(fname));
 i=i+1;
% plot(cvs0,fave-min(fave), char(styles(i)));
 plot(cvs0,fave, char(styles(i)), 'markersize',8,'linewidth',0.75);
end
errorbar(cvs0,fave,fstd,char(styles(i)));
usleg={'$\delta$=2.5806','$\delta$=1.9355','$\delta$=1.2903','$\delta$=0.6452'};
legend([{'TAMD 6000K'} usleg(inds)], 'location','southeast');

%legend(leg,'location','northwest');
box on;
ylabel('\it F (kcal/mol)', 'fontsize',14);
xlabel('\it x', 'fontsize',14);
%
set(gcf, 'paperpositionmode', 'auto');
print(gcf, '-depsc2', 'tfe2.eps');
print(gcf, '-dpng', 'tfe2.png');
%