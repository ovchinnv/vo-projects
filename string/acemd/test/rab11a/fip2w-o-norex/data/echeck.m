% compute fts restraint energy from projection variables
% and compare to restraint energy output from sm plugin
%

ibeg=0;
iend=31;

pbase='proj';
ebase='ftse';

% distance units of D -- need to convert to distance units of 2D

kpar_=1.63 ;
kprp_=1.63 ;
dperp_=7;

de=[];
for i = 1 : 31
% get projections :
 fp=[pbase,num2str(i),'.dat'];
 proj=load(fp);
% take the last 6 entries, which correspond to parallel and perpendicular positions
% of this replica, left neighbor replica and right neighbor replica
 proj=proj(end-5:end);
 pthis=proj(1:2);
 pleft=proj(3:4);
 pright=proj(5:6);

% get plugin energy
 ep=[ebase,num2str(i),'.dat'];
 ene=load(ep);
 ethis_=ene(1:2);
 eleft_=ene(3:4);
 eright_=ene(5:6);
%
% compare energies :
 [ethis(1), ethis(2)]=ftse(i,ibeg,iend,kpar_,kprp_,pthis(1),pthis(2),dperp_) ;
 dethis = ethis(:) - ethis_(:);
 if (i>ibeg)
  [eleft(1), eleft(2)]=ftse(i-1,ibeg,iend,kpar_,kprp_,pleft(1),pleft(2),dperp_) ;
  deleft = eleft(:) - eleft_(:);
 end
 if (i<iend)
  [eright(1), eright(2)]=ftse(i+1,ibeg,iend,kpar_,kprp_,pright(1),pright(2),dperp_) ;
  deright = eright(:) - eright_(:);
 end
%
 de(i) = sum ( abs(dethis) + abs(deleft) +abs(deright) ) ; 
%
%return
end

['maximum error:', num2str(max(de))]

