% read & check coords
format long ;
%r=load('../positions0.dat');r=reshape(r,3,[]); % no difference when positions are saved in the main wshell routine - check plugin master ?
r=load('../positions.dat');r=reshape(r,3,[]); % fortran indexing (w/o transpose)
natom=numel(r)/3;

surf_inds=load('../surf_inds.dat');
solv_inds=load('../solv_inds.dat');

r_surf=load('../r_surf.dat');r_surf=reshape(r_surf,3,[]);
r_solv=load('../r_solv.dat');r_solv=reshape(r_solv,3,[]);

nsurf=numel(r_surf)/3;
nsolv=numel(r_solv)/3;

fprintf(':%d total atoms\n',natom)
fprintf(':%d surface atoms\n',nsurf)
fprintf(':%d solvent atoms\n',nsolv)

fprintf(':checking rsurf against inst positions:') ;
 max ( max( abs ( r_surf-r(:,surf_inds) )) ) % do _not_ match exactly -- why -- r(:,5) contaminated ... even at beginning of wshell_main !

fprintf(':checking rsolv against inst positions:')
 max ( max( abs ( r_solv-r(:,solv_inds) )) )

