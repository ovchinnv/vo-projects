% combine multiple force files from serially-run replicas into one file
%

brep=0;
erep=15;

basename='force';
ext='.dat';

for i=brep:erep
 fname=[basename, num2str(i),ext];
 if (i==brep)
  d=load(fname);
 else
  d=[d load(fname)];
 end
end

% write one file
fname=[basename,ext];
fout=fopen(fname,'w');

nlines=size(d,1);
isamp=0
for i=1:2:nlines
 isamp=isamp+1;
 fprintf(fout,'%% %10d\n',isamp) ;
 fprintf(fout,[repmat('%15.10f ',[1,erep-brep+1]),'\n'],d(i,:));
 fprintf(fout,[repmat('%15.10f ',[1,erep-brep+1]),'\n'],d(i+1,:));

end

fclose(fout)