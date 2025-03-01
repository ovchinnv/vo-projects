% combine multiple force files from serially-run replicas into one file
%

brep=0;
erep=15;

nsamples=1300 ; % number of force samples to take from file (useful if files do not have an equal number of samples)

basename='force';
ext='.dat';

for i=brep:erep
 fname=[basename, num2str(i),ext];
 if (i==brep)
  dnew=load(fname);
  if exist('nsamples')
   d=dnew(1:nsamples,:);
  else
   d=dnew;
  end
 else
  dnew=load(fname);
  if exist('nsamples')
   d=[d dnew(1:nsamples,:)];
  else
   d=[d dnew];
  end
 end
end

% write one file
fname=[basename,ext];
fout=fopen(fname,'w');

nlines=size(d,1);
isamp=0;
for i=1:2:nlines
 isamp=isamp+1;
 fprintf(fout,'%% %10d\n',isamp) ;
 fprintf(fout,[repmat('%15.10f ',[1,erep-brep+1]),'\n'],d(i,:));
 fprintf(fout,[repmat('%15.10f ',[1,erep-brep+1]),'\n'],d(i+1,:));

end

fclose(fout)