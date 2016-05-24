% plot parts of ftsm string
% need to strip off formating, e.g. as in below :
%#!/bin/bash
%# remove extraneous text from string file so we can read it into matlab
%sfile='fip2rw-o_ftsm-zts.txt'; # initial path from zts
%sfile='fip2rw-o_ftsm.txt'; # path after optimization on odysey
%
%temp='string.dat';
%grep "^\sr[f,o]" $sfile | awk '{print $2" "$3" "$4}' > $temp

%

nrep=32; % string replicas


sfile='string.dat'

d=load(sfile);


[nrow,dim]=size(d);

npos=nrow/nrep;
x=d(:,1) ;
y=d(:,2) ;
z=d(:,3) ;

x=reshape(x,npos,nrep);
y=reshape(y,npos,nrep);
z=reshape(z,npos,nrep);

%now can plot paths of positions :

range=16:31 ;
clf

for i=range ;
 xx=x(i,:);
 yy=y(i,:);
 zz=z(i,:);
 plot3(xx,yy,zz,'r.-') ; hold on

end