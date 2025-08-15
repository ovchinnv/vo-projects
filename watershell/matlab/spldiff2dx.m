% differentiate analytically a 2D spline in the x- (first) dimension
function pp2dx=spldiff2dx(pp2d)
% implicitly differentiate `constituent' 1D splines
 diffcoefs=zeros(pp2d.order{1}-1, pp2d.order{2}, pp2d.pieces{1}, pp2d.pieces{2});
 orderx=pp2d.order{1};
 for i=1:pp2d.order{2};
  for j=1:pp2d.pieces{2};
%   coefs=squeeze(pp2d.coefs(1:orderx-1,i,:,j));
   coefs=(pp2d.coefs(1:orderx-1,i,:,j));
   newcoefs = bsxfun( @times, coefs , [orderx-1:-1:1]' );
   diffcoefs(:,i,:,j)=newcoefs;
  end
 end
% now, populate 2D derivative spline :
 pp2dx=pp2d; % copy everything from lower order spline
 pp2dx.order{1}=pp2dx.order{1}-1; % decrease order in x
 pp2dx.coefs=diffcoefs; % replace coefs
end
