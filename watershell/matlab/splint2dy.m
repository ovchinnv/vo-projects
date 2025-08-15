% integrate analytically a 2D spline in the y- (second) dimension
function ppi2dy=splint2dy(pp2d)
% integrate `constituent' splines in pp form segment-by-segment
% initialize 1D spline:
 ppy.form = pp2d.form;
 ppy.order = pp2d.order{2};
 ppy.breaks = pp2d.breaks{2};
 ppy.pieces = pp2d.pieces{2};
 ppy.dim = 1 ;
% initialize new coeffs :
 oldcoefs=pp2d.coefs;
 newcoefs=zeros(pp2d.order{1}, pp2d.order{2}+1, pp2d.pieces{1}, pp2d.pieces{2});
 for i=1:pp2d.order{1};
  for j=1:pp2d.pieces{1};
   ppy.coefs=squeeze(pp2d.coefs(i,:,j,:))'; % note that the 1D spline uses the transpose of the coeff matrix, e.g. (nx-1,1:order)
%   ppy.coefs=pp2d.coefs(i,:,j,:)';
   ppiy=splint(ppy);
   newcoefs(i,:,j,:)=ppiy.coefs';
  end
 end
% now, populate 2D integral spline :
 ppi2dy=pp2d; % copy everything from lower order spline
 ppi2dy.order{2}=ppi2dy.order{2}+1; % increase order in y
 ppi2dy.coefs=newcoefs; % replace coefs
end
