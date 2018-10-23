function pdir = mpdir(p,direction)
%MPDIR Directional derivative of a polynomial in (shifted normalized) power form
%
%   MPDIR(P,DIRECTION)  returns the (shifted power form of the) derivative in 
%   the given DIRECTION(s) of the function contained in F (and of the same 
%   order).  
%   If the polynomial in P is d-variate, then DIRECTION must be a (list of)
%   d-vector(s), i.e., of size [d,nd] for some nd.  
%
%   Assuming the polynomial in P to be d-variate and (m-vector)-valued, PDIR
%   describes the (m*nd)-valued polynomial whose value at a point, reshaped
%   as a matrix of size [m,nd], provides in its j-th column the directional
%   derivative, of the polynomial in P, in the direction DIRECTION(:,j), j=1:nd.
%
%   For example, if  p  describes a d-variate m-vector-valued polynomial and
%    x  is some point in its domain, then
%
%      reshape(mpval(mpdir(p,eye(d)),x),m,d)
%
%   is the Jacobian of that function at that point.
%
%   As another example, a contour plot with gradients superimposed is obtained
%   by the following statements:
%     x = -1:.1:1; [xx,yy] = meshgrid(x,x); z = zeros(size(xx));
%     z(:) = mpval(p,[xx(:) yy(:)].');
%     contour(x,x,z), hold on, axis equal
%     dp = mpdir(p,eye(2));
%     grads = reshape(mpval(dp,[xx(:) yy(:)].'),[2,length(xx),length(yy)]);
%     quiver(x,x,squeeze(grads(1,:,:)),squeeze(grads(2,:,:))), hold off

%  cb: 23jan99, 07feb99

[dd,nd] = size(direction);

% if p is in plain power form, convert it first to normalized form:
if iscell(p)&p{1}==95, p = mp2fm(p); end

[coefs,d,k,center,blaise,mm] = mpbrk(p);

[df,nc] = size(coefs);

if dd~=d
   error(['The given polynomial is ',int2str(d),'-variate', ...
            ', while DIRECTION is in R^',int2str(dd),'.'])
end

%  Let  p(x) = sum_a c(a) x^a (|a| choose a) . Then
%   D_y p(x) = sum_a c_y(a) x^a (|a| choose a) , with
%
%   c_y(a) = (|a|+1) y'* (c(a+e_i): i=1:d), all  a .

if k==0  % the polynomial is a constant, hence its derivative is zero
   pdir = mpmak(zeros(df*nd,1),d,k,center,blaise,mm);
   return
end

dcoefs = zeros(df*nd,blaise(end,end-1));
for j=1:k
   rangej = blaise(end,j)+1:blaise(end,j+1); l = length(rangej);
   % we want sum_{j=1:d} coefs(:,j,:)*direction(j,:), with
   %  coefs(:,:,:) = reshape(coefs(:,blaise(end,j+1)+mm(1:d,rangej)),[df,d,l])

   dcoefs(:,rangej) = j*permute(...
       reshape(reshape( ...
       permute(reshape(coefs(:,blaise(end,j+1)+mm(1:d,rangej)),[df,d,l]),[3 1 2]),...
       [l*df,d])*direction, [l,df*nd]),[2 1]);

%   dcoefs(:,rangej) = j*reshape(...
%         sum(...
%         reshape(repmat(coefs(:,blaise(end,j+1)+mm(1:d,rangej)),nd,1),...
%                                        [df*nd,d,length(rangej)]).* ...
%         reshape(repmat(direction(:).',df,1),[df*nd,d,length(rangej)]),...
%             2), [df*nd, rangej]);
end

% Since the derivative is of one degree less, we drop the last entry of
% dimpjd, and also modify the last part of mm, mindful of the fact that
% this means dropping mm(1:d,rangek) (with rangek the highest range used
% in the preceding loop). It also means changes in the next range,
% rangekp1 = blaise(end,end-1)+1:blaise(end,end).
% All the rows of  mm(1:d,rangekp1)  contain the most recent m, whose last
% part, namely the range 
topend = (blaise(end,end)-blaise(end,end-1)) +rangej-blaise(end,end-2);
% contains the previous m.
% Altogether, this means removing from  mm  the following columns:
mm(:,blaise(end,end-2)+[1:blaise(end,end)-blaise(end,end-1)]) = [];

pdir = mpmak(dcoefs,d,k-1,center,...
       blaise(:,1:end-1),mm(:,[1:blaise(end,end-2) end-length(rangej)+1:end]));
