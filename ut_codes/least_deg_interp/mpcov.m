function mpnew = mpcov(mp,A)
%MPCOV linear change of variables
%
%        MPNEW = MPCOV(MP,A)
%
% Returns the mp-form  of  y\mapsto p(Ay)  for the  p  in  MP
%
% See MPMAK for details concerning the mp-form 
%   p(x) = sum_a (x-center)^a (|a| choose a) coefs(:,nn(a))
%
% At present, the brute-force method of interpolation is used.

% cb: 25jan99,08feb99
%   Copyright (c) 1999 by C. de Boor

[coefs,d,k,center,blaise,mm] = mpbrk(mp);
[df,nc] = size(coefs);
[ra,newd] = size(A);
if ra~=d, error(['The matrix needs to have the same number of rows ',...
          'as the polynomial has variables.'])
end

newcenter = zeros(newd,1);
if any(center), newcenter = A\center; end

pts = zeros(d,1);

for j=1:k

   % we need A_j(d) = {|a|=j}, in lexicographical order but with its elements
   % stored in the columns of a d-rowed matrix, and obtain it from
   % its last row and the permutation that orders that row, as stored
   % in the last two rows of mm:
   Ajd = zeros(d,blaise(d,j+2)); 
   Ajd(end,:) = mm(d+1,blaise(end,j+1)+1:blaise(end,j+2));
   pj = mm(d+2,blaise(end,j+1)+1:blaise(end,j+2));
   for jj=d-1:-1:1
      Ajd(jj,pj) = Ajd(jj+1,:);
   end
   pts = [pts Ajd];

end

mpnew = mpapi(pts,mpval(mp,A*pts),1e-10,newcenter);

