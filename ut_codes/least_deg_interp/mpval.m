function v = mpval(p,x)
%MPVAL evaluate multivariate polynomial in (shifted, normalized) power form
%
%        V = MPVAL(P,X)
%
% Returns the value at  X  of the (d-variate) polynomial whose (shifted,
% normalized) power form 
%
%   p(x) = sum_a (x-center)^a (|a| choose a) coefs(:,nn(a))
%
% is contained in  P , with  X  of size  d-by-nx .
%
% The powers are so normalized that a multivariate version of Horner's scheme
% (aka the de Casteljau algorithm) works simply; see below.
%
% See MPMAK for details on the form.

% cb: 26jul97, 9aug97
%   Copyright (c) 1997 by C. de Boor

[dx,nx] = size(x);
[coefs,d,k,center,blaise,mm] = mpbrk(p);
[df,nc] = size(coefs);

if dx~=d 
   error(['The given polynomial is ',int2str(d),'-variate', ...
            ', while  X  is in R^',int2str(dx),'.'])
end

x = x - repmat(center,1,nx);

switch mpbrk(p,'form')

case 94  % normalized shifted power form, calls for
%   Symmetric Horner's scheme (aka Nested Multiplication) consists of the 
%   following loop:
%
%     vk <-- (coefs(nn(a)) : |a|=k)
%     for j=k-1:-1:0
%        vj <-- (coefs(nn(a)) : |a|=j) + sum_{i=1:d} x(i)*vjp1(n(a+e_i))
%     end
%
%   Hence  v0 = sum_a coefs(nn(a)) x^a * # sequences (0=b_0,b_1,...,b_|a|=a) 
%   with  b_{j+1}-b_j in {e_1,...,e_d}, all j . In other words,
%       v0 = sum_a x^a (|a| choose a) coefs(nn(a)) = p(x) .
%
%   For this, it is important to know that
%   BLAISE(end,2+j) = dim Pi_j(R^d) = ( j + d choose d ) is the number of 
%   d-variate monomials of degree  <=j , while   
%   M(i,nn(a)) = n(a+e_i), i=1:d , and
%   nn(a) = dim Pi_{<|a|} + n(a) .

   v = coefs(:,blaise(end,k+1)+1:blaise(end,k+2),ones(1,nx));
   if k>0
      x = reshape(x,[1,d,1,nx]);
      for j=k-1:-1:0
         rangej = blaise(end,j+1)+1:blaise(end,j+2); lj = length(rangej);
         v = coefs(:,rangej,ones(1,nx)) + ...
                 reshape(sum(reshape(v(:,mm(1:d,rangej),:),[df,d,lj,nx]).* ...
                     x(ones(df,1),:,ones(lj,1),:),2),[df,lj,nx]);                
      end
   end

case 95 % plain shifted power form, calls for

% plain (unsymmetric) Horner's scheme (aka Nested Multiplication)
%     vk <-- (coefs(nn(a)) : |a|=k)
%     for j=k-1:-1:0
%        vj <-- (coefs(nn(a)) : |a|=j)
%        for i=1:d
%           vj(1:#A_j(i)) = vj(1:#A_j(i)) + ...
%                             x(d+1-i)*vjp1(#A_{j+1}(i-1)+1:#A_{j+1}(i))
%        end
%     end
%
% Here, #A_j(i) is the number of monomials of exact degree j  in  i  arguments,
% i.e., #A_j(i) = BLAISE(i+1,j+2)
   v = coefs(:,blaise(end,k+1)+1:blaise(end,k+2),ones(1,nx));
   if k>0

     x = reshape(x,[1,d,nx]);
     for j=k-1:-1:0
        temp = coefs(:,blaise(end,j+1)+1:blaise(end,j+2), ones(1,nx));
        blaises = [0;blaise(:,j+3)];
        for i=d:-1:1
           rangei = blaises(i)+1:blaises(i+1); li = length(rangei);
           temp(:,1:li,:) = temp(:,1:li,:) + ...
              x(ones(df,1),repmat(i,1,li),:).*v(:,rangei,:);
        end
        v = temp;
     end
   end

end
v = reshape(v,[df,nx]);
