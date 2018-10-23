function dmp = mpder(mp,xi)
%MPDER derivative of a multivariate polynomial in mp-form
%
%        DMP = MPDER(MP,XI)
%
% applies the directional derivatives D_XI(:,j) to the polynomial in MP.
%
%   See also: MPMAK

% cb: 9aug97, 07feb99
%   Copyright (c) 1997 by C. de Boor

% MP contains (see MPMAK), among other things, the COEFS for the representation
% 
%   p(x) = sum_a x^a (|a| choose a) COEFS(nn(a)) 
%
% Differentiation in the direction  xi  gives  
%   D_xi ()^a/a! = sum_i xi(i) ()^{a-e_i}/(a-e_i)!
% hence
%
%  D_xi p(x) = sum_a x^a (|a| choose a) (|a|+1) \sum_i xi(i)*COEFS(nn(a+e_i))

% if mp is in plain power form, convert it first to normalized form:
if iscell(mp)&mp{1}==95, mp = mp2fm(mp); end

[coefs,d,k,center,blaise,mm] = mpbrk(mp);

[dxi,nxi] = size(xi);
if dxi~=d
   error(['The direction(s) must be in R^',num2str(d),'.']), end

[df,nc] = size(coefs);
for ixi=1:nxi
   x = reshape(xi(:,ixi),[1,d,1]);
   if k==0
      coefs(:,1) = zeros(df,1);
   else
      for j=0:k-1
         ranj = blaise(end,1+j)+1:blaise(end,2+j); lj = length(ranj);
         coefs(:,ranj) = (j+1)*reshape(sum( ...
                  reshape(coefs(:,blaise(end,j+2)+mm(1:d,ranj)),[df,d,lj]) ...
                                 .*x(ones(df,1),:,ones(1,lj)),2),[df,lj]);
      end
      % adjust MM to the lowered degree
      mm(:,blaise(end,1+k)+1:blaise(end,2+k)-lj) = []; 
      k = k-1;
   end
end

dmp = mpmak(coefs(:,1:blaise(end,k+2)),d,k,center,blaise(:,1:k+2),mm);
