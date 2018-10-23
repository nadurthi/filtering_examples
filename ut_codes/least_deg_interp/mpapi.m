function [mp,mpfactor] = mpapi(t,f,tol,center)
%MPAPI multivariate polynomial interpolation
%
%        [MP,MPFACTOR] = MPAPI(T,F[,TOL[,CENTER]])
%
% returns the least polynomial interpolant to the data  
%
%          (T(:,i),F(:,i)), i=1,...,n,
%
% and, optionally, the cell array MPFACTOR containing all the information
% generated during factorization for the construction of an interpolant
% at the same data sites to given data F, using the subsequent call
%
%        MP = MPAPI(MPFACTOR,F) 
%
% The program is a close copy of the informal program in [de Boor and Ron,
% Computational aspects of multivariate polynomial interpolation, Math.Comp.
% 58 (198); 1992; 705--727], but is explicit about the linear ordering of the
% sets A_jd := {a in Z_+^d : norm(a)_1 = j} and the resulting computation of
% each block of the blocked Vandermonde from its preceding block. The same
% linear ordering is used in MPVAL to evaluate the interpolating polynomial
% by Horner's method. Also, the coefficients are returned in that ordering.

% cb: 28jul97, 9aug97, 24jan99 (made BLAISES a d-vector only)
%   Copyright (c) 1997 by C. de Boor

if iscell(t),   % no need to repeat the factorization, just extract all the
                % needed information from the cell array T = MPFACTOR
                % generated earlier by the statement 
                %        mpfactor = {center LU p K blaise mm gstring};
   center = t{1}; d = length(center);
   LU = t{2}; p = t{3}; n = length(p); K = t{4}; k = K(end);
   blaise = t{5}; mm = t{6}; gstring = t{7};
else,           % need to carry out factorization
   [d,n] = size(t);
   
   if nargin<4, center = sum(t,2)/n; end
   % subtract the CENTER from the data sites ...
   nones = ones(n,1);
   t = t - center(:,nones);
   % and make the point closest to the origin the first point to be considered
   % (to avoid a zero row later in the interesting part of  Vk  for  k>0 )
   % by setting  P(1)  equal to its index, where  P  is the n-sequence
   % keeping track of pivoting. Precisely, at the end of the j-th elimination 
   % step, row P(i) is known to be pivot row for step i=1:j, while rows P(i), 
   % i=j+1:n, are not yet used as pivot rows.  
   p = 1:n; [ignored,i] = min(sum(abs(t),1)); p(1)=i; p(i)=1;

   % Use the mechanism for generating  Vkp1  from  Vk  to generate the weight
   % vector for the scalar product, by adjoining the point  e := ones(d,1) :
   t = [ones(d,1) t];
   
   if nargin<3, tol = 1.e-14; end % Elimination needs a tolerance for
                       % recognizing when a prospective pivot element is zero.
   
   % initialize Gauss elimination items: 
   LU = zeros(n,n);  % to contain the interesting parts of  L  and  U ;
   K = zeros(1,n);   % to record the pivot degree for each unknown;
                     % the vector  P  for recording interchanges was set earlier
   % initialize block Vandermonde:
   % Keep around only the current block of the Vandermonde, but accumulate the 
   % the various  (g_j)_least, in GSTRING.
   Vk = ones(n+1,1); np=2:n+1; Wk = Vk(np,:); gstring = [];
   % also compute the initial size of each entry of Vk:
   iVk = Vk(np,:).^2./Vk(nones,:);
   % initialize indexing scheme for A_jd = {a in Z_+^d : |a| = j}
   blaise = zeros(d+1,2);
   [mm,blaise(:,2),mpm] = mpnext(d);
   
   k = 0; kfact = 1;
   for j=1:n
      while 1 % determine the next pivot row, generating the next block of
              % the Vandermonde if necessary.
         [pivot,ipivot] = ...
            max(sum(Wk(p(j:n),:).^2./Vk(ones(n+1-j,1),:),2)./iVk(p(j:n)));
         % [pivot,ipivot] = max(sum(abs(Wk(p(j:n),:)),2)./iVk(p(j:n)));
         if pivot>tol, break, end
         k = k+1; kfact= kfact*k;
         % generate the next block
         % first, get next piece of  MM  and
         % update the other combinatorial information, specifically
         nakd = blaise(d,end);              % = (k-1 + d-1 choose d-1)
         [mk,blaise(:,k+2),mpm] = mpnext(blaise(:,k+1),mpm);
         mm = [mm mk];
         %  now  BLAISE(:,end)  = ( (k + h-1 choose h-1) : h=1:d+1 )
         %  and  BLAISE(end,:) = ( dim Pi_h(R^d) : h=-1:k )
         %  Initialize the next V-block to be zero ...
         Vkp1 = zeros(n+1,blaise(d,end));
         %  ... then add to it a properly ordered t(i,:)*Vk, i=1:d
         for i=1:d
            Vkp1(:,mk(i,:)) = Vkp1(:,mk(i,:)) + ...
                                 reshape(t(i,:,ones(1,nakd)),n+1,nakd).*Vk;
         end
         Vk = Vkp1; 
         iVk = sum(Vk(np,:).^2./Vk(nones,:),2);
         % iVk = sum(abs(Vk(np,:)),2);
         Wk = Vk(np,:);
         % apply the previous elimination steps to the new block:
         % for jj=1:j-1
         %    Wk(p(jj+1:n),:) = Wk(p(jj+1:n),:) - LU(jj+1:n,jj)*Wk(p(jj),:);
         % end
         Wk(p,:) = (tril(LU,-1)+eye(n,n))\Wk(p,:);
      end
   
      if ipivot>1 % exchange the relevant row indices:
         p(j-1+[ipivot,1]) = p(j-1+[1,ipivot]);
         LU(j-1+[ipivot,1],:) = LU(j-1+[1,ipivot],:);
      end
         
      K(j) = k;  % record the degree of g_jleast
   
      i=1:n;  (Wk(p(j),:)./Vk(1,:))/kfact; % = coefficients of g_jleast
      LU(i,j) = sum(Wk(p(i),:).*ans(nones,:),2);
      i=j+1:n;
      LU(i,j) = LU(i,j)/LU(j,j);
      Wk(p(i),:) = Wk(p(i),:) - LU(i,j)*Wk(p(j),:);
      
      gstring = [gstring ans];
   end
   mm(d+1:d+2,blaise(end,end-1)+1:blaise(end,end))= mpm; 
                % the most recent M is stored in the last part of MM, for
                % possible use, in degree-raising, integration, etc.
   

   if nargout>1    % store the information obtained during factorization
                   % in the cell array MPFACTOR for later use.
      mpfactor = {center LU p K blaise mm gstring n d k};
   end
end

[df,nf] = size(f);
if nf~=n
   error('The number of data values does not equal the number of data sites.')
end

% for j=1:n-1
%    i=j+1:n;
%    a(i) = a(i) - LU(i,j)*a(j);
% end

% a(n) = a(n)/LU(n,n);
% for i=n-1:-1:1
%   a(i) = (a(i) - LU(i,i+1:n)*a(i+1:n))/LU(i,i);
% end

a = (triu(LU)\((tril(LU,-1)+eye(n,n))\f(:,p).')).';

%  a = diag(LU).*a; % having suppressed division of Wk(p(j),:) by LU(j,j)
                    % earlier, no need now to multiply  a  by  diag(LU) .

coefs = zeros(df,blaise(end)); gend = 0;
for j=1:n
   rangej = blaise(end,1+K(j))+1:blaise(end,2+K(j)); 
   gbeg = gend+1; gend = gend + length(rangej);
   coefs(:,rangej) = coefs(:,rangej) + a(:,j)*gstring(gbeg:gend);
end

mp = mpmak(coefs,d,k,center,blaise,mm);
