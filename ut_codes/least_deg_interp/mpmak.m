function mp = mpmak(coefs,d,k,center,blaise,mm)
%MPMAK multivariate polynomial in (shifted normalized) power form
%
%        MP = MPMAK(COEFS,D[,K,CENTER,BLAISE,MM])
%
% Returns in MP the given (and other helpful) information about the polynomial
%
%   p(x) = sum_a (x-center)^a (|a| choose a) coefs(nn(a)) ,
%
% with  nn(a) := dim Pi_{<|a|}(R^d) + n(a) ,  and  n(a)  the position of  a  in 
% the lexicographic ordering of  A_{|a|}(d) , where  |a| := sum_{i=1:d} a(i) 
% and  A_j(d) := {a in Z_+^d : |a| = j} .
% The combinatorial factor is just right to make the evaluation of  p  by 
% Horner's symmetric scheme very simple (see MPVAL).
%
% Further, BLAISE(i,j) = #A_{j-2}(i-1), hence, in particular,
%  BLAISE(d+1,:) = (dim Pi_j(R^d): j=-1:k) ,  and  MM  contains indexing 
% information related to the graded lexicographic ordering. Precisely,
%
%    MM(i,nn(a)) = n(a+e_i), i=1:d, 
%
% while MM(d+1,BLAISE(d+1,j+1)+1:BLAISE(d+1,j+2)) is the last column of the 
% matrix that contains the elements of A_j(d) in lexicographic order in its 
% rows, and MM(d+2,BLAISE(d+1,j+1)+1:BLAISE(d+1,j+2)) is the permutation that 
% sorts that column.  

% todo: add row MM(d+3,:) so that nn(a+e_i) = MM(i,nn(a))+ MM(d+3,nn(a))
% cb 26jul97; 24jan99 (made BLAISE matrix to carry both former blaises and
%                      dimpjd)
%          07feb99 (combined m and mp and stored it in mm)
%   Copyright (c) 1997 by C. de Boor

[df,nc] = size(coefs);
if nargin<3 % we must deduce  K  from  NC 
   index = ones(1,d+1); k = 0;
   while index(d+1)<nc, index = cumsum(index); k = k+1; end
   if index(d+1)>nc
      mp = []; error('For every k, length(coefs(1,:)) neq dim Pi_k(R^d) .'), end
end

if nargin < 4, center = zeros(d,1); end

if nargin < 5 % we must generate  BLAISE and  MM 
   blaise = zeros(d+1,k+2);
   [mm,blaise(:,2),mpm] = mpnext(d);
   for j=1:k
      [mj,blaise(:,j+2),mpm] = mpnext(blaise(:,j+1),mpm);
      mm = [mm mj];
   end
   mm(d+1:d+2,blaise(end,end-1)+1:blaise(end))= mpm; 
                % the most recent M is stored in the last part of MM, for
                % possible use, in degree-raising, integration, etc.
end

% % of use on MATLAB versions prior to v5
% mp = [94, df, d, reshape(center,1,d), k, blaise(:), ...
%        reshape(coefs,1,df*nc), reshape(mm,1,d*nc)];

mp = {94, center, coefs, blaise, mm};
