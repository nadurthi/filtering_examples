function [mk,blaises,mp1] = mpnext(blaises,mp)
%MPNEXT next section of MM
%
%        [MK,BLAISES,MP1] = MPNEXT(BLAISES,MP)
%
% input:  BLAISES := ( (k+h-1 choose h-1) : h=1:d ), hence k+1 = BLAISES(2);
%        MP(1,:) := the last column of the matrix containing in its rows the
%                     elements of A_k(d) := { a in Z_+^d : |a| = k } in 
%                     lexicographic order;
%        MP(2,:) := the result of  [ignored, MP(2,:)] = SORT(MP(1,:)) ;
% output:  MK  for which  MK(i,n(a)) = n(a+e_i), i=1:d , with  n(a)  the place
%             at which  a  occurs in the lexicographic ordering of  A_{|a|}(d),
%             while MK(d+1:d+2,:) equals MP;
%          as well as  BLAISES  and  MP  for  k+1  rather than  k .
%       
% This is of use in working with multivariate polynomials in power form.
%
%   See also: MPMAK, MPAPI

% cb: 9aug97; 24jan99 (made BLAISES a d-vector only); 
%     7feb99 (combined m and p, saved them in mm)
%   Copyright (c) 1997 by C. de Boor

if nargin==1  % initialize the process
   d = blaises(1);
   blaises = ones(d+1,1);
   mk = zeros(d+2,0);
   mp1 = [0;1];
else          % update it
   d = length(blaises)-1;
   % Ready to generate, in MP1, the last column of the matrix which has
   % in its rows the elements of A_{k+1}(d), in lexicographic order. This
   % column starts off with  k+1 , followed by increasingly longer initial
   % segments of the last column of the same matrix for  A_k(d) , presumably
   % available in  M .
   m1 = mp(1,1)+1;
   for jj=2:d
      m1 =[m1 mp(1,1:blaises(jj))];
   end
   [ignored,p1] = sort(m1);  %  gets the permutation that orders  M1 
   nakd = blaises(d);    % = # A_k(d) 
   blaises = cumsum(blaises); % raises BLAISES to the next level
   rangek = (blaises(d)-nakd) + [1:nakd]; % the tail-end of interest
   % Ready to generate  MK  for which MK(i,n(a)) = n(a+e_i), i=1:d, a in A_k(d).
   % This is done as follows. It is easily seen that the map  a |--> a+e_1
   % carries  A_k(d)  without any order change to the tail-end of the ordered
   %  A_{k+1}(d) , hence  MK(1,:) = RANGEK := (BLAISES(d)-NAKD) + [1:NAKD] .
   % For  a+e_j , let  Ak  be the matrix whose rows contain the elements of
   %  A_k(d)  in lexicographic order. Then, in particular,  Ak(:,d) = M
   % and  Ak(P,[d 1:d-1]) = Ak  (since  P  is the permutation that puts  
   %  M  into increasing order without unnecessary interchanges). Hence  
   %        Ak(P^{-1},[2:d 1]) = Ak  and  Akp1(P1^{-1},[2:d,1]) = Akp1 .
   % Therefore,
   %          MK(2,P^{-1}) = P1^{-1}(RANGEK) ,
   % and, in general
   %    MK(i,P^{-i+1}) = P1^{-i+1}(RANGEK), i=1:d
   % and, of course,  P^{-d} = id, P1^{-d} = id .
   % This is realized below by starting with  NEXT = 1:BLAISES(d), and
   % repeatedly permuting the columns of  MK  by  P^{-1}  and the columns of
   % NEXT, in concert, by P1^{-1}, and then setting the appropriate row of
   % the permuted  MK  equal to the tail-end of the analogously permuted  NEXT .
   
   mm = zeros(d,nakd); next = 1:blaises(d);
   for i=1:d
      mk(i,:) = next(rangek);
      mk(:,mp(2,:)) = mk; next(p1) = next;
   end
   mk = [mk;mp]; mp1 = [m1;p1];
end

