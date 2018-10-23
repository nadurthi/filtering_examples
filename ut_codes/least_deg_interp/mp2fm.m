function mpnew = mp2fm(mp)
%MP2FM converts between the two polynomial forms
%
%        MPNEW = MP2FM(MP)
%
% Converts between the two forms (normalized and plain)
%
% See MPMAK for details concerning the normalized mp-form 
%   p(x) = sum_a (x-center)^a (|a| choose a) coefs(:,nn(a))
% and the plain mp-form
%   p(x) = sum_a (x-center)^a coefs(:,nn(a))

% cb: 25jan99,08feb99
%   Copyright (c) 1999 by C. de Boor

[coefs,d,k,center,blaise,mm] = mpbrk(mp);
df = size(coefs,1);

% Need to generate the relevant multinomial coefficients (|a| choose a).
% Do it with the aid of mm, using the fact that
% mcjp1 := ( (|a| choose a) : |a|=j+1 )  can be obtained from  mcj  by the loop
%
%   mcjp1 = zeros(1,dimAjp1);
%   for i=1:d
%      mcjp1(mm(i,:)) = mcjp1(mm(i,:)) + mcj;
%   end

mcs = zeros(1,blaise(end)); mcs(1) = 1;
for j=1:k
  range = blaise(end,j)+1:blaise(end,j+1); mcjm1 = mcs(range); 
  for i=1:d
     temp = blaise(end,j+1)+mm(i,range);
     mcs(temp) = mcs(temp) + mcjm1;
  end
end  

switch mp{1}
case 94
   mpnew = mpmak(coefs.*repmat(mcs,df,1),d,k,center,blaise,mm);
   mpnew{1}=95;
case 95
   mpnew = mpmak(coefs./repmat(mcs,df,1),d,k,center,blaise,mm);
otherwise
   error('No conversion from given form available at present.')
end
