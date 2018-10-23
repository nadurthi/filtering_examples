function mpnew = mpshft(mp,newcenter)
%MPSHIFT Shift to a different center
%
%        MPNEW = MPSHFT(MP,NEWCENTER)
%
% Recenters the mp-form at the given NEWCENTER. 
%
% See MPMAK for details concerning the mp-form 
%   p(x) = sum_a (x-center)^a (|a| choose a) coefs(:,nn(a))
%

% cb: 01aug97, 9aug97
%   Copyright (c) 1997 by C. de Boor

[coefs,d,k,center,blaise,mm] = mpbrk(mp);
[df,nc] = size(coefs);

if size(newcenter)~=[d,1]
   error(['NEWCENTER should be of size [',int2str(d),',1].'])
end
x = reshape(newcenter-center,[1,d,1]);

for i=0:k-1
   for j=k-1:-1:i 
      ranj = blaise(end,j+1)+1:blaise(end,j+2); lj = length(ranj);
      coefs(:,ranj) = coefs(:,ranj) + ...
        reshape(sum(reshape(coefs(:,blaise(end,j+2)+mm(1:d,ranj)),[df,d,lj]) ...
                           .*x(ones(df,1),:,ones(1,lj)),2),[df,lj]);
   end
end

mpnew = mpmak(coefs,d,k,newcenter,blaise,mm);
