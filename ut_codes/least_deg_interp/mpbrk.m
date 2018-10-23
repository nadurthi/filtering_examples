function [out1,d,k,center,blaise,mm] = mpbrk(mp,part) 
%MPBRK parts of a multivariate polynomial in (shifted normalized) power form 
%
%        [COEFS,D,K,CENTER,BLAISE,M] = MPBRK(MP)
%
% breaks the mp-form in MP into its parts and returns as many of them as
% are specified by the output arguments.
%
%        OUT1 = MPBRK(MP, PART)
%
% returns the part, specified by PART, of the multivariate polynomial in 
% (shifted normalized) power form described in MP. 
%
% The string PART may be (the beginning character(s) of) one of the following 
% strings: 'coefs', 'dimension', 'k' or 'degree', 'center', 'blaise', 'mm',
%          'form'..
%
% If no output it specified, prints all the details instead.
%
%   See also: MPMAK

% cb: 26jul97, 10aug97, 01feb99 (add form), 07feb99 (adjust to change in mm)
%   Copyright (c) 1997 by C. de Boor

if iscell(mp)        % mp = {94, center, coefs, blaise, mm};
   formi = mp{1};
   if formi~=94&formi~=95
      error(['The input cell array does not seem to describe', ...
             ' a multivariate polynomial in shifted power form.'])
   end
   centeri = mp{2}; coefsi = mp{3}; blaisei = mp{4}; mmi = mp{5};
   di = length(centeri); [dfi,nci] = size(coefsi); ki = size(blaisei,2)-2;
else                 % dealing with pre-v5 MATLAB
   if mp(1)~=94
      error(['The input array does not seem to describe', ...
             ' a multivariate polynomial in shifted power form.'])
   end
   formi=94;
   dfi = mp(2); di = mp(3); centeri = mp(3+[1:di]).'; ki = mp(4+di);
   sofar = 4+di; kp2dp1 = (ki+2)*(di+1);
   blaisei = reshape(mp(sofar+[1:kp2dp1]),di+1,ki+2); sofar = sofar+kp2dp1;
   nci = blaisei(end);
   coefsi = reshape(mp(sofar+[1:(dfi*nci)]),dfi,nci); sofar = sofar + dfi*nci;
   mmi = reshape(mp(sofar+1:end),di+2,nci);
end

if nargin>1 % a specific part is called for
   if nargout>1, error('Too many output arguments for the given input.')
   else
      switch part(1)
      case 'c'
         if length(part)<2, error('Second argument is ambiguous'), end
         switch part(2)
         case 'o'
            out1 = coefsi; return
         case 'e' 
            out1 = centeri; return
         otherwise
            error(sprintf('''%s'' is not part of a mp-form.',part))
         end
      case 'd'
         if length(part)<2, error('Second argument is ambiguous'), end
         switch part(2)
         case 'i', out1 = di; return
         case 'e', out1 = ki; return
         otherwise
            error(sprintf('''%s'' is not part of a mp-form.',part))
         end
      case 'f', out1 = formi; return
      case 'k', out1 = ki; return
      case 'b', out1 = blaisei; return
      case 'm', out1 = mmi; return
      otherwise
         error(sprintf('''%s'' is not part of a mp-form.',part))
      end
   end
else
   if nargout==0 % print out all the parts
   kind = '(normalized)'; if mp{1}==95, kind = '(plain)'; end
      fprintf(['The input is a ',num2str(di),'-variate'])
      if dfi>1 fprintf([' R^',num2str(dfi),'-valued']), end
      fprintf([' polynomial of degree ',num2str(ki), ...
               ',\n centered at the point\n'])
      disp(centeri)
      % coefsi = [coefsi;mmi];
      fprintf(['with ' kind ' coefficients of degree 0 :\n'])
      disp(coefsi(:,1))
      for j=1:ki
         fprintf(['and ' kind ' coefficients of degree ',num2str(j),' :\n'])
         disp(coefsi(:,blaisei(end,j+1)+1:blaisei(end,j+2)))
      end
   else
      out1=coefsi; d=di; k=ki; center=centeri; blaise=blaisei; mm=mmi;
   end
end
