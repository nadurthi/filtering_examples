function [y,lam,xl,xu]=MaxEntPdf(y,M,xl,xu,lam0,method)
% xl and xu are the lower and uppre bounds for the 
% y contains all the combinations for the moments
% M starts from the second moment

ns=size(y,2);
nm=size(y,1);

% %quad points or samples
% II=2*eye(ns);
% P=zeros(1,ns);
% for i=1:1:ns
%      ind=find(sum(abs(y-repmat(II(i,:),nm,1)),2)==0);
%      P(i)=M(ind);
% end
% 
% xl=-bndThres*sqrt(P);
% xu=bndThres*sqrt(P);
% 
if strcmp(method,'GL')
[X,W] = GLeg_pts(10*ones(size(xl)), xl, xu);
W=prod(abs(xu-xl))*W;
elseif strcmp(method,'GH')
    pp=0.5;
[X,W] = GH_points(zeros(size(y,2),1),pp*eye(size(y,2)),10);
W=W./mvnpdf(X, zeros(size(y,2),1)', pp*eye(size(y,2)));
elseif strcmp(method,'GLag')
a=0;
b=4;
alpha=0;
[X,W]=gen_laguerre_rule ( 20, alpha, a, b, 'test' );
W=W./((abs(X-a).^alpha).*exp(-b*(X-a)));
end

options=optimset('disp','iter','TolCon',1e-8,'MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-8,'TolX',1e-8,'Jacobian','on');
% lam=fmincon(@(lam)1,zeros(nm,1),[],[],[],[],[],[],@(lam)maxentfcn(lam,y,M,X,W),options);
lam=fsolve(@(lam)maxentFSOLVE(lam,y,M,X,W),lam0,options);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [c,ceq]=maxentfcn(lam,y,M,X,W)
% 
% nm=size(y,1);
% ceq=0;
% for i=1:1:length(W)
%     ceq=ceq+W(i)*prod(repmat(X(i,:),nm,1).^y,2)*pdf_MaxEnt(X(i,:),lam,y);
%     
% end
% ceq=ceq-M;
% c=[];
% 
% end

