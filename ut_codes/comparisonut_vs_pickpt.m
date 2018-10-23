clc
clear
% P=[4,2,1;2,9,1;1,1,16];
 P=10*eye(7);
%  P=[4,1,2,1;1,9,2,3;2,2,16,4;1,3,4,25];
mu=zeros(7,1);
% F=@(x1,x2,x3)(x1^4+x2^4+x3^4+x1^3*x2+x1^2*x2^2+x3^2*x2^2+x1^2*x3^2+x1^3*x3+x2^3*x3++x3^3*x2);
% F=@(x1,x2,x3,x4,x5,x6,x7,x8)(x1^4+x8^2*x2^2+x3^2+x4^2*x5^2);
% F=@(x1,x2,x3,x4)(x1^6+x2^6+x3^6+x4^6+x1^4+x2^4+x3^4+x4^4+x1^2*x2^2*x3^2+...
% +x1^3*x3+x1^2*x2^4+x3^4*x2^2+x1^2*x3^4+x1^3*x3+x2^3*x3+x3^3*x2);
F=@(x1,x2,x3,x4,x5,x6,x7)(x1^6+x1^2*x2^2*x3^2+x1^4*x6^2+x5^4+x3^2*x4^2);
%***********************************************************
[xgh,wgh]=GenerateQuadPoints(P,mu,5);
sumgh=0;
for i=1:1:length(xgh)
    sumgh=sumgh+wgh(i)*F(xgh(i,1),xgh(i,2),xgh(i,3),xgh(i,4),xgh(i,5),xgh(i,6),xgh(i,7));
end

%***********************************************************
% [xgh,wgh10]=GenerateQuadPoints(P,mu,9);
% sumgh10=0;
% for i=1:1:length(xgh)
%     sumgh10=sumgh10+wgh10(i)*F(xgh(i,1),xgh(i,2),xgh(i,3),xgh(i,4),xgh(i,5),xgh(i,6));
% end

%***********************************************************
[X,wpk]=conjugate_dir_gausspts_till_6moment_scheme2(mu,P);
sumutpk=0;

for i=1:1:length(X)
    sumutpk=sumutpk+wpk(i)*F(X(i,1),X(i,2),X(i,3),X(i,4),X(i,5),X(i,6),X(i,7));
end

% %***************************************************
[sumgh,sumutpk]
% [100*abs(sumgh-sumgh10)/sumgh10,100*abs(sumutpk-sumgh10)/sumgh10]
% Ngh10=length(wgh10)
% Ngh=length(wgh)
% Nnm=length(wpk)
% [mgh';mpk']%;mut';m6n1']


