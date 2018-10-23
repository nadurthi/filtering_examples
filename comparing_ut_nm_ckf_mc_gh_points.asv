% comparing the UT,GH,MC,CKF,NM
clc
clear
%% function to be integrated
F=@(x)(sqrt(1+sum(x.*x,2)).^(6));
n=9;% Dimension of system
mu=zeros(n,1);
% P=diag(randi([1,5],n,1).^2)
P=100*eye(n);
%% firstly monte carlo
% i=1;
% Int_mc=0;
% for N=10000:1000:50000
% [x,w]=monte_carlo_int_normal(mu,P,N);
% Int_mc(i)=sum(w.*F(x));
% i=i+1;
% end
% figure(1)
% plot(10000:1000:50000,Int_mc)
% Int_mc_mean=mean(Int_mc(floor(3*end/4):1:end))
%% second is CKF points
% [xcb,wcb]=cubature_KF_points(mu,P);
% Int_cb=sum(wcb.*F(xcb))

%% third GH points
i=1;
Int_gh=0;
for N=3:1:5
[xgh,wgh]=GenerateQuadPoints(P,mu,N);
Int_gh(i)=sum(wgh.*F(xgh));
i=i+1;
end
% figure(3)
% plot(3:1:15,Int_gh)
% Int_gh_mean=mean(Int_gh(floor(3*end/4):1:end))
%% fourth UT sigma points
% [x2,w2]=UT_sigmapoints(mu,P,2);
% Int_ut2=sum(w2.*F(x2))
% 
% [x4,w4]=UT_sigmapoints(mu,P,4);
% Int_ut4=sum(w4.*F(x4))
% 
% [x6,w6]=UT_sigmapoints(mu,P,6);
% Int_ut6=sum(w6.*F(x6))
%% fifth NM sigma points
% [xnm,wnm]=conjugate_dir_gausspts_till_8moment(mu,P);
[xnm,wnm]=conjugate_dir_gausspts_till_6moment_scheme2(mu,P);
Int_nm=sum(wnm.*F(xnm));

100*abs([Int_gh(1),Int_gh(2),Int_nm]-Int_gh(3))./Int_gh(3)

% [100*abs(Int_nm-Int_gh_mean)/Int_gh_mean,100*abs(Int_cb-Int_gh_mean)/Int_gh_mean,100*abs(Int_ut2-Int_gh_mean)/Int_gh_mean,100*abs(Int_ut4-Int_gh_mean)/Int_gh_mean,100*abs(Int_ut6-Int_gh_mean)/Int_gh_mean]
