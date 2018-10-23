% integrating different functions
% function int_testing_NM

clc
clear
global kappa
%% weighting gaussian kernel  properties
n=2;
mu=ones(n,1);
% A=randn(n,n).^2;
%  P=10*A*A';
%  P=[20.5058   69.2898
%    69.2898  234.3919];

 P=eye(n);
%% function to be integrated
% F=@(x)(sqrt(1+x'*x)^3);
 F=@(x)(sqrt(1+sum(x.^2,2)).^(-3));
% F=@(x,n)(normpdf(x,ones(n,1),eye(n)));
% F=@(x)((+20*x(:,2).^2.*x(:,1).^2+x(:,1)-x(:,1).^2).^5);
%% montecarlo integration
j=1;
for i=10^4:10^4:10^5
[x,w]=monte_carlo_int_normal(mu,P,i);
m_mc(j)=sum(w.*F(x));
j=j+1;
end
figure(1)
plot(10^4:10^4:10^5,m_mc)
set(gca,'FontSize',16)
title('MC')
xlabel('No. of points')
ylabel('integral value')
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 16)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
% [m_mc,P_mc]=prop_mean_cov_points_discr(x,w,F);
%% GH integration
for i=4:1:21
[x,w]=GH_points(mu,P,i);
m_gh(i-3)=sum(w.*F(x));
end
figure(2)
plot(4:1:21,m_gh);
set(gca,'FontSize',16)
title('GH')
xlabel('No. of points in 1D')
ylabel('integral value')
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 16)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)



% [m_gh,P_gh]=prop_mean_cov_points_discr(x,w,F);

%% CKF integration
[x,w]=cubature_KF_points(mu,P);
m_ckf=sum(w.*F(x));
% [m_ckf,P_ckf]=prop_mean_cov_points_discr(x,w,F);
% 
%% Cheng ckf integration
[x,w]=cubature_KF_points(mu,P);
m_ckf=sum(w.*F(x));
% [m_ckf,P_ckf]=prop_mean_cov_points_discr(x,w,F);
% 
%% UT integration
kappa=1;
[x,w]=UT_sigmapoints(mu,P,2);
m_ut=sum(w.*F(x));
% [m_ut,P_ut]=prop_mean_cov_points_discr(x,w,F);

%% NM4 integration
[x,w]=conjugate_dir_gausspts(mu,P);
m_nm4=sum(w.*F(x));
% [m_nm4,P_nm4]=prop_mean_cov_points_discr(x,w,F);

%% NM6 integration
[x,w]=conjugate_dir_gausspts_till_6moment_scheme2(mu,P);
m_nm6=sum(w.*F(x));

%% NM8 integration
[x,w]=conjugate_dir_gausspts_till_8moment(mu,P);
m_nm6=sum(w.*F(x));


m_mc=m_mc(end);
m_gh=m_gh(end);
[m_mc,m_gh,m_ckf,m_ut,m_nm4]
100*[abs(m_mc-m_gh),abs(m_ckf-m_gh),abs(m_ut-m_gh),abs(m_nm4-m_gh)]./abs(m_gh)


% [x1,x2]=meshgrid(-10:0.1:10);
% for i=1:1:length(x1)
%     for j=1:1:length(x2)
%         f(i,j)=F([x1(i,j),x2(i,j)]);
%     end
% end
% mesh(x1,x2,f)