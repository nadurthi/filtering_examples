% integrating different functions
% function int_testing_NM

clc
clear
global kappa
%% function to be integrated
% F=@(x)(sqrt(1+x'*x)^3);
  F=@(x,mu,P)(sqrt(1+sum(x.^2,2)).^(-3));
% F=@(x,mu,P)(1./(1+exp(-sum(x,2))));
% F=@(x,mu,P)(normpdf(x,mu,P));
% F=@(x,mu,P)((1/(2*pi)^(n/2))*(1/sqrt(det(P)))*exp(-0.5*sum((x).^2,2)));
% F=@(x,mu,P)(exp(-sum(x,2)));
% F=@(x)((+20*x(:,2).^2.*x(:,1).^2+x(:,1)-x(:,1).^2).^5);
%% weighting gaussian kernel  properties
m_gh_n=0;
m_ckf_n=0;
m_ut_n=0;
m_ut_n_k1=0;
m_ut_n_k2=0;
m_nm4_n=0;
m_nm6_n=0;
m_nm8_n=0;
Ds=2;
Dl=6;
s=0;
m_gh=0;
for n=Ds:1:Dl
    s=s+1;
mu=zeros(n,1);
% A=10*rand(n,n).^2;
% P=A*A';
P=0.1*eye(n);

%% montecarlo integration
% j=1;
% for i=10^4:10^4:10^5
% [x,w]=monte_carlo_int_normal(mu,P,i);
% m_mc(j)=sum(w.*F(x));
% j=j+1;
% end
% figure(1)
% plot(10^4:10^4:10^5,m_mc)
% set(gca,'FontSize',16)
% title('MC')
% xlabel('No. of points')
% ylabel('integral value')
% h = get(gca, 'title');
% k = get(gca, 'xlabel');
% l = get(gca, 'ylabel');
% set(h, 'FontName', 'Helvetica', 'FontSize', 16)
% set(k, 'FontName', 'Helvetica', 'FontSize', 16)
% set(l, 'FontName', 'Helvetica', 'FontSize', 16)
% % [m_mc,P_mc]=prop_mean_cov_points_discr(x,w,F);
%% GH integration
j=1;
Ns=2;
Nl=15;
for i=Ns:1:Nl
[x,w]=GH_points(mu,P,i);
m_gh(s,j)=sum(w.*F(x,mu,P));
j=j+1;
end
m_gh_n(s)=m_gh(s,end);
figure(1)
hold on
plot(Ns:1:Nl,m_gh,'LineWidth',2);

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
clear x w
[x,w]=cubature_KF_points(mu,P);
m_ckf=sum(w.*F(x,mu,P));

m_ckf_n(s)=m_ckf;
N_ckf(s)=length(w);
% [m_ckf,P_ckf]=prop_mean_cov_points_discr(x,w,F);
% 
%% UT integration k0
kappa=0;
clear x w
[x,w]=UT_sigmapoints(mu,P,2);
m_ut=sum(w.*F(x,mu,P));

m_ut_n(s)=m_ut;
N_ut(s)=length(w);

%% UT integration k1
kappa=1;
clear x w
[x,w]=UT_sigmapoints(mu,P,2);
m_ut=sum(w.*F(x,mu,P));

m_ut_n_k1(s)=m_ut;
N_ut_k1(s)=length(w);

%% UT integration k2
kappa=2;
clear x w
[x,w]=UT_sigmapoints(mu,P,2);
m_ut=sum(w.*F(x,mu,P));

m_ut_n_k2(s)=m_ut;
N_ut_k2(s)=length(w);
%% NM4 integration
clear x w
[x,w]=conjugate_dir_gausspts(mu,P);
m_nm4=sum(w.*F(x,mu,P));

m_nm4_n(s)=m_nm4;
N_nm4(s)=length(w);
% [m_nm4,P_nm4]=prop_mean_cov_points_discr(x,w,F);

%% NM6 integration
clear x w
[x,w]=conjugate_dir_gausspts_till_6moment_scheme2(mu,P);
m_nm6=sum(w.*F(x,mu,P));

m_nm6_n(s)=m_nm6;
N_nm6(s)=length(w);

%% NM8 integration
clear x w
[x,w]=conjugate_dir_gausspts_till_8moment(mu,P);
m_nm8=sum(w.*F(x,mu,P));

m_nm8_n(s)=m_nm8;
N_nm8(s)=length(w);

end
%%
PP=size(m_gh(:,2:6));
m_gh_n_mat=[];
for i=1:1:PP(2)
m_gh_n_mat=[m_gh_n_mat,m_gh_n'];
end
A=[(Ds:1:Dl)',m_gh(:,2:6),m_gh_n',m_ckf_n',m_ut_n',m_ut_n_k1',m_ut_n_k2',m_nm4_n',m_nm6_n',m_nm8_n'];
AA={'Dim','GH3','GH4','GH5','GH6','GH7','GH9','CKF','UT_K0','UT_K1','UT_K2','CUT4','CUT6','CUT8'};
B=[(Ds:1:Dl)',100*abs(m_gh(:,2:6)-m_gh_n_mat)./m_gh_n_mat,100*abs(m_ckf_n'-m_gh_n')./m_gh_n',100*abs(m_ut_n'-m_gh_n')./m_gh_n',100*abs(m_ut_n_k1'-m_gh_n')./m_gh_n',100*abs(m_ut_n_k2'-m_gh_n')./m_gh_n',100*abs(m_nm4_n'-m_gh_n')./m_gh_n',100*abs(m_gh_n'-m_nm6_n')./m_gh_n',100*abs(m_nm8_n'-m_gh_n')./m_gh_n'];
C=[(Ds:1:Dl)',3.^((Ds:1:Dl)'),4.^((Ds:1:Dl)'),5.^((Ds:1:Dl)'),6.^((Ds:1:Dl)'),7.^((Ds:1:Dl)'),9.^((Ds:1:Dl)'),N_ckf',N_ut',N_nm4',N_nm6',N_nm8'];
BB={'Dim','GH3','GH4','GH5','GH6','GH7','CKF','UT_K0','UT_K1','UT_K2','CUT4','CUT6','CUT8'};
CC={'Dim','GH3','GH4','GH5','GH6','GH7','GH9','CKF','UT','CUT4','CUT6','CUT8'};
size(A)
size(B)
size(C)
xlswrite('ck_eg_p_m3.xls',AA,num2str(P(1,1)),'A1:N1')
xlswrite('ck_eg_p_m3.xls',BB,num2str(P(1,1)),'A8:M8')
xlswrite('ck_eg_p_m3.xls',CC,num2str(P(1,1)),'A15:L15')
xlswrite('ck_eg_p_m3.xls',A,num2str(P(1,1)),'A2:N6')
xlswrite('ck_eg_p_m3.xls',B,num2str(P(1,1)),'A9:M13')
xlswrite('ck_eg_p_m3.xls',C,num2str(P(1,1)),'A16:L20')
xlswrite('ck_eg_p_m3.xls','P',num2str(P(1,1)),'A22:A22')
xlswrite('ck_eg_p_m3.xls',P,num2str(P(1,1)),'A23:F28')


% figure(1)
% legend('2D','3D','4D','5D','6D')
% figure(2)
% plot(2:1:6,B(:,2),'+',2:1:6,B(:,3),'o',2:1:6,B(:,4),'*',2:1:6,B(:,5),'s',2:1:6,B(:,6),'d',2:1:6,B(:,7),'.',2:1:6,B(:,8),'x',2:1:6,B(:,9),'p',2:1:6,B(:,10),'h',2:1:6,B(:,11),'^',2:1:6,B(:,12),'<',2:1:6,B(:,13),'>')
% legend('GH3','GH4','GH5','GH6','GH7','CKF','UT_K0','UT_K1','UT_K2','CUT4','CUT6','CUT8')
% figure(3)
% bar3(2:1:6,B(:,2:13))
% set(gca,'XTickLabel',{'GH3','GH4','GH5','GH6','GH7','CKF','UT0','UT1','UT2','CUT4','CUT6','CUT8'})
% colormap(gray)
% figure(4)
% bar3(2:1:6,C(:,[2,3,4,10,11,12]))
% set(gca,'XTickLabel',{'GH3','GH4','GH5','CUT4','CUT6','CUT8'})
% colormap(gray)
% figure(5)
% bar3(2:1:6,B(:,2:6))
% set(gca,'XTickLabel',{'GH3','GH4','GH5','GH6','GH7','CKF','UT0','UT1','UT2','CUT4','CUT6','CUT8'})
% colormap(gray)
% figure(6)
% plot(2:1:6,B(:,2),'^--',2:1:6,B(:,3),'o--',2:1:6,B(:,4),'*--',2:1:6,B(:,5),'s--',2:1:6,B(:,6),'d--')
% legend('GH3','GH4','GH5','GH6','GH7')
figure(7)
bar3(2:1:6,C(:,[3,12]))
set(gca,'XTickLabel',{'GH4','CUT8'})
xlabel('Dimension')
ylabel('Method')
% title('Rel. error vs No. of points from CUT4,CUT6,CUT8 each in 2D,3D,4d,5D and 6D')
set(gca,'FontSize',16)
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 16)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
% colormap(gray)
% figure(8)
% % hold on
% % plot3(2:1:6,C(:,2),B(:,2),'bo')
% % plot3(2:1:6,C(:,3),B(:,3),'bd')
% % plot3(2:1:6,C(:,4),B(:,4),'bs')
% % legend('GH3','GH4','GH5')
% mesh(2:1:6,[3 4 5],B(:,2:4)')
% % set(gca,'XTickLabel',{'GH3'})%,'GH4','GH5','GH6','GH7','CKF','UT0','UT1','UT2','CUT4','CUT6','CUT8'})
% % colormap(gray)
% clear f x1 x2
% 
% F=@(x,mu,P)(sqrt(1+sum(x.^2,2)).^(-3));
% [x1,x2]=meshgrid(-5:0.1:5);
% for i=1:1:length(x1)
%     for j=1:1:length(x2)
%         f(i,j)=F([x1(i,j),x2(i,j)],mu,P);
%     end
% end
% figure(9)
% mesh(x1,x2,f)
% 
% for ndim=2:6
% for gh = 3:5
% np(ndim-1,gh-2)=gh^ndim;
% dime(ndim-1,gh-2) = ndim;
% end
% end
%  figure(10)
% surf(dime,np,B(:,2:4))


figure(11)
hold on
% for ct = 1:2
    plot((C(1,2:6)),B(1,2:6),'ko--','linewidth',2,'MarkerSize',15)
    plot((C(2,2:6)),B(2,2:6),'kd--','linewidth',2,'MarkerSize',15)
    plot((C(3,2:6)),B(3,2:6),'ks--','linewidth',2,'MarkerSize',15)
    plot((C(4,2:6)),B(4,2:6),'kh--','linewidth',2,'MarkerSize',15)
    plot((C(5,2:6)),B(5,2:6),'k^--','linewidth',2,'MarkerSize',15)

    plot((C(1,2)),B(1,2),'ko--','linewidth',2,'MarkerSize',15)
    plot((C(2,2)),B(2,2),'kd--','linewidth',2,'MarkerSize',15)
    plot((C(3,2)),B(3,2),'ks--','linewidth',2,'MarkerSize',15)
    plot((C(4,2)),B(4,2),'kh--','linewidth',2,'MarkerSize',15)
    plot((C(5,2)),B(5,2),'k^--','linewidth',2,'MarkerSize',15)
    set(gca,'YTick',0:0.25:3)
% end
legend('2D','3D','4D','5D','6D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')

% title('Rel error vs total no. points from GH3,GH4,GH5,GH6,GH7 each in 2D,3D,4d,5D and 6D')
set(gca,'FontSize',20)
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 20)
set(k, 'FontName', 'Helvetica', 'FontSize', 20)
set(l, 'FontName', 'Helvetica', 'FontSize', 20)

figure(12)
hold on
% for ct = 1:2
    plot((C(1,10:12)),B(1,11:13),'ko--','linewidth',2,'MarkerSize',15)
    plot((C(2,10:12)),B(2,11:13),'kd--','linewidth',2,'MarkerSize',15)
    plot((C(3,10:12)),B(3,11:13),'ks--','linewidth',2,'MarkerSize',15)
    plot((C(4,10:12)),B(4,11:13),'kh--','linewidth',2,'MarkerSize',15)
    plot((C(5,10:12)),B(5,11:13),'k^--','linewidth',2,'MarkerSize',15)

    plot((C(1,10)),B(1,11),'ko--','linewidth',2,'MarkerSize',15)
    plot((C(2,10)),B(2,11),'kd--','linewidth',2,'MarkerSize',15)
    plot((C(3,10)),B(3,11),'ks--','linewidth',2,'MarkerSize',15)
    plot((C(4,10)),B(4,11),'kh--','linewidth',2,'MarkerSize',15)
    plot((C(5,10)),B(5,11),'k^--','linewidth',2,'MarkerSize',15)
    set(gca,'YTick',0:0.25:3)
   % end
legend('2D','3D','4D','5D','6D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')

% title('Rel. error vs No. of points from CUT4,CUT6,CUT8 each in 2D,3D,4d,5D and 6D')
set(gca,'FontSize',20)
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 20)
set(k, 'FontName', 'Helvetica', 'FontSize', 20)
set(l, 'FontName', 'Helvetica', 'FontSize', 20)