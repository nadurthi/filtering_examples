% integrating different functions
% function int_testing_NM

% clc
clear
global kappa
%% function to be integrated

  F=@(x,mu,P)(sqrt(1+sum(x.^2,2)).^(-3));

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
N_spar=0;
    for n=Ds:1:Dl
        n
    s=s+1;
mu=zeros(n,1);
P=0.1*eye(n);


%% GH integration
j=1;
Ns=2;
Nl=13;
for i=Ns:1:Nl
[x,w]=GH_points(mu,P,i);
m_gh(s,j)=sum(w.*F(x,mu,P));
j=j+1;
end
m_gh_n(s)=m_gh(s,end);


%% GH-Sparse Grid integration
clear x w
j=1;
Ns=2;
Nl=6;

for i=Ns:1:Nl
%    [x,w] = nwspgr('GQN',n, i);
   [x,w]=smolyak_sparse_grid(n,i,'GH');

   for pp=1:1:length(w)
   x(pp,:)=mu+sqrtm(P)*x(pp,:)';
   end
m_gh_spar(s,j)=sum(w.*F(x,mu,P));
N_spar(s,j)=length(w);
j=j+1;
end
% m_gh_spar_n(s)=m_gh_spar(s,end);
%% CKF integration
clear x w
[x,w]=cubature_KF_points(mu,P);
m_ckf=sum(w.*F(x,mu,P));

m_ckf_n(s)=m_ckf;
N_ckf(s)=length(w);
% [m_ckf,P_ckf]=prop_mean_cov_points_discr(x,w,F);
% 

%% UT integration k1
kappa=1;
clear x w
[x,w]=UT_sigmapoints(mu,P,2);
m_ut=sum(w.*F(x,mu,P));

m_ut_n_k1(s)=m_ut;
N_ut_k1(s)=length(w);


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
% 
% PP=size(m_gh(:,2:6));
% m_gh_n_mat=[];
% for i=1:1:PP(2)
% m_gh_n_mat=[m_gh_n_mat,m_gh_n'];
% end
% A=[(Ds:1:Dl)',m_gh(:,2:6),m_gh_n',m_ckf_n',m_ut_n',m_ut_n_k1',m_ut_n_k2',m_nm4_n',m_nm6_n',m_nm8_n'];
% AA={'Dim','GH3','GH4','GH5','GH6','GH7','GH9','CKF','UT_K0','UT_K1','UT_K2','CUT4','CUT6','CUT8'};
% B=[(Ds:1:Dl)',100*abs(m_gh(:,2:6)-m_gh_n_mat)./m_gh_n_mat,100*abs(m_ckf_n'-m_gh_n')./m_gh_n',100*abs(m_ut_n'-m_gh_n')./m_gh_n',100*abs(m_ut_n_k1'-m_gh_n')./m_gh_n',100*abs(m_ut_n_k2'-m_gh_n')./m_gh_n',100*abs(m_nm4_n'-m_gh_n')./m_gh_n',100*abs(m_gh_n'-m_nm6_n')./m_gh_n',100*abs(m_nm8_n'-m_gh_n')./m_gh_n'];
%  C=[(Ds:1:Dl)',3.^((Ds:1:Dl)'),4.^((Ds:1:Dl)'),5.^((Ds:1:Dl)'),6.^((Ds:1:Dl)'),7.^((Ds:1:Dl)'),9.^((Ds:1:Dl)'),N_ckf',N_ut',N_nm4',N_nm6',N_nm8'];
% BB={'Dim','GH3','GH4','GH5','GH6','GH7','CKF','UT_K0','UT_K1','UT_K2','CUT4','CUT6','CUT8'};
% CC={'Dim','GH3','GH4','GH5','GH6','GH7','GH9','CKF','UT','CUT4','CUT6','CUT8'};
% size(A)
% size(B)
% size(C)
% xlswrite('ck_eg_p_m3.xls',AA,num2str(P(1,1)),'A1:N1')
% xlswrite('ck_eg_p_m3.xls',BB,num2str(P(1,1)),'A8:M8')
% xlswrite('ck_eg_p_m3.xls',CC,num2str(P(1,1)),'A15:L15')
% xlswrite('ck_eg_p_m3.xls',A,num2str(P(1,1)),'A2:N6')
% xlswrite('ck_eg_p_m3.xls',B,num2str(P(1,1)),'A9:M13')
% xlswrite('ck_eg_p_m3.xls',C,num2str(P(1,1)),'A16:L20')
% xlswrite('ck_eg_p_m3.xls','P',num2str(P(1,1)),'A22:A22')
% xlswrite('ck_eg_p_m3.xls',P,num2str(P(1,1)),'A23:F28')
% 

%% calculating the error
truth=m_gh_n'; % for all dimensions
% GH2:10
Egh=100*abs((m_gh(:,1:9)-repmat(truth,1,9))./repmat(truth,1,9));
Ngh=[(2:10).^2;(2:10).^3;(2:10).^4;(2:10).^5;(2:10).^6];
% Sparse GH2:10

Egh_spar=100*abs((m_gh_spar(:,1:end)-repmat(truth,1,5))./repmat(truth,1,5));
% CUT4
Ecut=[100*abs((m_nm4_n'-truth)./truth),100*abs((m_nm6_n'-truth)./truth),100*abs((m_nm8_n'-truth)./truth)];
Ncut=[N_nm4',N_nm6',N_nm8'];
% CUT6
% Ecut6=100*abs((m_nm6_n'-truth)./truth);
% % CUT8
% Ecut8=100*abs((m_nm8_n'-truth)./truth);
%% %%%%%%%%%%%%
% figure(7)
% bar3(2:1:6,C(:,[3,12]))
% set(gca,'XTickLabel',{'GH4','CUT8'})
% xlabel('Dimension')
% ylabel('Method')
% % title('Rel. error vs No. of points from CUT4,CUT6,CUT8 each in 2D,3D,4d,5D and 6D')
% set(gca,'FontSize',16)
% h = get(gca, 'title');
% k = get(gca, 'xlabel');
% l = get(gca, 'ylabel');
% set(h, 'FontName', 'Helvetica', 'FontSize', 16)
% set(k, 'FontName', 'Helvetica', 'FontSize', 16)
% set(l, 'FontName', 'Helvetica', 'FontSize', 16)
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


figure
    plot(Ngh(1,2:6),Egh(1,2:6),'ko--','linewidth',2,'MarkerSize',15)
        set(gca,'XScale','log')

legend('2D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
axis([1,10^6,0,2])
plot_prop_paper
mov(1) = getframe(gcf);
pause(2)
figure
hold on
plot(Ngh(1,2:6),Egh(1,2:6),'ko--','linewidth',2,'MarkerSize',15)
    plot(Ngh(2,2:6),Egh(2,2:6),'kd--','linewidth',2,'MarkerSize',15)
    set(gca,'XScale','log')
axis([1,10^6,0,2])
legend('2D','3D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper
mov(2) = getframe(gcf);
    figure
    hold on
    plot(Ngh(1,2:6),Egh(1,2:6),'ko--','linewidth',2,'MarkerSize',15)
    plot(Ngh(2,2:6),Egh(2,2:6),'kd--','linewidth',2,'MarkerSize',15)
    plot(Ngh(3,2:6),Egh(3,2:6),'ks--','linewidth',2,'MarkerSize',15)
    set(gca,'XScale','log')
axis([1,10^6,0,2])
legend('2D','3D','4D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper
mov(3) = getframe(gcf);
    figure
    hold on
        plot(Ngh(1,2:6),Egh(1,2:6),'ko--','linewidth',2,'MarkerSize',15)
    plot(Ngh(2,2:6),Egh(2,2:6),'kd--','linewidth',2,'MarkerSize',15)
    plot(Ngh(3,2:6),Egh(3,2:6),'ks--','linewidth',2,'MarkerSize',15)
    plot(Ngh(4,2:6),Egh(4,2:6),'kh--','linewidth',2,'MarkerSize',15)
    set(gca,'XScale','log')
axis([1,10^6,0,2])
legend('2D','3D','4D','5D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper
mov(4) = getframe(gcf);

    figure
    hold on
            plot(Ngh(1,2:6),Egh(1,2:6),'ko--','linewidth',2,'MarkerSize',15)
    plot(Ngh(2,2:6),Egh(2,2:6),'kd--','linewidth',2,'MarkerSize',15)
    plot(Ngh(3,2:6),Egh(3,2:6),'ks--','linewidth',2,'MarkerSize',15)
    plot(Ngh(4,2:6),Egh(4,2:6),'kh--','linewidth',2,'MarkerSize',15)
    plot(Ngh(5,2:6),Egh(5,2:6),'k^--','linewidth',2,'MarkerSize',15)

axis([1,10^6,0,2])
    set(gca,'XScale','log')

legend('2D','3D','4D','5D','6D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper
mov(5) = getframe(gcf);

movie2avi(mov, 'GHerrm3.avi', 'compression', 'None');











figure
hold on
    plot(Ncut(1,:),Ecut(1,:),'ko--','linewidth',2,'MarkerSize',15)
 set(gca,'XScale','log')
axis([1,10^3,0,2.5])
legend('2D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper
mov(1) = getframe(gcf);   
  figure
hold on
 plot(Ncut(1,:),Ecut(1,:),'ko--','linewidth',2,'MarkerSize',15)
    plot(Ncut(2,:),Ecut(2,:),'kd--','linewidth',2,'MarkerSize',15)
set(gca,'XScale','log')
axis([1,10^3,0,2.5])
legend('2D','3D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper
mov(2) = getframe(gcf);    
   figure
hold on 
 plot(Ncut(1,:),Ecut(1,:),'ko--','linewidth',2,'MarkerSize',15)
    plot(Ncut(2,:),Ecut(2,:),'kd--','linewidth',2,'MarkerSize',15)
    plot(Ncut(3,:),Ecut(3,:),'ks--','linewidth',2,'MarkerSize',15)
 set(gca,'XScale','log')
axis([1,10^3,0,2.5])
legend('2D','3D','4D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper
mov(3) = getframe(gcf);   
   figure
hold on 
 plot(Ncut(1,:),Ecut(1,:),'ko--','linewidth',2,'MarkerSize',15)
    plot(Ncut(2,:),Ecut(2,:),'kd--','linewidth',2,'MarkerSize',15)
    plot(Ncut(3,:),Ecut(3,:),'ks--','linewidth',2,'MarkerSize',15)
    plot(Ncut(4,:),Ecut(4,:),'kh--','linewidth',2,'MarkerSize',15)
 set(gca,'XScale','log')
axis([1,10^3,0,2.5])
legend('2D','3D','4D','5D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper
mov(4) = getframe(gcf);   
    figure
hold on
 plot(Ncut(1,:),Ecut(1,:),'ko--','linewidth',2,'MarkerSize',15)
    plot(Ncut(2,:),Ecut(2,:),'kd--','linewidth',2,'MarkerSize',15)
    plot(Ncut(3,:),Ecut(3,:),'ks--','linewidth',2,'MarkerSize',15)
    plot(Ncut(4,:),Ecut(4,:),'kh--','linewidth',2,'MarkerSize',15)
    plot(Ncut(5,:),Ecut(5,:),'k^--','linewidth',2,'MarkerSize',15)
set(gca,'XScale','log')
axis([1,10^3,0,2.5])
legend('2D','3D','4D','5D','6D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper
mov(5) = getframe(gcf);


movie2avi(mov, 'CUTerrm3.avi', 'compression', 'None');

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on
% for ct = 1:2
    plot(N_spar(1,2:5),Egh_spar(1,2:5),'ko--','linewidth',2,'MarkerSize',15)
    plot(N_spar(2,2:5),Egh_spar(2,2:5),'kd--','linewidth',2,'MarkerSize',15)
    plot(N_spar(3,2:5),Egh_spar(3,2:5),'ks--','linewidth',2,'MarkerSize',15)
    plot(N_spar(4,2:5),Egh_spar(4,2:5),'kh--','linewidth',2,'MarkerSize',15)
    plot(N_spar(5,2:5),Egh_spar(5,2:5),'k^--','linewidth',2,'MarkerSize',15)


%     set(gca,'YTick',0:0.5:10)
    set(gca,'XScale','log')
   % end
legend('2D','3D','4D','5D','6D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper

figure(4)
bar3(2:1:6,[Ngh(:,3),Ncut(:,3)])
set(gca,'XTickLabel',{'GH4','CUT8'})
xlabel('Dimension')
ylabel('Method')
zlabel('Number of points')

figure(5)
bar3(2:1:6,[Ngh(:,2),N_spar(:,2),N_nm4'])
set(gca,'XTickLabel',{'GH3','SPARSE-GH3','CUT4'})
xlabel('Dimension')
ylabel('Method')
zlabel('Number of points')
plot_prop_paper

figure(6)
bar3(2:1:6,[Ngh(:,5),N_spar(:,5)])
set(gca,'XTickLabel',{'GH6','SPARSE-GH6'})
xlabel('Dimension')
ylabel('Method')
zlabel('Number of points')
plot_prop_paper

%absolute value of the sum of weights
% A=0;
% for n=3:1:9
%     for i=2:1:5
%         [n,i]
%      [x,w]=smolyak_sparse_grid(n,i,'GH');
%      A(n-2,i-1)=sum(abs(w));
%     end
% end
% figure
% hold on
% % for ct = 1:2
%     plot(2:1:5,A(1,:),'ko--','linewidth',2,'MarkerSize',15)
%     plot(2:1:5,A(2,:),'kd--','linewidth',2,'MarkerSize',15)
%     plot(2:1:5,A(3,:),'ks--','linewidth',2,'MarkerSize',15)
%     plot(2:1:5,A(4,:),'kh--','linewidth',2,'MarkerSize',15)
%     plot(2:1:5,A(5,:),'k^--','linewidth',2,'MarkerSize',15)
%     plot(2:1:5,A(6,:),'k*--','linewidth',2,'MarkerSize',15)
%     plot(2:1:5,A(7,:),'kp--','linewidth',2,'MarkerSize',15)
% 
% legend('3D','4D','5D','6D','7D','8D','9D')
% xlabel('level k')
% ylabel('Sum of absolute value of weights')
% plot_prop_paper