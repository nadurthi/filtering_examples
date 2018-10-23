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
Dl=8;
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
Nl=6;
for i=Ns:1:Nl
% [x,w]=GH_points(mu,P,i);
[x,w] = GLeg_pts(i*ones(1,n),-1*ones(1,n), ones(1,n));
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
   [x,w]=Sigma__Quadrature_points(-1*ones(1,n),ones(1,n),i,'SGLgn','uniform');
%    [x,w]=smolyak_sparse_grid(n,i,'GLgn');

%    for pp=1:1:length(w)
%    x(pp,:)=mu+sqrt(P)*x(pp,:)';
%    end
m_gh_spar(s,j)=sum(w.*F(x,mu,P));
N_spar(s,j)=length(w);
j=j+1;
end
% m_gh_spar_n(s)=m_gh_spar(s,end);

%% CUT-U integration
clear x w
if n<=8
[x4,w4]=uniform_sigma_pts(-1*ones(1,n), ones(1,n),4);
m_nm4=sum(w4.*F(x4,mu,P));
m_nm4_n(s)=m_nm4;
N_nm4(s)=length(w4);
else
    m_nm4=NaN;
m_nm4_n(s)=NaN;
N_nm4(s)=NaN;
end

[x6,w6]=uniform_sigma_pts(-1*ones(1,n), ones(1,n),6);
m_nm6=sum(w6.*F(x6,mu,P));
m_nm6_n(s)=m_nm6;
N_nm6(s)=length(w6);

if n<=5
[x8,w8]=uniform_sigma_pts(-1*ones(1,n), ones(1,n),8);
m_nm8=sum(w8.*F(x8,mu,P));
m_nm8_n(s)=m_nm8;
N_nm8(s)=length(w8);
else
    m_nm8=NaN;
m_nm8_n(s)=NaN;
N_nm8(s)=NaN;
end



end


%% calculating the error
[m_gh_n;m_gh_spar(:,end)';m_nm4_n;m_nm6_n;m_nm8_n]
truth=m_gh_n'; % for all dimensions
% GH2:10
Egh=100*abs((m_gh(:,1:4)-repmat(truth,1,4))./repmat(truth,1,4));
Ngh=[(2:6).^2;(2:6).^3;(2:6).^4;(2:6).^5;(2:6).^6;(2:6).^7;(2:6).^8;(2:6).^9];
% Sparse GH2:10

Egh_spar=100*abs((m_gh_spar(:,:)-repmat(truth,1,5))./repmat(truth,1,5));
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


figure(1)
hold on
% for ct = 1:2
    plot(Ngh(1,1:4),Egh(1,:),'ko--','linewidth',2,'MarkerSize',15)
    plot(Ngh(2,1:4),Egh(2,:),'kd--','linewidth',2,'MarkerSize',15)
    plot(Ngh(3,1:4),Egh(3,:),'ks--','linewidth',2,'MarkerSize',15)
    plot(Ngh(4,1:4),Egh(4,:),'kh--','linewidth',2,'MarkerSize',15)
    plot(Ngh(5,1:4),Egh(5,:),'k^--','linewidth',2,'MarkerSize',15)
    plot(Ngh(6,1:4),Egh(6,:),'k*--','linewidth',2,'MarkerSize',15)
    plot(Ngh(7,1:4),Egh(7,:),'kp--','linewidth',2,'MarkerSize',15)
    plot(Ngh(8,1:4),Egh(8,:),'k>--','linewidth',2,'MarkerSize',15)

%     set(gca,'YTick',0:0.25:3)
    set(gca,'XScale','log')
% end
legend('2D','3D','4D','5D','6D','7D','8D','9D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper


figure(2)
hold on
% for ct = 1:2
    plot(Ncut(1,:),Ecut(1,:),'ko--','linewidth',2,'MarkerSize',15)
    plot(Ncut(2,:),Ecut(2,:),'kd--','linewidth',2,'MarkerSize',15)
    plot(Ncut(3,:),Ecut(3,:),'ks--','linewidth',2,'MarkerSize',15)
    plot(Ncut(4,:),Ecut(4,:),'kh--','linewidth',2,'MarkerSize',15)
    plot(Ncut(5,:),Ecut(5,:),'k^--','linewidth',2,'MarkerSize',15)
plot(Ncut(6,:),Ecut(6,:),'k*--','linewidth',2,'MarkerSize',15)
plot(Ncut(7,:),Ecut(7,:),'kp--','linewidth',2,'MarkerSize',15)
plot(Ncut(8,:),Ecut(8,:),'k>--','linewidth',2,'MarkerSize',15)


%     set(gca,'YTick',0:0.25:3)
    set(gca,'XScale','log')
   % end
legend('2D','3D','4D','5D','6D','7D','8D','9D')
xlabel('No. of points (log-scale)')
ylabel('% relative error')
plot_prop_paper

figure(3)
hold on
% for ct = 1:2
    plot(N_spar(1,:),Egh_spar(1,:),'ko--','linewidth',2,'MarkerSize',15)
    plot(N_spar(2,:),Egh_spar(2,:),'kd--','linewidth',2,'MarkerSize',15)
    plot(N_spar(3,:),Egh_spar(3,:),'ks--','linewidth',2,'MarkerSize',15)
    plot(N_spar(4,:),Egh_spar(4,:),'kh--','linewidth',2,'MarkerSize',15)
    plot(N_spar(5,:),Egh_spar(5,:),'k^--','linewidth',2,'MarkerSize',15)
    plot(N_spar(6,:),Egh_spar(6,:),'k*--','linewidth',2,'MarkerSize',15)
    plot(N_spar(7,:),Egh_spar(7,:),'kp--','linewidth',2,'MarkerSize',15)
    plot(N_spar(8,:),Egh_spar(8,:),'k>--','linewidth',2,'MarkerSize',15)

    set(gca,'XScale','log')

legend('2D','3D','4D','5D','6D','7D','8D','9D')
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%CUT4-Uniform
clc
clear
for n=2:1:5
 [x,w]=smolyak_sparse_grid(n,5,'GLgn');
 [xint,wint] = GLeg_pts( 5*ones(1,n), -1*ones(1,n), ones(1,n));
 [Xcut,wcut]=uniform_sigma_pts(-1*ones(1,n), ones(1,n),8);
Nspar(n-1)=length(w);
Ngl(n-1)=length(wint);
Ncut(n-1)=length(wcut);
end
figure
bar3(2:1:5,[Ngl',Nspar',Ncut'])
set(gca,'XTickLabel',{'GLgn5','SPARSE-GLgn5','CUT8-U'})
view(119,24)
xlabel('Dimension')
ylabel('Method')
zlabel('Number of points')
plot_prop_paper
% figure
% bar3(2:1:4,[Ncut',Nspar'])
% set(gca,'XTickLabel',{'CUT6-U','SPARSE-GLgn4'})
% view(119,24)
% xlabel('Dimension')
% ylabel('Method')
% zlabel('Number of points')
% plot_prop_paper
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot to show the abs sum of weights
% for n=2:1:8
%     for k=2:1:5
% [x,w]=smolyak_sparse_grid(n,k,'GLgn');
% W(n-1,k-1)=sum(abs(w));
%     end
% end
% 
% plot(2:1:5,W(1,:),'o--',2:1:5,W(2,:),'s--',2:1:5,W(3,:),'d--',2:1:5,W(4,:),'p--',2:1:5,W(5,:),'*--',2:1:5,W(6,:),'^--',2:1:5,W(7,:),'h--','linewidth',2,'MarkerSize',15)
% legend('2D','3D','4D','5D','6D','7D','8D')
% xlabel('level k')
% ylabel('Sum of absolute value of weights')
% plot_prop_paper

figure
bar3(2:1:9,[NNgh',NNcut'])
set(gca,'XTickLabel',{'GLgn','CUT'})
xlabel('Dimension')
ylabel('Method')
zlabel('Number of points')
plot_prop_paper