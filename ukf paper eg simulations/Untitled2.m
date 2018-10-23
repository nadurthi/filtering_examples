pp=1;
plot(t,est_fin_ckf(:,1),'rs-',t,est_fin_ukf(:,1),'bo-',t,est_fin_cut6(:,1),'m*-',t,est_fin_cut8(:,pp),'kd-',t,est_fin_mopf(:,pp),'g--',t,est_fin_mupf(:,pp),'g','linewidth',2,'MarkerSize',10)
legend('ckf','ukf','cut6','cut8','pf-mo','pf-mu')
ylabel('RMSE in x_1')
xlabel('t')
plot_prop_paper

figure
pp=2;
plot(t,est_fin_ckf(:,pp),'rs-',t,est_fin_ukf(:,pp),'bo-',t,est_fin_cut6(:,pp),'m*-',t,est_fin_cut8(:,pp),'kd-',t,est_fin_mopf(:,pp),'g--',t,est_fin_mupf(:,pp),'g','linewidth',2,'MarkerSize',10)
legend('ckf','ukf','cut6','cut8','pf-mo','pf-mu')
ylabel('RMSE in x_2')
xlabel('t')
plot_prop_paper

figure
pp=3;
plot(t,est_fin_ckf(:,pp),'rs-',t,est_fin_ukf(:,pp),'bo-',t,est_fin_cut6(:,pp),'m*-',t,est_fin_cut8(:,pp),'kd-',t,est_fin_mopf(:,pp),'g--',t,est_fin_mupf(:,pp),'g','linewidth',2,'MarkerSize',10)
legend('ckf','ukf','cut6','cut8','pf-mo','pf-mu')
ylabel('RMSE in x_3')
xlabel('t')
plot_prop_paper

figure
pp=4;
plot(t,est_fin_ckf(:,pp),'rs-',t,est_fin_ukf(:,pp),'bo-',t,est_fin_cut6(:,pp),'m*-',t,est_fin_cut8(:,pp),'kd-',t,est_fin_mopf(:,pp),'g--',t,est_fin_mupf(:,pp),'g','linewidth',2,'MarkerSize',10)
legend('ckf','ukf','cut6','cut8','pf-mo','pf-mu')
ylabel('RMSE in \sigma')
xlabel('t')
plot_prop_paper

figure
pp=5;
plot(t,est_fin_ckf(:,pp),'rs-',t,est_fin_ukf(:,pp),'bo-',t,est_fin_cut6(:,pp),'m*-',t,est_fin_cut8(:,pp),'kd-',t,est_fin_mopf(:,pp),'g--',t,est_fin_mupf(:,pp),'g','linewidth',2,'MarkerSize',10)
legend('ckf','ukf','cut6','cut8','pf-mo','pf-mu')
ylabel('RMSE in \rho')
xlabel('t')
plot_prop_paper