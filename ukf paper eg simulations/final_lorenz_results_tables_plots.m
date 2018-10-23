
%% % Averaging over all the runs
est_fin_ukf=zeros(time.nSteps,3);
est_fin_ckf=zeros(time.nSteps,3);
est_fin_cut4=zeros(time.nSteps,3);
est_fin_cut6=zeros(time.nSteps,3);
est_fin_cut8=zeros(time.nSteps,3);
est_fin_gh=zeros(time.nSteps,3);
est_fin_ekf=zeros(time.nSteps,3);
est_fin_mupf=zeros(time.nSteps,3);
est_fin_mopf=zeros(time.nSteps,3);

avg_traj_ukf=zeros(time.nSteps,model.fn);
avg_traj_ckf=zeros(time.nSteps,model.fn);
avg_traj_cut4=zeros(time.nSteps,model.fn);
avg_traj_cut6=zeros(time.nSteps,model.fn);
avg_traj_cut8=zeros(time.nSteps,model.fn);
avg_traj_gh=zeros(time.nSteps,model.fn);
avg_traj_ekf=zeros(time.nSteps,model.fn);
avg_traj_mupf=zeros(time.nSteps,model.fn);
avg_traj_mopf=zeros(time.nSteps,model.fn);
avg_traj_mc=zeros(time.nSteps,model.fn);
avg_traj_meas=zeros(time.nSteps,model.hn);

Pest_fin_ukf=zeros(time.nSteps,model.hn^2);
Pest_fin_ckf=zeros(time.nSteps,model.hn^2);
Pest_fin_cut4=zeros(time.nSteps,model.hn^2);
Pest_fin_cut6=zeros(time.nSteps,model.hn^2);
Pest_fin_cut8=zeros(time.nSteps,model.hn^2);
Pest_fin_gh=zeros(time.nSteps,model.hn^2);
Pest_fin_pf=zeros(time.nSteps,model.hn^2);

% save(model.name,'model','time','pf','filter','PF_complete_data','avg_traj_ukf','avg_traj_ckf','avg_traj_cut4','avg_traj_cut6','avg_traj_cut8','avg_traj_gh','avg_traj_ekf','avg_traj_mopf','avg_traj_mupf','avg_traj_mc','est_fin_ukf','est_fin_ckf','est_fin_cut4','est_fin_cut6','est_fin_cut8','est_fin_gh','est_fin_ekf','est_fin_mupf','est_fin_mopf','xNNN_mc','YNNN_mc','xNNN_ukf','xNNN_ckf','xNNN_cut4','xNNN_cut6','xNNN_cut8','xNNN_gh','xNNN_ekf','xNNN_mupf','xNNN_mopf','PNNN_ukf','PNNN_ckf','PNNN_cut4','PNNN_cut6','PNNN_cut8','PNNN_gh','PNNN_ekf','PNNN_pf')
  for r=1:1:time.nSteps
        for c=1:1:model.fn
           est_fin_ukf(r,c)= sqrt(mean((xNNN_ukf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_ckf(r,c)= sqrt(mean((xNNN_ckf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_cut4(r,c)= sqrt(mean((xNNN_cut4(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_cut6(r,c)= sqrt(mean((xNNN_cut6(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_cut8(r,c)= sqrt(mean((xNNN_cut8(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_gh(r,c)= sqrt(mean((xNNN_gh(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_mupf(r,c)= sqrt(mean((xNNN_mupf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_mopf(r,c)= sqrt(mean((xNNN_mopf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
        end
        for c=1:1:model.fn^2
             Pest_fin_ukf(r,c)= sqrt(mean(PNNN_ukf(r,c,:)));
             Pest_fin_ckf(r,c)= sqrt(mean(PNNN_ckf(r,c,:)));
             Pest_fin_cut4(r,c)= sqrt(mean(PNNN_cut4(r,c,:)));
             Pest_fin_cut6(r,c)= sqrt(mean(PNNN_cut6(r,c,:)));
             Pest_fin_cut8(r,c)= sqrt(mean(PNNN_cut8(r,c,:)));
             Pest_fin_gh(r,c)= sqrt(mean(PNNN_gh(r,c,:)));
             Pest_fin_pf(r,c)= sqrt(mean(PNNN_pf(r,c,:)));
        end
  end

for r=1:1:time.nSteps
        for c=1:1:model.fn
avg_traj_ukf(r,c)=mean(xNNN_ukf(r,c,:));
avg_traj_ckf(r,c)=mean(xNNN_ckf(r,c,:));
avg_traj_cut4(r,c)=mean(xNNN_cut4(r,c,:));
avg_traj_cut6(r,c)=mean(xNNN_cut6(r,c,:));
avg_traj_cut8(r,c)=mean(xNNN_cut8(r,c,:));
avg_traj_gh(r,c)=mean(xNNN_gh(r,c,:));
avg_traj_ekf(r,c)=mean(xNNN_ekf(r,c,:));
avg_traj_mupf(r,c)=mean(xNNN_mupf(r,c,:));
avg_traj_mopf(r,c)=mean(xNNN_mopf(r,c,:));
avg_traj_mc(r,c)=mean(xNNN_mc(r,c,:));

        end
end

%  save(model.name,'model','time','pf','filter','PF_complete_data','avg_traj_ukf','avg_traj_ckf','avg_traj_cut4','avg_traj_cut6','avg_traj_cut8','avg_traj_gh','avg_traj_ekf','avg_traj_mopf','avg_traj_mupf','avg_traj_mc','est_fin_ukf','est_fin_ckf','est_fin_cut4','est_fin_cut6','est_fin_cut8','est_fin_gh','est_fin_ekf','est_fin_mupf','est_fin_mopf','xNNN_mc','YNNN_mc','xNNN_ukf','xNNN_ckf','xNNN_cut4','xNNN_cut6','xNNN_cut8','xNNN_gh','xNNN_ekf','xNNN_mupf','xNNN_mopf','PNNN_ukf','PNNN_ckf','PNNN_cut4','PNNN_cut6','PNNN_cut8','PNNN_gh','PNNN_ekf','PNNN_pf','Pest_fin_ukf','Pest_fin_ckf','Pest_fin_cut4','Pest_fin_cut6','Pest_fin_cut8','Pest_fin_gh','Pest_fin_pf')

  
% Plotting
 t=time.tspan;
%   figure(1)
% % subplot(2,1,1)
% plot(t,est_fin_mupf(:,1),'-.',t,est_fin_mopf(:,1),':',t,est_fin_ckf(:,1),t,est_fin_ukf(:,1),t,est_fin_cut4(:,1),t,est_fin_cut6(:,1),t,est_fin_cut8(:,1),t,est_fin_gh(:,1))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% % subplot(2,1,2)
% % plot(t,Pest_fin_pf(:,1),t,Pest_fin_ckf(:,1),t,Pest_fin_ukf(:,1),t,Pest_fin_cut4(:,1),t,Pest_fin_cut6(:,1),t,Pest_fin_cut8(:,1),t,Pest_fin_gh(:,1))
% % legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% 
% figure(2)
% % subplot(2,1,1)
% plot(t,est_fin_mupf(:,2),'-.',t,est_fin_mopf(:,2),':',t,est_fin_ckf(:,2),t,est_fin_ukf(:,2),t,est_fin_cut4(:,2),t,est_fin_cut6(:,2),t,est_fin_cut8(:,2),t,est_fin_gh(:,2))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% % subplot(2,1,2)
% % plot(t,Pest_fin_pf(:,7),t,Pest_fin_ckf(:,7),t,Pest_fin_ukf(:,7),t,Pest_fin_cut4(:,7),t,Pest_fin_cut6(:,7),t,Pest_fin_cut8(:,7),t,Pest_fin_gh(:,7))
% % legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% 
% figure(3)
% % subplot(2,1,1)
% plot(t,est_fin_mupf(:,3),'-.',t,est_fin_mopf(:,3),':',t,est_fin_ckf(:,3),t,est_fin_ukf(:,3),t,est_fin_cut4(:,3),t,est_fin_cut6(:,3),t,est_fin_cut8(:,3),t,est_fin_gh(:,3))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% % subplot(2,1,2)
% % plot(t,Pest_fin_pf(:,13),t,Pest_fin_ckf(:,13),t,Pest_fin_ukf(:,13),t,Pest_fin_cut4(:,13),t,Pest_fin_cut6(:,13),t,Pest_fin_cut8(:,13),t,Pest_fin_gh(:,13))
% % legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% 
% figure(4)
% % subplot(2,1,1)
% plot(t,est_fin_mupf(:,4),'-.',t,est_fin_mopf(:,4),':',t,est_fin_ckf(:,4),t,est_fin_ukf(:,4),t,est_fin_cut4(:,4),t,est_fin_cut6(:,4),t,est_fin_cut8(:,4),t,est_fin_gh(:,4))
% % legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% % subplot(2,1,2)
% % plot(t,Pest_fin_pf(:,19),t,Pest_fin_ckf(:,19),t,Pest_fin_ukf(:,19),t,Pest_fin_cut4(:,19),t,Pest_fin_cut6(:,19),t,Pest_fin_cut8(:,19),t,Pest_fin_gh(:,19))
% % legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% 
% figure(5)
% % subplot(2,1,1)
% plot(t,est_fin_mupf(:,5),'-.',t,est_fin_mopf(:,5),':',t,est_fin_ckf(:,5),t,est_fin_ukf(:,5),t,est_fin_cut4(:,5),t,est_fin_cut6(:,5),t,est_fin_cut8(:,5),t,est_fin_gh(:,5))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% % subplot(2,1,2)
% % plot(t,Pest_fin_pf(:,25),t,Pest_fin_ckf(:,25),t,Pest_fin_ukf(:,25),t,Pest_fin_cut4(:,25),t,Pest_fin_cut6(:,25),t,Pest_fin_cut8(:,25),t,Pest_fin_gh(:,25))
% % legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% 
% figure(6)
% plot(avg_traj_mupf(:,1),avg_traj_mupf(:,2),'-.',avg_traj_mopf(:,1),avg_traj_mopf(:,2),':',avg_traj_ckf(:,1),avg_traj_ckf(:,2),avg_traj_ukf(:,1),avg_traj_ukf(:,2),avg_traj_cut4(:,1),avg_traj_cut4(:,2),avg_traj_cut6(:,1),avg_traj_cut6(:,2),avg_traj_cut8(:,1),avg_traj_cut8(:,2),avg_traj_gh(:,1),avg_traj_gh(:,2),avg_traj_mc(:,1),avg_traj_mc(:,2),'+')
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh','MC')
% figure(7)
% plot(t,avg_traj_mupf(:,1),'-.',t,avg_traj_mopf(:,1),':',t,avg_traj_ckf(:,1),t,avg_traj_ukf(:,1),t,avg_traj_cut4(:,1),t,avg_traj_cut6(:,1),t,avg_traj_cut8(:,1),t,avg_traj_gh(:,1))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% matlabpool close


A=[sqrt(sum(sqrt(sum(est_fin_mupf.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_mopf.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_ckf.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_ukf.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_cut4.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_cut6.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_cut8.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_gh.^2/model.fn,2)).^2)/time.nSteps)] 
B=[max(max(est_fin_mupf)),max(max(est_fin_mopf)),max(max(est_fin_ckf)),max(max(est_fin_ukf)),max(max(est_fin_cut4)),max(max(est_fin_cut6)),max(max(est_fin_cut8)),max(max(est_fin_gh))]

% plot3(avg_traj_mupf(:,1),avg_traj_mupf(:,2),avg_traj_mupf(:,3),'-.',avg_traj_ckf(:,1),avg_traj_ckf(:,2),avg_traj_ckf(:,3),avg_traj_ukf(:,1),avg_traj_ukf(:,2),avg_traj_ukf(:,3),avg_traj_cut4(:,1),avg_traj_cut4(:,2),avg_traj_cut4(:,3),avg_traj_cut6(:,1),avg_traj_cut6(:,2),avg_traj_cut6(:,3),avg_traj_cut8(:,1),avg_traj_cut8(:,2),avg_traj_cut8(:,3))
% legend('PF-mu','ckf','ukf','cut4','cut6','cut8')