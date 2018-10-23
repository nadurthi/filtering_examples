%%  19th feb 2012 sim of ckf
time.t0 = 0;
time.dt = 1;
time.tf = 100;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);

   [f,r]=meshgrid(1:1:10,1:2:20);
   Z1=zeros([size(f),7]);
   Z2=zeros([size(f),7]);
   Z3=zeros([size(f),7]);
   Z4=zeros([size(f),7]);
   Z5=zeros([size(f),7]);
   Z6=zeros([size(f),7]);
   Z7=zeros([size(f),7]);
   Z8=zeros([size(f),7]);

for i=1:1:10
    for j=1:1:10
[i,j]
       load(strcat('FEB19_CKF_allfilters_data_set_no_',num2str(i),'_',num2str(j),'_200'));
 est_fin_ukf=zeros(time.nSteps,7);
  est_fin_ckf=zeros(time.nSteps,7);  
   est_fin_cut4=zeros(time.nSteps,7);  
    est_fin_cut6=zeros(time.nSteps,7);  
     est_fin_cut8=zeros(time.nSteps,7);  
      est_fin_gh=zeros(time.nSteps,7); 
       est_fin_mupf=zeros(time.nSteps,7);  
        est_fin_mopf=zeros(time.nSteps,7); 
        
         for c=1:1:5  
      parfor r=1:time.nSteps
        
           est_fin_ukf(r,c)= sqrt(mean((xNNN_ukf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_ckf(r,c)= sqrt(mean((xNNN_ckf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_cut4(r,c)= sqrt(mean((xNNN_cut4(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_cut6(r,c)= sqrt(mean((xNNN_cut6(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_cut8(r,c)= sqrt(mean((xNNN_cut8(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_gh(r,c)= sqrt(mean((xNNN_gh(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_mupf(r,c)= sqrt(mean((xNNN_mupf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
           est_fin_mopf(r,c)= sqrt(mean((xNNN_mopf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
        end
      end
      for r=1:time.nSteps
           est_fin_ukf(r,6)= sqrt(mean((xNNN_ukf(r,1,:)-xNNN_mc(r,1,:)).^2+(xNNN_ukf(r,3,:)-xNNN_mc(r,3,:)).^2))^2;
           est_fin_ckf(r,6)= sqrt(mean((xNNN_ckf(r,1,:)-xNNN_mc(r,1,:)).^2+(xNNN_ckf(r,3,:)-xNNN_mc(r,3,:)).^2))^2;
           est_fin_cut4(r,6)= sqrt(mean((xNNN_cut4(r,1,:)-xNNN_mc(r,1,:)).^2+(xNNN_cut4(r,3,:)-xNNN_mc(r,3,:)).^2))^2;
           est_fin_cut6(r,6)= sqrt(mean((xNNN_cut6(r,1,:)-xNNN_mc(r,1,:)).^2+(xNNN_cut6(r,3,:)-xNNN_mc(r,3,:)).^2))^2;
           est_fin_cut8(r,6)= sqrt(mean((xNNN_cut8(r,1,:)-xNNN_mc(r,1,:)).^2+(xNNN_cut8(r,3,:)-xNNN_mc(r,3,:)).^2))^2;
           est_fin_gh(r,6)= sqrt(mean((xNNN_gh(r,1,:)-xNNN_mc(r,1,:)).^2+(xNNN_gh(r,3,:)-xNNN_mc(r,3,:)).^2))^2;
           est_fin_mupf(r,6)= sqrt(mean((xNNN_mupf(r,1,:)-xNNN_mc(r,1,:)).^2+(xNNN_mupf(r,3,:)-xNNN_mc(r,3,:)).^2))^2;
           est_fin_mopf(r,6)= sqrt(mean((xNNN_mopf(r,1,:)-xNNN_mc(r,1,:)).^2+(xNNN_mopf(r,3,:)-xNNN_mc(r,3,:)).^2))^2;
      
           est_fin_ukf(r,7)= sqrt(mean((xNNN_ukf(r,2,:)-xNNN_mc(r,2,:)).^2+(xNNN_ukf(r,4,:)-xNNN_mc(r,4,:)).^2))^2;
           est_fin_ckf(r,7)= sqrt(mean((xNNN_ckf(r,2,:)-xNNN_mc(r,2,:)).^2+(xNNN_ckf(r,4,:)-xNNN_mc(r,4,:)).^2))^2;
           est_fin_cut4(r,7)= sqrt(mean((xNNN_cut4(r,2,:)-xNNN_mc(r,2,:)).^2+(xNNN_cut4(r,4,:)-xNNN_mc(r,4,:)).^2))^2;
           est_fin_cut6(r,7)= sqrt(mean((xNNN_cut6(r,2,:)-xNNN_mc(r,2,:)).^2+(xNNN_cut6(r,4,:)-xNNN_mc(r,4,:)).^2))^2;
           est_fin_cut8(r,7)= sqrt(mean((xNNN_cut8(r,2,:)-xNNN_mc(r,2,:)).^2+(xNNN_cut8(r,4,:)-xNNN_mc(r,4,:)).^2))^2;
           est_fin_gh(r,7)= sqrt(mean((xNNN_gh(r,2,:)-xNNN_mc(r,2,:)).^2+(xNNN_gh(r,4,:)-xNNN_mc(r,4,:)).^2))^2;
           est_fin_mupf(r,7)= sqrt(mean((xNNN_mupf(r,2,:)-xNNN_mc(r,2,:)).^2+(xNNN_mupf(r,4,:)-xNNN_mc(r,4,:)).^2))^2;
           est_fin_mopf(r,7)= sqrt(mean((xNNN_mopf(r,2,:)-xNNN_mc(r,2,:)).^2+(xNNN_mopf(r,4,:)-xNNN_mc(r,4,:)).^2))^2;
      end
      
      for pp=1:1:7
    Z1(i,j,pp)=sqrt(sum(est_fin_ckf(:,pp).^2));
    Z2(i,j,pp)=sqrt(sum(est_fin_ukf(:,pp).^2));
    Z3(i,j,pp)=sqrt(sum(est_fin_cut4(:,pp).^2));
    Z4(i,j,pp)=sqrt(sum(est_fin_cut6(:,pp).^2));
    Z5(i,j,pp)=sqrt(sum(est_fin_cut8(:,pp).^2));
    Z6(i,j,pp)=sqrt(sum(est_fin_gh(:,pp).^2));
    Z7(i,j,pp)=sqrt(sum(est_fin_mopf(:,pp).^2));
    Z8(i,j,pp)=sqrt(sum(est_fin_mupf(:,pp).^2));
      end
    %  t=time.tspan;
%   figure(1)
% % subplot(2,1,1)
% plot(t,est_fin_mupf(:,1),'-.',t,est_fin_mopf(:,1),':',t,est_fin_ckf(:,1),t,est_fin_ukf(:,1),t,est_fin_cut4(:,1),t,est_fin_cut6(:,1),t,est_fin_cut8(:,1),t,est_fin_gh(:,1))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% % subplot(2,1,2)
% % plot(t,Pest_fin_pf(:,1),t,Pest_fin_ckf(:,1),t,Pest_fin_ukf(:,1),t,Pest_fin_cut4(:,1),t,Pest_fin_cut6(:,1),t,Pest_fin_cut8(:,1),t,Pest_fin_gh(:,1))
% legend('PF','ckf','ukf','cut4','cut6','cut8','gh')

% figure(2)
% % subplot(2,1,1)
% plot(t,est_fin_mupf(:,2),'-.',t,est_fin_mopf(:,2),':',t,est_fin_ckf(:,2),t,est_fin_ukf(:,2),t,est_fin_cut4(:,2),t,est_fin_cut6(:,2),t,est_fin_cut8(:,2),t,est_fin_gh(:,2))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% % subplot(2,1,2)
% % plot(t,Pest_fin_pf(:,7),t,Pest_fin_ckf(:,7),t,Pest_fin_ukf(:,7),t,Pest_fin_cut4(:,7),t,Pest_fin_cut6(:,7),t,Pest_fin_cut8(:,7),t,Pest_fin_gh(:,7))
% % legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% % 
% figure(3)
% % subplot(2,1,1)
% plot(t,est_fin_mupf(:,3),'-.',t,est_fin_mopf(:,3),':',t,est_fin_ckf(:,3),t,est_fin_ukf(:,3),t,est_fin_cut4(:,3),t,est_fin_cut6(:,3),t,est_fin_cut8(:,3),t,est_fin_gh(:,3))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% % subplot(2,1,2)
% % plot(t,Pest_fin_pf(:,13),t,Pest_fin_ckf(:,13),t,Pest_fin_ukf(:,13),t,Pest_fin_cut4(:,13),t,Pest_fin_cut6(:,13),t,Pest_fin_cut8(:,13),t,Pest_fin_gh(:,13))
% % legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% % 
% figure(4)
% % subplot(2,1,1)
% plot(t,est_fin_mupf(:,4),'-.',t,est_fin_mopf(:,4),':',t,est_fin_ckf(:,4),t,est_fin_ukf(:,4),t,est_fin_cut4(:,4),t,est_fin_cut6(:,4),t,est_fin_cut8(:,4),t,est_fin_gh(:,4))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% % subplot(2,1,2)
% % plot(t,Pest_fin_pf(:,19),t,Pest_fin_ckf(:,19),t,Pest_fin_ukf(:,19),t,Pest_fin_cut4(:,19),t,Pest_fin_cut6(:,19),t,Pest_fin_cut8(:,19),t,Pest_fin_gh(:,19))
% % legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% % 
% figure(5)
% % subplot(2,1,1)
% plot(t,est_fin_mupf(:,5),'-.',t,est_fin_mopf(:,5),':',t,est_fin_ckf(:,5),t,est_fin_ukf(:,5),t,est_fin_cut4(:,5),t,est_fin_cut6(:,5),t,est_fin_cut8(:,5),t,est_fin_gh(:,5))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% % subplot(2,1,2)
% % plot(t,Pest_fin_pf(:,25),t,Pest_fin_ckf(:,25),t,Pest_fin_ukf(:,25),t,Pest_fin_cut4(:,25),t,Pest_fin_cut6(:,25),t,Pest_fin_cut8(:,25),t,Pest_fin_gh(:,25))
% % legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% % 
% % figure(6)
% % plot(avg_traj_mupf(:,1),avg_traj_mupf(:,2),'-.',avg_traj_mopf(:,1),avg_traj_mopf(:,2),':',avg_traj_ckf(:,1),avg_traj_ckf(:,2),avg_traj_ukf(:,1),avg_traj_ukf(:,2),avg_traj_cut4(:,1),avg_traj_cut4(:,2),avg_traj_cut6(:,1),avg_traj_cut6(:,2),avg_traj_cut8(:,1),avg_traj_cut8(:,2),avg_traj_gh(:,1),avg_traj_gh(:,2),avg_traj_mc(:,1),avg_traj_mc(:,2),'+')
% % legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh','MC')
% %  figure(7)
% % plot(t,avg_traj_mupf(:,1),'-.',t,avg_traj_mopf(:,1),':',t,avg_traj_ckf(:,1),t,avg_traj_ukf(:,1),t,avg_traj_cut4(:,1),t,avg_traj_cut6(:,1),t,avg_traj_cut8(:,1),t,avg_traj_gh(:,1))
% % legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
%          
%      pause
      
     
    end
end
save('CKF_post_proc_results_save','Z1','Z2','Z3','Z4','Z5','Z6','Z7','Z8')
ir=sqrt(10)*1e-3;
pp=6;
for rr=1:1:10
    figure(rr)
   %    plot(1:1:10,Z1(1:1:10,rr,pp),'-ks',1:1:10,Z2(1:1:10,rr,pp),'-ko',1:1:10,Z3(1:1:10,rr,pp),'-k*',1:1:10,Z4(1:1:10,rr,pp),'-kd',1:1:10,Z5(1:1:10,rr,pp),'-k^',1:1:10,Z6(1:1:10,rr,pp),'-kp',1:1:10,Z7(1:1:10,rr,pp),'-k<',1:1:10,Z8(1:1:10,rr,pp),'-k>','LineWidth',2,'MarkerSize',15)
   plot(1:1:10,Z1(1:1:10,rr,pp),'-ks',1:1:10,Z2(1:1:10,rr,pp),'-k*',1:1:10,Z4(1:1:10,rr,pp),'-ko',1:1:10,Z7(1:1:10,rr,pp),'-k<',1:1:10,Z8(1:1:10,rr,pp),'-k>','LineWidth',2,'MarkerSize',17)
    xlabel('Time gap between meas updates')
    if pp==5
    ylabel('2-norm of the RMSE error in omg')
    elseif pp==6
        ylabel('2-norm of the RMSE error in pos')
    elseif pp==7
        ylabel('2-norm of the RMSE error in vel')
    end
    title(['Bearing std dev ',num2str((2*rr-1)*ir*180/pi), 'degrees'])
%     legend('ckf','ukf','cut4','cut6','cut8','gh','PF-Mode','PF-Mean','Location','NorthWest')
legend('ckf','ukf','cut6','PF-Mode','PF-Mean','Location','NorthWest')
    hh=num2str((2*rr-1)*ir*180/pi);
    hh(find(hh=='.'))='_';
    plot_prop_paper
    if pp==5
    saveas(gcf,strcat('2_norm_error3_in_omg_std_dev_',hh,'.eps'),'psc2') 
    elseif pp==6
      saveas(gcf,strcat('2_norm_error3_in_pos_std_dev_',hh,'.eps'),'psc2')   
    elseif pp==7
    saveas(gcf,strcat('2_norm_error3_in_vel_std_dev_',hh,'.eps'),'psc2') 
    end
    
close
end