clc
clear
scrsz = get(0,'ScreenSize');
   [f,r]=meshgrid(1:1:10,1:2:20);
   Z1=zeros(size(f));
   Z2=zeros(size(f));
   Z3=zeros(size(f));
   Z4=zeros(size(f));
   Z5=zeros(size(f));
   Z6=zeros(size(f));
   
   ff=1;
   rr=1;
for para_set=1:1:100

    load(strcat('ckf_eg_run_nos_',num2str(10),'paraset_',num2str(para_set)));
    HHH=num2str(sqrt(CKFeg_simulation_para.R(2,2))*(180/pi));
    HHH(find(HHH=='.'))=',';
 
    if para_set==1;
    ir=sqrt(CKFeg_simulation_para.R(2,2));
    end
    
    Z1(ff,rr)=sqrt(sum(est_fin_ckf(:,3).^2));
    Z2(ff,rr)=sqrt(sum(est_fin_ukf(:,3).^2));
    Z3(ff,rr)=sqrt(sum(est_fin_cut4(:,3).^2));
    Z4(ff,rr)=sqrt(sum(est_fin_cut6(:,3).^2));
    Z5(ff,rr)=sqrt(sum(est_fin_cut8(:,3).^2));
    Z6(ff,rr)=sqrt(sum(est_fin_gh(:,3).^2));
    
    ff=ff+1;
    if ff==length(f)+1
        rr=rr+1;
        ff=1;
    end
    
% figure('Name','RMSE_pos(m)','Position',[15 1 1240 1000])
% plot(t,est_fin_ckf(:,1),t,est_fin_ukf(:,1),t,est_fin_cut4(:,1),t,est_fin_cut6(:,1),t,est_fin_cut8(:,1),t,est_fin_gh(:,1),'LineWidth',1.5)
% legend('ckf','ukf','cut4','cut6','cut8','gh')
% xlabel('t')
% ylabel('RMSE_{pos}(m)')
% set(gca,'FontSize',14)
% title(strcat('Meas-updt-freq = ',num2str(CKFeg_simulation_para.freq),':  :Std-dev-bearing = ',num2str(sqrt(CKFeg_simulation_para.R(2,2))*(180/pi)),' deg'))
% h = get(gca, 'title');
% k = get(gca, 'xlabel');
% l = get(gca, 'ylabel');
% set(h, 'FontName', 'Helvetica', 'FontSize', 14)
% set(k, 'FontName', 'Helvetica', 'FontSize', 16)
% set(l, 'FontName', 'Helvetica', 'FontSize', 16)
% saveas(gcf,strcat('RMSE_pos',' f = ',num2str(CKFeg_simulation_para.freq),' th = ',HHH),'jpg') 
% 
% close
% 
% figure('Name','RMSE_vel(m)','Position',[15 1 1240 1000])
% plot(t,est_fin_ckf(:,2),t,est_fin_ukf(:,2),t,est_fin_cut4(:,2),t,est_fin_cut6(:,2),t,est_fin_cut8(:,2),t,est_fin_gh(:,2),'LineWidth',1.5)
% legend('ckf','ukf','cut4','cut6','cut8','gh')
% xlabel('t')
% ylabel('RMSE_{vel}(ms^{-1})')
% set(gca,'FontSize',14)
% title(strcat('Meas-updt-freq = ',num2str(CKFeg_simulation_para.freq),':  :Std-dev-bearing = ',num2str(sqrt(CKFeg_simulation_para.R(2,2))*(180/pi)),' deg'))
% h = get(gca, 'title');
% k = get(gca, 'xlabel');
% l = get(gca, 'ylabel');
% set(h, 'FontName', 'Helvetica', 'FontSize', 14)
% set(k, 'FontName', 'Helvetica', 'FontSize', 16)
% set(l, 'FontName', 'Helvetica', 'FontSize', 16)
% saveas(gcf,strcat('RMSE_vel',' f = ',num2str(CKFeg_simulation_para.freq),' th = ',HHH),'jpg') 
% 
% close
% 
% figure('Name','RMSE_omega(m)','Position',[15 1 1240 1000])
% plot(t,(180/pi)*est_fin_ckf(:,3),t,(180/pi)*est_fin_ukf(:,3),t,(180/pi)*est_fin_cut4(:,3),t,(180/pi)*est_fin_cut6(:,3),t,(180/pi)*est_fin_cut8(:,3),t,(180/pi)*est_fin_gh(:,3),'LineWidth',1.5)
% legend('ckf','ukf','cut4','cut6','cut8','gh')
% xlabel('t')
% ylabel('RMSE_{omega}(deg)')
% set(gca,'FontSize',14)
% title(strcat('Meas-updt-freq = ',num2str(CKFeg_simulation_para.freq),':  :Std-dev-bearing = ',num2str(sqrt(CKFeg_simulation_para.R(2,2))*(180/pi)),' deg'))
% h = get(gca, 'title');
% k = get(gca, 'xlabel');
% l = get(gca, 'ylabel');
% set(h, 'FontName', 'Helvetica', 'FontSize', 14)
% set(k, 'FontName', 'Helvetica', 'FontSize', 16)
% set(l, 'FontName', 'Helvetica', 'FontSize', 16)
% saveas(gcf,strcat('RMSE_omega',' f = ',num2str(CKFeg_simulation_para.freq),' th = ',HHH),'jpg') 
% close
% 
% figure('Name','Traj-Sensor-Meas','Position',[15 1 1240 1000])
% plot(x100_ckf(:,1,end),x100_ckf(:,3,end),x100_ukf(:,1,end),x100_ukf(:,3,end),x100_cut4(:,1,end),x100_cut4(:,3,end),x100_cut6(:,1,end),x100_cut6(:,3,end),x100_cut8(:,1,end),x100_cut8(:,3,end),x100_gh(:,1,end),x100_gh(:,3,end),x_mc(:,1),x_mc(:,3),ym(:,1).*cos(ym(:,2)),ym(:,1).*sin(ym(:,2)),'k*','LineWidth',1.5)
% legend('ckf','ukf','cut4','cut6','cut8','gh','MC','sensor')
% xlabel('x')
% ylabel('y')
% set(gca,'FontSize',14)
% title(strcat('Meas-updt-freq = ',num2str(CKFeg_simulation_para.freq),':  :Std-dev-bearing = ',num2str(sqrt(CKFeg_simulation_para.R(2,2))*(180/pi)),' deg'))
% h = get(gca, 'title');
% k = get(gca, 'xlabel');
% l = get(gca, 'ylabel');
% set(h, 'FontName', 'Helvetica', 'FontSize', 14)
% set(k, 'FontName', 'Helvetica', 'FontSize', 16)
% set(l, 'FontName', 'Helvetica', 'FontSize', 16)
% saveas(gcf,strcat('Trajectories and Sensor Meas',' f = ',num2str(CKFeg_simulation_para.freq),' th = ',HHH),'jpg') 
% close
end

for rr=1:1:10
    figure(rr)
    plot(1:1:10,Z1(1:1:10,rr),'-ks',1:1:10,Z2(1:1:10,rr),'-ko',1:1:10,Z3(1:1:10,rr),'-k*',1:1:10,Z4(1:1:10,rr),'-kd',1:1:10,Z5(1:1:10,rr),'-k^',1:1:10,Z6(1:1:10,rr),'-kp')
    xlabel('Time gap between meas updates')
    ylabel('2-norm of the RMSE error in omg')
    title(['Bearing std dev ',num2str((2*rr-1)*ir*180/pi), 'degrees'])
    legend('ckf','ukf','cut4','cut6','cut8','gh')
    hh=num2str((2*rr-1)*ir*180/pi);
    hh(find(hh=='.'))=',';
    saveas(gcf,strcat('2 norm error in omg',' std dev = ',hh),'jpg') 
close
end
regexprep('abcd', 'd', '2')
% figure(1)
% surf(f,r,Z1)
% xlabel('freq')
% ylabel('mag std dev bearing')
% title('ckf')
% axis([1 10 1 20 0 12000])
% 
% figure(2)
% surf(f,r,Z2)
% xlabel('freq')
% ylabel('mag std dev bearing')
% title('ukf')
% axis([1 10 1 20 0 12000])
% 
% figure(3)
% surf(f,r,Z3)
% xlabel('freq')
% ylabel('mag std dev bearing')
% title('cut4')
% axis([1 10 1 20 0 12000])
% 
% figure(4)
% surf(f,r,Z4)
% xlabel('freq')
% ylabel('mag std dev bearing')
% title('cut6')
% axis([1 10 1 20 0 12000])
% 
% figure(5)
% surf(f,r,Z5)
% xlabel('freq')
% ylabel('mag std dev bearing')
% title('cut8')
% axis([1 10 1 20 0 12000])
% 
% figure(6)
% surf(f,r,Z6)
% xlabel('freq')
% ylabel('mag std dev bearing')
% title('gh2')
% axis([1 10 1 20 0 12000])