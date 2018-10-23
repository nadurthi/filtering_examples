%% Plot results of ukf paper

clc
clear
scrsz = get(0,'ScreenSize');
    
 
for f=2:40:360
    for rr=1:1:4

%     load(strcat('ukf_set2_eg_run_nos_',num2str(1),'_paraset_no_',num2str(paraset)));
load(strcat('ukf_set3_eg_run_nos_',num2str(1),'_paraset_no_',num2str(rr),num2str(f)));
    HHH=num2str(sqrt(ukf_eg_para_set.R(2,2))*(180/pi));
    HHH(find(HHH=='.'))=',';
    %%%%%%%%%%%%%%%%   x1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','RMSE_x1(km)','Position',[15 1 1240 1000])
plot(t,est_fin_ckf(:,1),t,est_fin_ukf(:,1),t,est_fin_cut4(:,1),t,est_fin_cut6(:,1),t,est_fin_cut8(:,1),t,est_fin_gh(:,1),'LineWidth',1.5)
legend('ckf','ukf','cut4','cut6','cut8','gh')
xlabel('t')
ylabel('RMSE_{x_1}(km)')
set(gca,'FontSize',14)
title(strcat('Meas-updt-freq = ',num2str(ukf_eg_para_set.freq),':  :Std-dev-bearing = ',num2str(sqrt(ukf_eg_para_set.R(2,2))*(180/pi)),' deg'))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
saveas(gcf,strcat('1New RMSE_x1',' f = ',num2str(ukf_eg_para_set.freq),' th = ',HHH),'jpg') 

close

figure('Name','Pcov_x1(km)','Position',[15 1 1240 1000])
plot(t,Pest_fin_ckf(:,1),t,Pest_fin_ukf(:,1),t,Pest_fin_cut4(:,1),t,Pest_fin_cut6(:,1),t,Pest_fin_cut8(:,1),t,Pest_fin_gh(:,1),'LineWidth',1.5)
legend('ckf','ukf','cut4','cut6','cut8','gh')
xlabel('t')
ylabel('Pcov_{x_1}^{1/2}(km)')
set(gca,'FontSize',14)
title(strcat('Meas-updt-freq = ',num2str(ukf_eg_para_set.freq),':  :Std-dev-bearing = ',num2str(sqrt(ukf_eg_para_set.R(2,2))*(180/pi)),' deg'))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
saveas(gcf,strcat('New Pcov_x1',' f = ',num2str(ukf_eg_para_set.freq),' th = ',HHH),'jpg') 

close
%%%%%%%%%%%%%%%%   x3   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','RMSE_x3(kmps)','Position',[15 1 1240 1000])
plot(t,est_fin_ckf(:,3),t,est_fin_ukf(:,3),t,est_fin_cut4(:,3),t,est_fin_cut6(:,3),t,est_fin_cut8(:,3),t,est_fin_gh(:,3),'LineWidth',1.5)
legend('ckf','ukf','cut4','cut6','cut8','gh')
xlabel('t')
ylabel('RMSE_{x_3}(km/s)')
set(gca,'FontSize',14)
title(strcat('Meas-updt-freq = ',num2str(ukf_eg_para_set.freq),':  :Std-dev-bearing = ',num2str(sqrt(ukf_eg_para_set.R(2,2))*(180/pi)),' deg'))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
saveas(gcf,strcat('2New RMSE_x3',' f = ',num2str(ukf_eg_para_set.freq),' th = ',HHH),'jpg') 

close

figure('Name','Pcov_x3(kmps)','Position',[15 1 1240 1000])
plot(t,Pest_fin_ckf(:,13),t,Pest_fin_ukf(:,13),t,Pest_fin_cut4(:,13),t,Pest_fin_cut6(:,13),t,Pest_fin_cut8(:,13),t,Pest_fin_gh(:,13),'LineWidth',1.5)
 legend('ckf','ukf','cut4','cut6','cut8','gh')
xlabel('t')
ylabel('Pcov_{x_3}^{1/2}(km/s)')
set(gca,'FontSize',14)
title(strcat('Meas-updt-freq = ',num2str(ukf_eg_para_set.freq),':  :Std-dev-bearing = ',num2str(sqrt(ukf_eg_para_set.R(2,2))*(180/pi)),' deg'))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
saveas(gcf,strcat('New Pcov_x3',' f = ',num2str(ukf_eg_para_set.freq),' th = ',HHH),'jpg') 

close
%%%%%%%%%%%%%%%%   x5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','RMSE_x5','Position',[15 1 1240 1000])
plot(t,est_fin_ckf(:,5),t,est_fin_ukf(:,5),t,est_fin_cut4(:,5),t,est_fin_cut6(:,5),t,est_fin_cut8(:,5),t,est_fin_gh(:,5),'LineWidth',1.5)
legend('ckf','ukf','cut4','cut6','cut8','gh')
xlabel('t')
ylabel('RMSE_{x_5}')
set(gca,'FontSize',14)
title(strcat('Meas-updt-freq = ',num2str(ukf_eg_para_set.freq),':  :Std-dev-bearing = ',num2str(sqrt(ukf_eg_para_set.R(2,2))*(180/pi)),' deg'))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
saveas(gcf,strcat('3New RMSE_x5',' f = ',num2str(ukf_eg_para_set.freq),' th = ',HHH),'jpg') 

close

figure('Name','Pcov_x5','Position',[15 1 1240 1000])
plot(t,Pest_fin_ckf(:,25),t,Pest_fin_ukf(:,25),t,Pest_fin_cut4(:,25),t,Pest_fin_cut6(:,25),t,Pest_fin_cut8(:,25),t,Pest_fin_gh(:,25),'LineWidth',1.5)
legend('ckf','ukf','cut4','cut6','cut8','gh')
xlabel('t')
ylabel('Pcov_{x_5}^{1/2}')
set(gca,'FontSize',14)
title(strcat('Meas-updt-freq = ',num2str(ukf_eg_para_set.freq),':  :Std-dev-bearing = ',num2str(sqrt(ukf_eg_para_set.R(2,2))*(180/pi)),' deg'))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
saveas(gcf,strcat('New Pcov_x5',' f = ',num2str(ukf_eg_para_set.freq),' th = ',HHH),'jpg') 
close

figure('Name','traj','Position',[15 1 1240 1000])
p=-15:0.1:15;
xear=6374*cosd(p);
year=6374*sind(p);
plot(x100_ckf(:,1,1),x100_ckf(:,2,1),x100_ukf(:,1,1),x100_ukf(:,2,1),x100_cut4(:,1,1),x100_cut4(:,2,1),x100_cut6(:,1,1),x100_cut6(:,2,1),x100_cut8(:,1,1),x100_cut8(:,2,1),x100_gh(:,1,1),x100_gh(:,2,1),xear,year,'.','LineWidth',1.5)
legend('ckf','ukf','cut4','cut6','cut8','gh')
axis([6350 6500 -200 500])
xlabel('x_1(km)')
ylabel('x_2(km)')
set(gca,'FontSize',14)
title(strcat('Meas-updt-freq = ',num2str(ukf_eg_para_set.freq),':  :Std-dev-bearing = ',num2str(sqrt(ukf_eg_para_set.R(2,2))*(180/pi)),' deg'))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
saveas(gcf,strcat('New traj',' f = ',num2str(ukf_eg_para_set.freq),' th = ',HHH),'jpg') 
close
    end
end