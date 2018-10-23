%%ckf_paper_all_filters_call
clear all
clc 
close all
% global T
%   matlabpool open 3
warning off;
%% ------------------------------------------------------------------------
% time
%--------------------------------------------------------------------------
time.t0 = 0;
time.dt = 1;
time.tf = 100;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);
%__________________________________________________________________________

%% ------------------------------------------------------------------------
% model
T=time.dt;
    

    omg=-3*pi/180;
    q1=0.1;
    q2=1.75*10^(-4);
    
    sigr=10;
    sigth=0.1833*(pi/180);

    M=[T^3/3,T^2/2;T^2/2,T];
    Q=blkdiag(q1*M,q1*M,q2*T);
    
    R=diag([sigr^2,sigth^2]);
    
model.fn = 5;               % state space dimensionality
model.fx = @CKF_eg_dyn_disc;
model.fx_jac=@KIRB_CT_eg_dyn_jac_disc;

model.hn = 2;               % measurement dimensionality
model.hx =@CKF_eg_meas_disc;
model.hx_jac=@KIRB_eg_meas_jac_disc;

model.Q = Q;
model.sQ =sqrtm(Q);
model.R = R;
model.sR=sqrtm(R);
model.Qtruth=1e-200; % noise to generate the truth
model.Qt_sq=1e-200; % the corresponding sqrt


model.x0tr=[1000,300,1000,0,omg]'; % initial conditons
model.P0tr=diag([100,10,100,10,100e-6]);

%% ------------------------------------------------------------------------
% particle filter settings
%--------------------------------------------------------------------------
pf.no_particles = 5000;
pf.no_bins = 1000;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;                  % treshold for the number of effective samples

%% --------------------------------------------------------------
%filter props
% x0f=[6500.4,349.14,-1.8093,-6.7967,0.6932]';
% P0f=diag([1e-1,1e-1,2e-2,2e-2,2]);

filter.paras_ukf_kappa=1;
filter.paras_gh_pts=5;

filter.freq=1; %% This is actually the number of 'dt' steps after which a meas updt is done

% filter.x0_filt_start=x0f;
filter.P0_filt_start=diag([100,10,100,10,100e-6]);

filter.save_pf_data='false';

%switch on filters
filter.EKF='false';
filter.CKF='true';
filter.UKF='true';
filter.KF='false';
filter.GHKF='true';
filter.PF='true';
filter.CUT4KF='true';
filter.CUT6KF='true';
filter.CUT8KF='true';
filter.truth=[];
filter.ymeas=[];
%% Allocate space for NNN no. of iterations
NNN=2;

xNNN_mc=zeros(time.nSteps,model.fn,NNN);
YNNN_mc=zeros(time.nSteps,model.hn,NNN);
xNNN_ukf=zeros(time.nSteps,model.fn,NNN);
xNNN_ckf=zeros(time.nSteps,model.fn,NNN);
xNNN_cut4=zeros(time.nSteps,model.fn,NNN);
xNNN_cut6=zeros(time.nSteps,model.fn,NNN);
xNNN_cut8=zeros(time.nSteps,model.fn,NNN);
xNNN_gh=zeros(time.nSteps,model.fn,NNN);
xNNN_ekf=zeros(time.nSteps,model.fn,NNN);
xNNN_mupf=zeros(time.nSteps,model.fn,NNN);
xNNN_mopf=zeros(time.nSteps,model.fn,NNN);
PNNN_ukf=zeros(time.nSteps,model.fn^2,NNN);
PNNN_ckf=zeros(time.nSteps,model.fn^2,NNN);
PNNN_cut4=zeros(time.nSteps,model.fn^2,NNN);
PNNN_cut6=zeros(time.nSteps,model.fn^2,NNN);
PNNN_cut8=zeros(time.nSteps,model.fn^2,NNN);
PNNN_gh=zeros(time.nSteps,model.fn^2,NNN);
PNNN_ekf=zeros(time.nSteps,model.fn^2,NNN);
PNNN_pf=zeros(time.nSteps,model.fn^2,NNN);
PF_complete_data=cell([1,NNN]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the NNN data sets in a cell to sent to parfor lloop
[f,rr]=meshgrid(1:1:10,1:2:20);
    data_set=cell([10,10,NNN,4]);
%    nnn=1; 
for i=1:1:10
    for j=1:1:10
     for k=1:1:NNN
       filter.freq=f(i,j);
       model.R=diag([sigr^2,(rr(i,j)*sigth)^2]);
    [t,x_mc]=ode45_discc(model.fx,time.t0,time.dt,time.tf,model.x0tr,1e-200);
    ym=zeros(time.nSteps,model.hn);
    for ii=1:1:length(t)
        ym(ii,:)=model.hx(x_mc(ii,:))+(model.sR*randn(model.hn,1));
    end
   filter.x0_filt_start = mvnrnd(model.x0tr', model.P0tr)';
    
    filter.ymeas=ym;
    filter.truth=x_mc;
    
    model.name=strcat('CKF_allfilters_data_set_no_',num2str(i),'_',num2str(j),'_',num2str(k));
    
    data_set{i,j,k,1}=model;
    data_set{i,j,k,2}=filter;
    data_set{i,j,k,3}=pf;
    data_set{i,j,k,4}=time;
%     nnn=nnn+1;
     end
    end
end
% %% just for propagation the montecarlo runs
% X_mc=zeros(time.nSteps,model.fn,2000);
% x0trr = mvnrnd(model.x0tr', model.P0tr,2000);
% parfor i=1:2000
% i
%     [t,x_mc]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0trr(i,:)',model.Qt_sq);
%     X_mc(:,:,i)=x_mc;
%        
% end
%  mean_of_mc_runs=zeros(time.nSteps,model.fn);
%   for r=1:1:time.nSteps
%         for c=1:1:model.fn
%             mean_of_mc_runs(r,c)=mean(X_mc(r,c,:));
%         end
%   end
% %   plot(mean_of_mc_runs(:,1),mean_of_mc_runs(:,2))
%   save X_mc
%   clear X_mc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run in parallel

for i=1:1:10
    for j=1:1:10
        [i,j]
        parfor k=1:NNN
k
%  tic
results_filter=Allfilters(data_set{i,j,k,1},data_set{i,j,k,2},data_set{i,j,k,3},data_set{i,j,k,4});
% toc
  xNNN_mc(:,:,k)=results_filter.mc;
  YNNN_mc(:,:,k)=results_filter.ymeas;
  xNNN_mopf(:,:,k)=results_filter.mopf;
  xNNN_mupf(:,:,k)=results_filter.mupf;
  PNNN_pf(:,:,k)=results_filter.Ppf;
  xNNN_gh(:,:,k)=results_filter.mu_gh;
  xNNN_ekf(:,:,k)=results_filter.mu_ekf;
  xNNN_cut8(:,:,k)=results_filter.mu_cut8;
  xNNN_cut6(:,:,k)=results_filter.mu_cut6;
  xNNN_cut4(:,:,k)=results_filter.mu_cut4;
  xNNN_ckf(:,:,k)=results_filter.mu_ckf;
  xNNN_ukf(:,:,k)=results_filter.mu_ukf;
  PNNN_gh(:,:,k)=results_filter.P_gh;
  PNNN_ekf(:,:,k)=results_filter.P_ekf;
  PNNN_cut8(:,:,k)=results_filter.P_cut8;
  PNNN_cut6(:,:,k)=results_filter.P_cut6;
  PNNN_cut4(:,:,k)=results_filter.P_cut4;
  PNNN_ckf(:,:,k)=results_filter.P_ckf;
  PNNN_ukf(:,:,k)=results_filter.P_ukf; 
  PF_complete_data{k}=results_filter.PFdata;

        end
  save(strcat('FEB19_CKF_allfilters_data_set_no_',num2str(i),'_',num2str(j),'_200'),'xNNN_mc','YNNN_mc','xNNN_ukf','xNNN_ckf','xNNN_cut4','xNNN_cut6','xNNN_cut8','xNNN_gh','xNNN_ekf','xNNN_mupf','xNNN_mopf','PNNN_ukf','PNNN_ckf','PNNN_cut4','PNNN_cut6','PNNN_cut8','PNNN_gh','PNNN_ekf','PNNN_pf')

    end
end
matlabpool close
% %% % Averaging over all the runs
% est_fin_ukf=zeros(time.nSteps,3);
% est_fin_ckf=zeros(time.nSteps,3);
% est_fin_cut4=zeros(time.nSteps,3);
% est_fin_cut6=zeros(time.nSteps,3);
% est_fin_cut8=zeros(time.nSteps,3);
% est_fin_gh=zeros(time.nSteps,3);
% est_fin_ekf=zeros(time.nSteps,3);
% est_fin_mupf=zeros(time.nSteps,3);
% est_fin_mopf=zeros(time.nSteps,3);
% 
% avg_traj_ukf=zeros(time.nSteps,model.fn);
% avg_traj_ckf=zeros(time.nSteps,model.fn);
% avg_traj_cut4=zeros(time.nSteps,model.fn);
% avg_traj_cut6=zeros(time.nSteps,model.fn);
% avg_traj_cut8=zeros(time.nSteps,model.fn);
% avg_traj_gh=zeros(time.nSteps,model.fn);
% avg_traj_ekf=zeros(time.nSteps,model.fn);
% avg_traj_mupf=zeros(time.nSteps,model.fn);
% avg_traj_mopf=zeros(time.nSteps,model.fn);
% avg_traj_mc=zeros(time.nSteps,model.fn);
% avg_traj_meas=zeros(time.nSteps,model.hn);
% 
% Pest_fin_ukf=zeros(time.nSteps,model.hn^2);
% Pest_fin_ckf=zeros(time.nSteps,model.hn^2);
% Pest_fin_cut4=zeros(time.nSteps,model.hn^2);
% Pest_fin_cut6=zeros(time.nSteps,model.hn^2);
% Pest_fin_cut8=zeros(time.nSteps,model.hn^2);
% Pest_fin_gh=zeros(time.nSteps,model.hn^2);
% Pest_fin_pf=zeros(time.nSteps,model.hn^2);
% 
% save(model.name,'model','time','pf','filter','PF_complete_data','avg_traj_ukf','avg_traj_ckf','avg_traj_cut4','avg_traj_cut6','avg_traj_cut8','avg_traj_gh','avg_traj_ekf','avg_traj_mopf','avg_traj_mupf','avg_traj_mc','est_fin_ukf','est_fin_ckf','est_fin_cut4','est_fin_cut6','est_fin_cut8','est_fin_gh','est_fin_ekf','est_fin_mupf','est_fin_mopf','xNNN_mc','YNNN_mc','xNNN_ukf','xNNN_ckf','xNNN_cut4','xNNN_cut6','xNNN_cut8','xNNN_gh','xNNN_ekf','xNNN_mupf','xNNN_mopf','PNNN_ukf','PNNN_ckf','PNNN_cut4','PNNN_cut6','PNNN_cut8','PNNN_gh','PNNN_ekf','PNNN_pf')
%   for r=1:1:time.nSteps
%         for c=1:1:model.fn
%            est_fin_ukf(r,c)= sqrt(mean((xNNN_ukf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
%            est_fin_ckf(r,c)= sqrt(mean((xNNN_ckf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
%            est_fin_cut4(r,c)= sqrt(mean((xNNN_cut4(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
%            est_fin_cut6(r,c)= sqrt(mean((xNNN_cut6(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
%            est_fin_cut8(r,c)= sqrt(mean((xNNN_cut8(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
%            est_fin_gh(r,c)= sqrt(mean((xNNN_gh(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
%            est_fin_mupf(r,c)= sqrt(mean((xNNN_mupf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
%            est_fin_mopf(r,c)= sqrt(mean((xNNN_mopf(r,c,:)-xNNN_mc(r,c,:)).^2))^2;
%         end
%         for c=1:1:model.fn^2
%              Pest_fin_ukf(r,c)= sqrt(mean(PNNN_ukf(r,c,:)));
%              Pest_fin_ckf(r,c)= sqrt(mean(PNNN_ckf(r,c,:)));
%              Pest_fin_cut4(r,c)= sqrt(mean(PNNN_cut4(r,c,:)));
%              Pest_fin_cut6(r,c)= sqrt(mean(PNNN_cut6(r,c,:)));
%              Pest_fin_cut8(r,c)= sqrt(mean(PNNN_cut8(r,c,:)));
%              Pest_fin_gh(r,c)= sqrt(mean(PNNN_gh(r,c,:)));
%              Pest_fin_pf(r,c)= sqrt(mean(PNNN_pf(r,c,:)));
%         end
%   end
% 
% for r=1:1:time.nSteps
%         for c=1:1:model.fn
% avg_traj_ukf(r,c)=mean(xNNN_ukf(r,c,:));
% avg_traj_ckf(r,c)=mean(xNNN_ckf(r,c,:));
% avg_traj_cut4(r,c)=mean(xNNN_cut4(r,c,:));
% avg_traj_cut6(r,c)=mean(xNNN_cut6(r,c,:));
% avg_traj_cut8(r,c)=mean(xNNN_cut8(r,c,:));
% avg_traj_gh(r,c)=mean(xNNN_gh(r,c,:));
% avg_traj_ekf(r,c)=mean(xNNN_ekf(r,c,:));
% avg_traj_mupf(r,c)=mean(xNNN_mupf(r,c,:));
% avg_traj_mopf(r,c)=mean(xNNN_mopf(r,c,:));
% avg_traj_mc(r,c)=mean(xNNN_mc(r,c,:));
% 
%         end
% end
% 
% save(model.name,'data_set','model','time','pf','filter','PF_complete_data','avg_traj_ukf','avg_traj_ckf','avg_traj_cut4','avg_traj_cut6','avg_traj_cut8','avg_traj_gh','avg_traj_ekf','avg_traj_mopf','avg_traj_mupf','avg_traj_mc','est_fin_ukf','est_fin_ckf','est_fin_cut4','est_fin_cut6','est_fin_cut8','est_fin_gh','est_fin_ekf','est_fin_mupf','est_fin_mopf','xNNN_mc','YNNN_mc','xNNN_ukf','xNNN_ckf','xNNN_cut4','xNNN_cut6','xNNN_cut8','xNNN_gh','xNNN_ekf','xNNN_mupf','xNNN_mopf','PNNN_ukf','PNNN_ckf','PNNN_cut4','PNNN_cut6','PNNN_cut8','PNNN_gh','PNNN_ekf','PNNN_pf','Pest_fin_ukf','Pest_fin_ckf','Pest_fin_cut4','Pest_fin_cut6','Pest_fin_cut8','Pest_fin_gh','Pest_fin_pf')
% %  matlabpool close
%   
% %% Plotting
%  t=time.tspan;
%   figure(1)
% subplot(2,1,1)
% plot(t,est_fin_mupf(:,1),'-.',t,est_fin_mopf(:,1),':',t,est_fin_ckf(:,1),t,est_fin_ukf(:,1),t,est_fin_cut4(:,1),t,est_fin_cut6(:,1),t,est_fin_cut8(:,1),t,est_fin_gh(:,1))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% subplot(2,1,2)
% plot(t,Pest_fin_pf(:,1),t,Pest_fin_ckf(:,1),t,Pest_fin_ukf(:,1),t,Pest_fin_cut4(:,1),t,Pest_fin_cut6(:,1),t,Pest_fin_cut8(:,1),t,Pest_fin_gh(:,1))
% legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% 
% figure(2)
% subplot(2,1,1)
% plot(t,est_fin_mupf(:,2),'-.',t,est_fin_mopf(:,2),':',t,est_fin_ckf(:,2),t,est_fin_ukf(:,2),t,est_fin_cut4(:,2),t,est_fin_cut6(:,2),t,est_fin_cut8(:,2),t,est_fin_gh(:,2))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% subplot(2,1,2)
% plot(t,Pest_fin_pf(:,7),t,Pest_fin_ckf(:,7),t,Pest_fin_ukf(:,7),t,Pest_fin_cut4(:,7),t,Pest_fin_cut6(:,7),t,Pest_fin_cut8(:,7),t,Pest_fin_gh(:,7))
% legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% 
% figure(3)
% subplot(2,1,1)
% plot(t,est_fin_mupf(:,3),'-.',t,est_fin_mopf(:,3),':',t,est_fin_ckf(:,3),t,est_fin_ukf(:,3),t,est_fin_cut4(:,3),t,est_fin_cut6(:,3),t,est_fin_cut8(:,3),t,est_fin_gh(:,3))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% subplot(2,1,2)
% plot(t,Pest_fin_pf(:,13),t,Pest_fin_ckf(:,13),t,Pest_fin_ukf(:,13),t,Pest_fin_cut4(:,13),t,Pest_fin_cut6(:,13),t,Pest_fin_cut8(:,13),t,Pest_fin_gh(:,13))
% legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% 
% figure(4)
% subplot(2,1,1)
% plot(t,est_fin_mupf(:,4),'-.',t,est_fin_mopf(:,4),':',t,est_fin_ckf(:,4),t,est_fin_ukf(:,4),t,est_fin_cut4(:,4),t,est_fin_cut6(:,4),t,est_fin_cut8(:,4),t,est_fin_gh(:,4))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% subplot(2,1,2)
% plot(t,Pest_fin_pf(:,19),t,Pest_fin_ckf(:,19),t,Pest_fin_ukf(:,19),t,Pest_fin_cut4(:,19),t,Pest_fin_cut6(:,19),t,Pest_fin_cut8(:,19),t,Pest_fin_gh(:,19))
% legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% 
% figure(5)
% subplot(2,1,1)
% plot(t,est_fin_mupf(:,5),'-.',t,est_fin_mopf(:,5),':',t,est_fin_ckf(:,5),t,est_fin_ukf(:,5),t,est_fin_cut4(:,5),t,est_fin_cut6(:,5),t,est_fin_cut8(:,5),t,est_fin_gh(:,5))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% subplot(2,1,2)
% plot(t,Pest_fin_pf(:,25),t,Pest_fin_ckf(:,25),t,Pest_fin_ukf(:,25),t,Pest_fin_cut4(:,25),t,Pest_fin_cut6(:,25),t,Pest_fin_cut8(:,25),t,Pest_fin_gh(:,25))
% legend('PF','ckf','ukf','cut4','cut6','cut8','gh')
% 
% figure(6)
% plot(avg_traj_mupf(:,1),avg_traj_mupf(:,2),'-.',avg_traj_mopf(:,1),avg_traj_mopf(:,2),':',avg_traj_ckf(:,1),avg_traj_ckf(:,2),avg_traj_ukf(:,1),avg_traj_ukf(:,2),avg_traj_cut4(:,1),avg_traj_cut4(:,2),avg_traj_cut6(:,1),avg_traj_cut6(:,2),avg_traj_cut8(:,1),avg_traj_cut8(:,2),avg_traj_gh(:,1),avg_traj_gh(:,2),avg_traj_mc(:,1),avg_traj_mc(:,2),'+')
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh','MC')
% figure(7)
% plot(t,avg_traj_mupf(:,1),'-.',t,avg_traj_mopf(:,1),':',t,avg_traj_ckf(:,1),t,avg_traj_ukf(:,1),t,avg_traj_cut4(:,1),t,avg_traj_cut6(:,1),t,avg_traj_cut8(:,1),t,avg_traj_gh(:,1))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')