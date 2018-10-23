% function lorenz_filtering_eg()
clear
clc
%close all
% global T
%distcomp.feature( 'LocalUseMpiexec', false );
%   matlabpool open 3
% warning off;
%% ------------------------------------------------------------------------
% time
%--------------------------------------------------------------------------
time.t0 = 0;
time.dt = 0.5;
time.tf = 30;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);
%__________________________________________________________________________

%% ------------------------------------------------------------------------
% model
T=time.dt;

Qf=0.00000001*diag([1,1,1,1,1]);

Qt=0*diag([0,0,2.4064*1e-5,2.4064*1e-5,0]);

% sigr=100*1e-3;
% sigth=10*1e-3;

% R=diag([sigr^2,sigth^2]);
% R=diag([sigth^2]);



model.fn = 5;               % state space dimensionality
model.fx = @lorenz_cont;
model.fx_jac=@ABC;
model.dynamics='continuous';

R=diag([2]);
model.hn =1;               % measurement dimensionality
model.hx =@(x)[x(1)];
model.hx_jac=@ABC;

model.Q =Qf;
model.sQ =sqrtm(Qf);
model.R = R;
model.sR=sqrtm(R);
model.Qtruth=Qt;
model.Qt_sq=sqrtm(Qt);
model.dt=time.dt;

model.para_dt=0;

% model.x0tr=[6500.4,349.14,-1.8093,-6.7967,0.6932]';
% model.P0tr=diag([1e-5,1e-5,1e-5,1e-5,0]);
% [1.50887,-1.531271,25.46091,12,25]';
model.x0tr=[1.50887,-1.531271,25.46091,10,28]';
model.P0tr=diag([4,4,4,2,4]);
%% ------------------------------------------------------------------------
% particle filter settings
%--------------------------------------------------------------------------
pf.no_particles =5000;
pf.no_bins = 100;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;                  % treshold for the number of effective samples

pf.plotHists='Yes';
pf.plotHiststagname='1meas_verylownoise_';

%% --------------------------------------------------------------
%filter props
x0f=model.x0tr;
P0f=model.P0tr;

filter.paras_ukf_kappa=1;
filter.paras_gh_pts=2;

filter.freq=1; %% This is actually the number of 'dt' steps after which a meas updt is done

filter.x0_filt_start=x0f;
filter.P0_filt_start=P0f;

filter.save_pf_data='false';

%switch on filters
filter.EKF='false';
filter.CKF='false';
filter.UKF='false';
filter.KF='false';
filter.GHKF='false';
filter.PF='true';
filter.CUT4KF='false';
filter.CUT6KF='false';
filter.CUT8KF='false';
filter.truth=[];
filter.ymeas=[];
%% Allocate space for NNN no. of iterations
NNN=1;

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

data_set=cell([NNN,4]);


for k=1:1:NNN
    
    x0trr = model.x0tr;%mvnrnd(model.x0tr', model.P0tr);
    [t,x_mc]=ode45(model.fx,time.tspan,x0trr',model.Qt_sq);
    ym=zeros(time.nSteps,model.hn);
    for ii=1:1:time.nSteps
        ym(ii,:)=(model.hx(x_mc(ii,:)')+(sqrtm(R)*randn(model.hn,1)))';
    end
    filter.x0_filt_start = mvnrnd(model.x0tr', model.P0tr)';
    %     filter.x0_filt_start =[1.50887,-1.531271,25.46091,12,25]';
    filter.ymeas=ym;
    filter.truth=x_mc;
    
    model.name=strcat('lorenz_22may2013_data');
    
    data_set{k,1}=model;
    data_set{k,2}=filter;
    data_set{k,3}=pf;
    data_set{k,4}=time;
    
end

%% just for propagation the montecarlo runs
%X_mc=zeros(time.nSteps,model.fn,2000);
%x0trr = mvnrnd(model.x0tr', model.P0tr,2000);
%%for i=1:2000
%i
%    [t,x_mc]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0trr(i,:)',model.Qt_sq);
%    X_mc(:,:,i)=x_mc;
%
%end
% mean_of_mc_runs=zeros(time.nSteps,model.fn);
%  for r=1:1:time.nSteps
%        for c=1:1:model.fn
%            mean_of_mc_runs(r,c)=mean(X_mc(r,c,:));
%        end
%  end
%  plot(mean_of_mc_runs(:,1),mean_of_mc_runs(:,2))
%  save X_mc
%   clear X_mc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run in parallel
for k=1:NNN
    k
    results_filter=Allfilters(data_set{k,1},data_set{k,2},data_set{k,3},data_set{k,4});
    
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
%       save('PF_4000_lorenz_','model','time','pf','filter','PF_complete_data','xNNN_mc','YNNN_mc','xNNN_ukf','xNNN_ckf','xNNN_cut4','xNNN_cut6','xNNN_cut8','xNNN_gh','xNNN_ekf','xNNN_mupf','xNNN_mopf','PNNN_ukf','PNNN_ckf','PNNN_cut4','PNNN_cut6','PNNN_cut8','PNNN_gh','PNNN_ekf','PNNN_pf')
% s=datestr(now);
% s(find(s==' '))='_';
% s(find(s==':'))='-';
%      save(strcat('PF_4000_lorenz_',s,'_data'),'model','time','pf','filter','PF_complete_data','xNNN_mc','YNNN_mc','xNNN_ukf','xNNN_ckf','xNNN_cut4','xNNN_cut6','xNNN_cut8','xNNN_gh','xNNN_ekf','xNNN_mupf','xNNN_mopf','PNNN_ukf','PNNN_ckf','PNNN_cut4','PNNN_cut6','PNNN_cut8','PNNN_gh','PNNN_ekf','PNNN_pf')


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


% % Plotting
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
% % matlabpool close


[sqrt(sum(sqrt(sum(est_fin_mupf.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_mopf.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_ckf.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_ukf.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_cut4.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_cut6.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_cut8.^2/model.fn,2)).^2)/time.nSteps),sqrt(sum(sqrt(sum(est_fin_gh.^2/model.fn,2)).^2)/time.nSteps)]
[max(max(est_fin_mupf)),max(max(est_fin_mopf)),max(max(est_fin_ckf)),max(max(est_fin_ukf)),max(max(est_fin_cut4)),max(max(est_fin_cut6)),max(max(est_fin_cut8)),max(max(est_fin_gh))]

for i=1:1:size(PF_complete_data{1},1)
    i
    w=PF_complete_data{1}{i,2};
    X=PF_complete_data{1}{i,1};
    [length(unique(round(X(1,:)*10000000)/10000000)),length(unique(round(X(2,:)*10000000)/10000000)),length(unique(round(X(3,:)*10000000)/10000000)),length(unique(round(w*1000000000)/1000000000))]
    sort(w)
    %  hist(w,linspace(0,0.05,1000));
    plot(X(1,:),X(2,:),'ro')
    pause
end
% end