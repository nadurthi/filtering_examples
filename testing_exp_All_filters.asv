close all;
clear all;
clc;
global T
tic

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
model.hn = 2;               % measurement dimensionality
model.hx = @CKF_eg_meas_disc;
model.Q = Q;
model.sQ = sqrtm(model.Q);
model.R = R;
model.sR=sqrtm(R);
model.t0=time.t0;
model.dt=time.dt;
model.tf=time.tf;
model.Qtruth=1e-200;

% n=model.fn;
%% ------------------------------------------------------------------------
% particle filter settings
%--------------------------------------------------------------------------
pf.no_particles = 5000;
pf.no_bins = 1000;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;                  % treshold for the number of effective samples

%% --------------------------------------------------------------
%filter props
filter_paras.ukf_kappa=1;
filter_paras.gh_pts=7;
filter.freq=5; %% This is actually the number of 'dt' steps after which a meas updt is done
%---------------------------------
%% __-------------------------------------------------------
%TRUTH INITIALIZATION of the model
  x0tr=[1000,300,1000,0,omg]';
  P0tr=diag([100,10,100,10,100e-6]);
  
 
%% %%%%%%%%%%%%START FILTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NNN=1; %% Simulate over NNN iterations


xNNN_mc=zeros(time.nSteps,model.fn,NNN);
YNNN_mc=zeros(time.nSteps,model.hn,NNN);
xNNN_ukf=zeros(time.nSteps,model.fn,NNN);
xNNN_ckf=zeros(time.nSteps,model.fn,NNN);
xNNN_cut4=zeros(time.nSteps,model.fn,NNN);
xNNN_cut6=zeros(time.nSteps,model.fn,NNN);
xNNN_cut8=zeros(time.nSteps,model.fn,NNN);
xNNN_gh=zeros(time.nSteps,model.fn,NNN);
xNNN_mupf=zeros(time.nSteps,model.fn,NNN);
xNNN_mopf=zeros(time.nSteps,model.fn,NNN);
PNNN_ukf=zeros(time.nSteps,model.fn^2,NNN);
PNNN_ckf=zeros(time.nSteps,model.fn^2,NNN);
PNNN_cut4=zeros(time.nSteps,model.fn^2,NNN);
PNNN_cut6=zeros(time.nSteps,model.fn^2,NNN);
PNNN_cut8=zeros(time.nSteps,model.fn^2,NNN);
PNNN_gh=zeros(time.nSteps,model.fn^2,NNN);
PNNN_pf=zeros(time.nSteps,model.fn^2,NNN);

histdata=cell(time.nSteps,NNN);

for j=1:1:NNN
    
    %%-------------------------------------------------------------
    %%monte carlo truth generation using dynamics
    %%_-----------------------------------------------------------
    [t,x_mc]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0tr,model.Qtruth);
    %generating the measurement
    ym=zeros(size(x_mc,1),model.hn);
    for i=1:1:time.nSteps
        ym(i,:)=(model.hx(x_mc(i,:)')+model.sR*randn(model.hn,1))';
    end
    xNNN_mc(:,:,j)=x_mc;
    YNNN_mc(:,:,j)=ym;
%% ------------------------------------------------------------------------
%FILTERS COnditions

    %%Randomly pick initial point as filter initial cond
    x_0f=mvnrnd(x0tr,P0tr)';
    P_0f=P0tr;
    
    %%Initiate all the filters to the same starting point
    
    % UKF filter
    mu_ut=x_0f;
    P_ut=P_0f;
    xNNN_ukf(1,:,j)=mu_ut';
    PNNN_ukf(1,:,j)=reshape(P_ut,1,model.fn^2);
    
    % CKF filter
    mu_ckf=x_0f;
    P_ckf=P_0f;
    xNNN_ckf(1,:,j)=mu_ckf';
    PNNN_ckf(1,:,j)=reshape(P_ckf,1,model.fn^2);
 
    % CUT4 filter
    mu_cut4=x_0f;
    P_cut4=P_0f;
    xNNN_cut4(1,:,j)=mu_cut4';
    PNNN_cut4(1,:,j)=reshape(P_cut4,1,model.fn^2);
    
    % CUT6 filter
    mu_cut6=x_0f;
    P_cut6=P_0f;
    xNNN_cut6(1,:,j)=mu_cut6';
    PNNN_cut6(1,:,j)=reshape(P_cut6,1,model.fn^2);
    
    % CUT8 filter
    mu_cut8=x_0f;
    P_cut8=P_0f;
    xNNN_cut8(1,:,j)=mu_cut8';
    PNNN_cut8(1,:,j)=reshape(P_cut8,1,model.fn^2);
    
    % ghKF filter
    mu_gh=x_0f;
    P_gh=P_0f;
    xNNN_gh(1,:,j)=mu_gh';
    PNNN_gh(1,:,j)=reshape(P_gh,1,model.fn^2);
       
       
       % PF
       mu_pf=x_0f;
       mo_pf=x_0f;
       P_pf=P_0f;
       PNNN_pf(1,:,j)=reshape(P_pf,1,model.fn^2);
       xNNN_mupf(1,:,j)=mu_pf';
       xNNN_mopf(1,:,j)=mo_pf';
       X_pf=repmat(mu_pf,1,pf.no_particles)+sqrtm(P_pf)*randn(model.fn,pf.no_particles);
       w_pf = ones(1, pf.no_particles) / pf.no_particles;
       [PF.x{1},PF.P{1},PF.minB{1},PF.maxB{1},tmp_X_pf] = getPFdata(X_pf,w_pf);
       histdata{1,j}=tmp_X_pf;
       
        
  for k=2:1:time.nSteps
disp([num2str(k),' of ',num2str(time.nSteps)] )
        
        if rem(k,filter.freq)==0
            zm=ym(k,:)';
        else
            zm=-1234;
        end
      [mu_ut,P_ut]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_ut,P_ut,zm,'ut',filter_paras.ukf_kappa);
      [mu_ckf,P_ckf]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_ckf,P_ckf,zm,'ckf',0);
      [mu_cut4,P_cut4]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_cut4,P_cut4,zm,'cut4',0);
      [mu_cut6,P_cut6]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_cut6,P_cut6,zm,'cut6',0);
      [mu_cut8,P_cut8]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_cut8,P_cut8,zm,'cut8',0);
      [mu_gh,P_gh]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_gh,P_gh,zm,'gh',filter_paras.gh_pts); 
      [X_pf,w_pf]=BOOTSTRAP_Pfilter_disc_UPDT_disc_MEAS(model,pf,X_pf,w_pf,zm); 
      
      
      %storing data for this run
    xNNN_ukf(k,:,j)=mu_ut';
    PNNN_ukf(k,:,j)=reshape(P_ut,1,model.fn^2);
      
    xNNN_ckf(k,:,j)=mu_ckf';
    PNNN_ckf(k,:,j)=reshape(P_ckf,1,model.fn^2);
   
    xNNN_cut4(k,:,j)=mu_cut4';
    PNNN_cut4(k,:,j)=reshape(P_cut4,1,model.fn^2);
    
    xNNN_cut6(k,:,j)=mu_cut6';
    PNNN_cut6(k,:,j)=reshape(P_cut6,1,model.fn^2);
    
    xNNN_cut8(k,:,j)=mu_cut8';
    PNNN_cut8(k,:,j)=reshape(P_cut8,1,model.fn^2);

    xNNN_gh(k,:,j)=mu_gh';
    PNNN_gh(k,:,j)=reshape(P_gh,1,model.fn^2);
      
    [PFx,P_pf,PFminB,PFmaxB,tmp_X_pf] = getPFdata(X_pf, w_pf);
    histdata{k,j}=tmp_X_pf;
    
    xNNN_mupf(k,:,j)=PFx';
    PNNN_pf(k,:,j)=reshape(P_pf,1,model.fn^2);
    
    MODE=zeros(model.fn,1);
    for mo=1:1:model.fn
    rang=PFminB(mo):(PFmaxB(mo)-PFminB(mo))/(pf.no_bins-1):PFmaxB(mo);
    nb=histc(tmp_X_pf(mo,:),rang);
    [val,ind]=max(nb);
    MODE(mo,1)=(rang(ind)+rang(ind+1))/2;
    end
    
    xNNN_mopf(k,:,j)=MODE';
    
  end

end
% Averaging over all the runs
est_fin_ukf=zeros(time.nSteps,3);
est_fin_ckf=zeros(time.nSteps,3);
est_fin_cut4=zeros(time.nSteps,3);
est_fin_cut6=zeros(time.nSteps,3);
est_fin_cut8=zeros(time.nSteps,3);
est_fin_gh=zeros(time.nSteps,3);
est_fin_mupf=zeros(time.nSteps,3);
est_fin_mopf=zeros(time.nSteps,3);

for r=1:1:length(t)
        for c=0:1:2
%            x_fin_mc(r,c)= mean(x100_mc(r,c,:));
%            x_fin_ekf(r,c)=mean(x100_ekf(r,c,:));
           if c<=1
           est_fin_ukf(r,c+1)= sqrt(mean((xNNN_ukf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_ukf(r,c+3,:)-xNNN_mc(r,c+3,:)).^2));
           est_fin_ckf(r,c+1) =sqrt(mean((xNNN_ckf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_ckf(r,c+3,:)-xNNN_mc(r,c+3,:)).^2));
           est_fin_cut4(r,c+1)=sqrt(mean((xNNN_cut4(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_cut4(r,c+3,:)-xNNN_mc(r,c+3,:)).^2));
           est_fin_cut6(r,c+1)=sqrt(mean((xNNN_cut6(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_cut6(r,c+3,:)-xNNN_mc(r,c+3,:)).^2));
           est_fin_cut8(r,c+1)=sqrt(mean((xNNN_cut8(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_cut8(r,c+3,:)-xNNN_mc(r,c+3,:)).^2));
           est_fin_gh(r,c+1)=  sqrt(mean((xNNN_gh(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_gh(r,c+3,:)-xNNN_mc(r,c+3,:)).^2));
           est_fin_mupf(r,c+1)=  sqrt(mean((xNNN_mupf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_mupf(r,c+3,:)-xNNN_mc(r,c+3,:)).^2));
           est_fin_mopf(r,c+1)=  sqrt(mean((xNNN_mopf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_mopf(r,c+3,:)-xNNN_mc(r,c+3,:)).^2));
           else
           est_fin_ukf(r,c+1)= sqrt(mean((xNNN_ukf(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_ckf(r,c+1) =sqrt(mean((xNNN_ckf(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_cut4(r,c+1)=sqrt(mean((xNNN_cut4(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_cut6(r,c+1)=sqrt(mean((xNNN_cut6(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_cut8(r,c+1)=sqrt(mean((xNNN_cut8(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_gh(r,c+1)=  sqrt(mean((xNNN_gh(r,5,:)-xNNN_mc(r,5,:)).^2)); 
           est_fin_mupf(r,c+1)=  sqrt(mean((xNNN_mupf(r,5,:)-xNNN_mc(r,5,:)).^2)); 
           est_fin_mopf(r,c+1)=  sqrt(mean((xNNN_mopf(r,5,:)-xNNN_mc(r,5,:)).^2)); 
           end
        end
end
save(strcat('ckf_eg_PF_also_run_nos_',num2str(NNN)),'model','pf','filter_paras','histdata','est_fin_ukf','est_fin_ckf','est_fin_cut4','est_fin_cut6','est_fin_cut8','est_fin_gh','est_fin_mupf','est_fin_mopf','xNNN_mc','YNNN_mc','xNNN_ukf','xNNN_ckf','xNNN_cut4','xNNN_cut6','xNNN_cut8','xNNN_gh','xNNN_mupf','xNNN_mopf','PNNN_ukf','PNNN_ckf','PNNN_cut4','PNNN_cut6','PNNN_cut8','PNNN_gh','PNNN_pf')





%% Plotting
t=time.tspan;
figure(1)
plot(t,est_fin_mupf(:,1),'-.',t,est_fin_mopf(:,1),':',t,est_fin_ckf(:,1),t,est_fin_ukf(:,1),t,est_fin_cut4(:,1),t,est_fin_cut6(:,1),t,est_fin_cut8(:,1),t,est_fin_gh(:,1))
legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
figure(2)
plot(t,est_fin_mupf(:,2),'-.',t,est_fin_mopf(:,2),':',t,est_fin_ckf(:,2),t,est_fin_ukf(:,2),t,est_fin_cut4(:,2),t,est_fin_cut6(:,2),t,est_fin_cut8(:,2),t,est_fin_gh(:,2))
legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
figure(3)
plot(t,(180/pi)*est_fin_mupf(:,3),'-.',t,(180/pi)*est_fin_mopf(:,3),':',t,(180/pi)*est_fin_ckf(:,3),t,(180/pi)*est_fin_ukf(:,3),t,(180/pi)*est_fin_cut4(:,3),t,(180/pi)*est_fin_cut6(:,3),t,(180/pi)*est_fin_cut8(:,3),t,(180/pi)*est_fin_gh(:,3))
legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
figure(4)
plot(xNNN_mupf(:,1,end),xNNN_mupf(:,3,end),'-.',xNNN_mopf(:,1,end),xNNN_mopf(:,3,end),':',xNNN_ckf(:,1,end),xNNN_ckf(:,3,end),xNNN_ukf(:,1,end),xNNN_ukf(:,3,end),xNNN_cut4(:,1,end),xNNN_cut4(:,3,end),xNNN_cut6(:,1,end),xNNN_cut6(:,3,end),xNNN_cut8(:,1,end),xNNN_cut8(:,3,end),xNNN_gh(:,1,end),xNNN_gh(:,3,end),x_mc(:,1),x_mc(:,3),'--',ym(:,1).*cos(ym(:,2)),ym(:,1).*sin(ym(:,2)),'k*','LineWidth',2)
legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh','MC','sensor')



