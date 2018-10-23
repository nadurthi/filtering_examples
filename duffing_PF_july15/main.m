%--------------------------------------------------------------------------
%
% Gabriel Terejanu (terejanu@buffalo.edu)
%--------------------------------------------------------------------------
close all;
clear all;
clc;
global T
tic
%% ------------------------------------------------------------------------
% time
%--------------------------------------------------------------------------
time.t0 = 0;
time.dt = 5;
time.tf = 125+90+125+30+125;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);

%% ------------------------------------------------------------------------
% model

    T=time.dt;
%     omg=-3*pi/180;
%     q1=0.1;
%     q2=1.75*10^(-4);
%     sigr=10;
%     sigth=0.1833*(pi/180)*1;
% %    sigr=10;
% % sigth=0.3*pi/180;
%     M=[T^3/3,T^2/2;T^2/2,T];
%     Q=blkdiag(q1*M,q1*M,q2*T);
%     R=diag([sigr^2,sigth^2]);
L1=0.16;
L2=0.01;
Q=L1*[T^3/3,0,T^2/2,0,0;
        0,T^3/3,0,T^2/2,0;
        T^2/2,0,T,0,0;
        0,T^2/2,0,T,0;
        0,0,0,0,T*L2/L1];
    R=diag([(1e3)^2,(1e3)^2]);
%--------------------------------------------------------------------------
% sigma=0.05;
model.fn = 5;               % state space dimensionality
model.fx = 'KIRB_CT_eg_dyn_disc';
% model.mx = 'CKF_eg_meas_disc';
model.hn = 2;               % measurement dimensionality
model.hx = 'KIRB_eg_meas_disc';
model.Q = Q;
model.sQ = sqrt(model.Q);
model.R = R;

%% ------------------------------------------------------------------------
% prior Gaussian Sum   !!!!!! make them cells (mu & sig) !!!!!!!!
%--------------------------------------------------------------------------
%prior.n = 1;
prior.mu{1} = [25000,10000,-120,0,0.000001];
%prior.sig{1} = 0.01*eye(4);
% prior.mu{2} = 1;
 prior.sig{1} = diag([50^2,50^2,100,100,1]);
% prior.weig = [.1 .9];
%prior.weig=1;

% prior dist. for params.

% prior.p = [0.9 0.5 -1.45 0.5];
% 
% prior.p_t = [1.3572  0 -1.0658  0];

%% ------------------------------------------------------------------------
% measurements available every n time steps
%--------------------------------------------------------------------------
measTs = 1;
measAvail = zeros(1,time.nSteps);
measAvail(1:measTs:time.nSteps) = 1;

%% ------------------------------------------------------------------------
% particle filter settings
%--------------------------------------------------------------------------
pf.no_particles = 5000;
pf.no_bins = 100;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;                  % treshold for the number of effective samples

%% ------------------------------------------------------------------------
% create truth & measurement
%--------------------------------------------------------------------------
fprintf('  - create truth & measurements\n');
xt = cell(1, time.nSteps);
y_meas = cell(1, time.nSteps);

% select a gaussian component
% gs_sel = 1;
% u = rand;
% u_total = 0;
% for j = 1 : prior.n
%     if ((u >= u_total) && (u < u_total + prior.weig(j)))
%         gs_sel = j;
%         break;
%     end
%     u_total = u_total + prior.weig(j);
% end

% draw a sample from the chosen gaussian component
% xt{1} = prior.mu{gs_sel}' + chol(prior.sig{gs_sel})' * randn(model.fn,1); 
xt{1} = prior.mu{1}';

% for i = 1 : model.fn-model.hn
%     xt_param = prior.p_t((i-1)*2+1) + prior.p_t(i*2)*rand(1,1);
%     xt{1} = [xt{1}; xt_param];
% end;

% y(1,:) = feval(model.fx, xt{1});
% for k=2:time.nSteps
%    y(k,:)=feval(model.fx,y(k-1,:)') ;
% end
    [t,x_mc1]=ode45_discc(@KIRB_UM_eg_dyn_disc,time.t0 ,time.dt,125,[25000,10000,-120,0]',1e-200);
    [t,x_mc2]=ode45_discc(@KIRB_CT_eg_dyn_disc,t(end),time.dt,t(end)+90,[x_mc1(end,:)';1*pi/180],1e-200);
    [t,x_mc3]=ode45_discc(@KIRB_UM_eg_dyn_disc,t(end),time.dt,t(end)+125,x_mc2(end,1:end-1)',1e-200);
    [t,x_mc4]=ode45_discc(@KIRB_CT_eg_dyn_disc,t(end),time.dt,t(end)+30,[x_mc3(end,:)';-3*pi/180],1e-200);
    [t,x_mc5]=ode45_discc(@KIRB_UM_eg_dyn_disc,t(end),time.dt,t(end)+125,x_mc4(end,1:end-1)',1e-200);

    y=[x_mc1(1:end-1,1:4);x_mc2(1:end-1,1:4);x_mc3(1:end-1,1:4);x_mc4(1:end-1,1:4);x_mc5(1:end,1:4)];
% y = feval(model.hx, xt{1},time.tspan);
% get the trajectory of the sample over time
for k = 2 : time.nSteps
%     Y1(k)=y(1);
%     Y2(k)=y(2);
    xt{k} = y(k,:)'; %+ model.sQ * randn(model.fn,1);
    y_meas{k} = feval(model.hx, xt{k}) + sqrtm(R) * randn(model.hn,1);
%     xt{k}
%     y_meas{k}
%     keyboard
end
ymeas=[y_meas{1:1:end}]';
plot(y(:,1),y(:,2),ymeas(:,1),ymeas(:,2),'k*')
%% ------------------------------------------------------------------------
% PF - for particle filters
%--------------------------------------------------------------------------
fprintf('  - create initial particles\n');

% % generate i.c. samples
% u_pf = rand(1,pf.no_particles);
% tu_pf = 0;
% X_pf = [];
% for i = 1 : prior.n
%     % for each gaussian component draw samples dictated by its weight
%     % magnitude
%     ind = find((u_pf >= tu_pf) & (u_pf < tu_pf + prior.weig(i)));
%     Xtmp_pf = repmat(prior.mu{i}',1,length(ind)) + chol(prior.sig{i})' * randn(model.fn,length(ind));
%     
%     % collect all the samples and advance
%     X_pf = [X_pf Xtmp_pf];
%     tu_pf = tu_pf + prior.weig(i);
% end
% X_pf_state = repmat(prior.mu{1}',1,pf.no_particles);
% X_pf = X_pf_state;
X_pf=repmat(prior.mu{1}',1,pf.no_particles)+sqrtm(prior.sig{1})*randn(model.fn,pf.no_particles);
% for i = 1 : model.fn-model.hn
%     X_pf_param = Pm(:,i)';
% %     X_pf_param = prior.p((i-1)*2+1) + prior.p(i*2)*rand(1,pf.no_particles);
%     X_pf = [X_pf; X_pf_param];
% end;

w_pf = ones(1, pf.no_particles) / pf.no_particles;

%% ------------------------------------------------------------------------
% compute initial estimates
%--------------------------------------------------------------------------

[PF.x{1},PF.P{1},PF.minB{1},PF.maxB{1},tmp_X_pf] = getPFdata(X_pf,w_pf);

%% ------------------------------------------------------------------------
% RUN BOOSTRAP PARTICLE FILTER
%--------------------------------------------------------------------------   
% M_x1=[];M_x2=[];M_p1=[];M_p2=[];
% Nm=10;
for k = measTs+1 : measTs:time.nSteps
    fprintf('  - run filters tStep %d / %d\n', k, time.nSteps);    


%% ------------------------------------------------------------------------
% PF - propagate particles
%--------------------------------------------------------------------------
    fprintf('     - Particle Filter\n'); 

    [X_pf] = PF_time_update1(X_pf, model, time, k, measTs);
   
    if (measAvail(k) == 1)

        [X_pf, w_pf] = PF_meas_update(X_pf, w_pf, model, pf, y_meas{k});
    end

    [PF.x{k},PF.P{k},PF.minB{k},PF.maxB{k},tmp_X_pf] = getPFdata(X_pf, w_pf);
    histdata{k}=tmp_X_pf;

end;
mean_pf=zeros(time.nSteps,model.fn);
mean_pf(1,:)=xt{1};
for k=	2:1:time.nSteps  
    mean_pf(k,:)=mean(histdata{k},2)';
end
ymeas=[y_meas{1:1:end}]';
% plot(mean_pf(:,1),mean_pf(:,2),ymeas(:,1),ymeas(:,2),'k*')
plot(y(:,1),y(:,2),'r',mean_pf(:,1),mean_pf(:,2),ymeas(:,1),ymeas(:,2),'k*')
% p=1;
% s=2;
% for k=2:1:time.nSteps
% if p==26
%     s=s+1;
%     p=1;
% end
%     figure(s)
% for i=1:1:5
% subplot(5,5,p)
% p=p+1;
% a=histdata{k};
% hist(a(i,:),500)
% title(['time: ',num2str(k),' state: ',num2str(i)])
% end
% end
