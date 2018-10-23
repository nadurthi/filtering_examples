%--------------------------------------------------------------------------
% 1D - sin(x) Discrete from Continuous example
%
% Gabriel Terejanu (terejanu@buffalo.edu)
%--------------------------------------------------------------------------
close all;
clear all;
clc;
load duff_in_july13
load Pm_5e4
tic
%% ------------------------------------------------------------------------
% time
%--------------------------------------------------------------------------
time.t0 = 0;
time.dt = 0.002;
time.tf = 10;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);

%% ------------------------------------------------------------------------
% model
%--------------------------------------------------------------------------
sigma=0.05;
model.fn = 4;               % state space dimensionality
model.fx = 'fx_sinx';
model.mx = 'fx_meas';
model.hn = 2;               % measurement dimensionality
model.hx = 'hx_squared';
model.Q = blkdiag(0*eye(2),zeros(2))*time.dt;
model.sQ = sqrt(model.Q);
model.R = sigma^2*eye(2);

%% ------------------------------------------------------------------------
% prior Gaussian Sum   !!!!!! make them cells (mu & sig) !!!!!!!!
%--------------------------------------------------------------------------
%prior.n = 1;
prior.mu{1} = [-1 -1];
%prior.sig{1} = 0.01*eye(4);
% prior.mu{2} = 1;
% prior.sig{2} = 1;
% prior.weig = [.1 .9];
%prior.weig=1;

% prior dist. for params.

prior.p = [0.9 0.5 -1.45 0.5];

prior.p_t = [1.3572  0 -1.0658  0];

%% ------------------------------------------------------------------------
% measurements available every n time steps
%--------------------------------------------------------------------------
measTs = stv;
measAvail = zeros(1,time.nSteps);
measAvail(1:measTs:time.nSteps) = 1;

%% ------------------------------------------------------------------------
% particle filter settings
%--------------------------------------------------------------------------
pf.no_particles = 50000;
pf.no_bins = 100;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = inf;                  % treshold for the number of effective samples

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

for i = 1 : model.fn-model.hn
    xt_param = prior.p_t((i-1)*2+1) + prior.p_t(i*2)*rand(1,1);
    xt{1} = [xt{1}; xt_param];
end;

y = feval(model.mx, xt{1},time.tspan);
% get the trajectory of the sample over time
for k = 2 : time.nSteps
%     Y1(k)=y(1);
%     Y2(k)=y(2);
    xt{k} = y(k,:)' + model.sQ * randn(model.fn,1);
    y_meas{k} = feval(model.hx, xt{k}) + noise(k,:)';
end

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
X_pf_state = repmat(prior.mu{1}',1,pf.no_particles);
X_pf = X_pf_state;

for i = 1 : model.fn-model.hn
    X_pf_param = Pm(:,i)';
%     X_pf_param = prior.p((i-1)*2+1) + prior.p(i*2)*rand(1,pf.no_particles);
    X_pf = [X_pf; X_pf_param];
end;

w_pf = ones(1, pf.no_particles) / pf.no_particles;

%% ------------------------------------------------------------------------
% compute initial estimates
%--------------------------------------------------------------------------

[PF.x{1},PF.P{1},PF.minB{1},PF.maxB{1},tmp_X_pf] = getPFdata(X_pf,w_pf);

%% ------------------------------------------------------------------------
% RUN BOOSTRAP PARTICLE FILTER
%--------------------------------------------------------------------------   
M_x1=[];M_x2=[];M_p1=[];M_p2=[];
Nm=10;
for k = measTs+1 : measTs:time.nSteps
    fprintf('  - run filters tStep %d / %d\n', k, time.nSteps);    


%% ------------------------------------------------------------------------
% PF - propagate particles
%--------------------------------------------------------------------------
    fprintf('     - Particle Filter\n'); 
    
    % PF - time update
    [X_pf Momx1 Momx2 Momp1 Momp2] = PF_time_update1(X_pf, model, time, k, measTs,Nm);
%     keyboard
    % PF - measurement update    
    if (measAvail(k) == 1)
%         keyboard
        [X_pf, w_pf] = PF_meas_update(X_pf, w_pf, model, pf, y_meas{k});
    end
    % Updating Moments
    for j=1
        Momx1(end,j)=sum(X_pf(1,:))/pf.no_particles;
        Momx2(end,j)=sum(X_pf(2,:))/pf.no_particles;
        Momp1(end,j)=sum(X_pf(3,:))/pf.no_particles;
        Momp2(end,j)=sum(X_pf(4,:))/pf.no_particles;
    end

    for j=2:Nm
        Momx1(end,j)=sum((X_pf(1,:)-Momx1(end,1)).^j)/pf.no_particles;
        Momx2(end,j)=sum((X_pf(2,:)-Momx2(end,1)).^j)/pf.no_particles;
        Momp1(end,j)=sum((X_pf(3,:)-Momp1(end,1)).^j)/pf.no_particles;
        Momp2(end,j)=sum((X_pf(4,:)-Momp2(end,1)).^j)/pf.no_particles;
    end     
                
    % PF - compute estimates
    [PF.x{k},PF.P{k},PF.minB{k},PF.maxB{k},tmp_X_pf] = getPFdata(X_pf, w_pf);
    histdata{k}=tmp_X_pf;
    
    
    M_x1=[M_x1;Momx1];
    M_x2=[M_x2;Momx2];
    M_p1=[M_p1;Momp1];
    M_p2=[M_p2;Momp2];
end;         
delta_t=toc
keyboard




Nm=6;
pweight=1/pf.no_particles;
for i=1:4
    counter=0;
    for k = 1 : time.nSteps
        tempdata=histdata{k};
        Mom{k} = ComputeMoments(tempdata(i,:),pweight,1,Nm);
        if k==1;
            figure(i);subplot(2,3,1);
            hist(tempdata(i,:),100);title(['PF at t=',num2str(time.tspan(k))]);
            if i==1 
                xlabel('x')
%                 xlim([-1.4 -0.4])
            elseif i==2 
                xlabel('dx/dt')
%                 xlim([-1.5 1])
            elseif i==3
                xlabel('\eta')
                xlim([0.5 1.5])
            elseif i==4
                xlabel('\alpha')
                xlim([-1.5 -0.5])
            end
            if i>2
                line([prior.p_t(2*(i-2)-1) prior.p_t(2*(i-2)-1)],[0,500],'LineWidth',2,'Color','r')
            end
        end
        if mod(k,1000)==0
            counter=1+counter;
            figure(i);subplot(2,3,counter+1);
            hist(tempdata(i,:),100);title(['PF at t=',num2str(time.tspan(k+1))]);
            if i==1 && k~=1
                xlabel('x')
                xlim([-1.4 -0.4])
            elseif i==2 && k~=1
                xlabel('dx/dt')
                xlim([-1.5 1])
            elseif i==3
                xlabel('\eta')
                xlim([0.5 1.5])
            elseif i==4
                xlabel('\alpha')
                xlim([-1.5 -0.5])
            end
            if i>2
                line([prior.p_t(2*(i-2)-1) prior.p_t(2*(i-2)-1)],[0,500],'LineWidth',2,'Color','r')
            end
        end
    end
%     Mom1=cell2mat(Mom');
%     for j=1:Nm
%         figure(model.fn+i);subplot(2,3,j);
%         plot(time.tspan(2:end),Mom1(:,j))
%     end
end


for i=1:4
    counter=0;
    for k = 1 : time.nSteps
        tempdata=histdata{k};
        if k==1;
            figure(model.fn+2+i);subplot(2,3,1);
            hist(tempdata(i,:),100);title(['PF before Updat at t=',num2str(time.tspan(k))]);
            if i==1 
                xlabel('x')
%                 xlim([-1.4 -0.4])
            elseif i==2 
                xlabel('dx/dt')
%                 xlim([-1.5 1])
            elseif i==3
                xlabel('\eta')
                xlim([0.5 1.5])
            elseif i==4
                xlabel('\alpha')
                xlim([-1.5 -0.5])
            end
            if i>2
                line([prior.p_t(2*(i-2)-1) prior.p_t(2*(i-2)-1)],[0,500],'LineWidth',2,'Color','r')
            end
        end
        if mod(k,1000)==999
            counter=1+counter;
            figure(model.fn+2+i);subplot(2,3,counter+1);
            hist(tempdata(i,:),100);title(['PF before Update at t=',num2str(time.tspan(k+2))])
            
            if i==1 && k~=1
                xlabel('x')
                xlim([-1.4 -0.4])
            elseif i==2 && k~=1
                xlabel('dx/dt')
                xlim([-1.5 1])
            elseif i==3
                xlabel('\eta')
                xlim([0.5 1.5])
            elseif i==4
                xlabel('\alpha')
                xlim([-1.5 -0.5])
            end
            if i>2
                line([prior.p_t(2*(i-2)-1) prior.p_t(2*(i-2)-1)],[0,500],'LineWidth',2,'Color','r')
            end
        end
    end
end





%% ------------------------------------------------------------------------
% PLOT RESULTS
%--------------------------------------------------------------------------

%% - Measurements, Truth, Estimates ---------------------------------------
XP = cell2mat(xt);
PFP = cell2mat(PF.x);
PFminB = cell2mat(PF.minB);
PFmaxB = cell2mat(PF.maxB);

fig1=figure; hold on;box on
plot(time.tspan,XP(3,:),'k','LineWidth',2);
plot(time.tspan,PFP(3,:),'b','LineWidth',2);
plot(time.tspan,PFminB(3,:),'g--','LineWidth',2);
plot(time.tspan,PFmaxB(3,:),'r--','LineWidth',2);
legend('\eta_{act}','\eta_{mean}','\eta_{min}','\eta_{max}','Location','SouthEast');
% title('Estimation of Parameter \eta using Particle Filter');
xlabel('Time (sec)');ylabel('\eta');
ylim([0.5 1.5])
saveas(fig1,'etest.fig')
saveas(fig1,'etest.jpg')

fig2=figure; hold on;box on
plot(time.tspan,XP(4,:),'k','LineWidth',2);
plot(time.tspan,PFP(4,:),'b','LineWidth',2);
plot(time.tspan,PFminB(4,:),'g--','LineWidth',2);
plot(time.tspan,PFmaxB(4,:),'r--','LineWidth',2);
legend('\alpha_{act}','\alpha_{mean}','\alpha_{min}','\alpha_{max}','Location','SouthEast');
% title('Estimation of Parameter \alpha using Particle Filter');
xlabel('Time (sec)');ylabel('\alpha');
ylim([-2 -0.8])
saveas(fig2,'alest.fig')
saveas(fig2,'alest.jpg')


fig3=figure; hold on;box on
ym=cell2mat(y_meas);
plot(time.tspan,ym(1,:),'g.');
% plot(time.tspan,XP(1,:),'b','LineWidth',2);
plot(time.tspan,PFP(1,:),'k','LineWidth',2);
plot(time.tspan,PFminB(1,:),'b--','LineWidth',2);
plot(time.tspan,PFmaxB(1,:),'r--','LineWidth',2);
legend('x_{meas}','x_{mean}','x_{min}','x_{max}','Location','SouthEast');
xlabel('Time (sec)');ylabel('x(t)');
% ylim([-2 -0.8])
% title('Estimation of x(t) Using Particle Filter');
saveas(fig3,'xest.fig')
saveas(fig3,'xest.jpg')




fig5=figure; hold on;box on
ym=cell2mat(y_meas);
plot(time.tspan,ym(2,:),'g.');
% plot(time.tspan,XP(2,:),'b','LineWidth',2);
plot(time.tspan,PFP(2,:),'k','LineWidth',2);
plot(time.tspan,PFminB(2,:),'b--','LineWidth',2);
plot(time.tspan,PFmaxB(2,:),'r--','LineWidth',2);
legend('dx/dt_{meas}','dx/dt_{mean}','dx/dt_{min}','dx/dt_{max}','Location','SouthEast');
xlabel('Time (sec)');ylabel('dx/dt(t)');
% title('Estimation of dx/dt(t) Using Particle Filter');
saveas(fig5,'xdest.fig')
saveas(fig5,'xdest.jpg')



fig4=figure; hold on;box on
ym=cell2mat(y_meas);
% plot(time.tspan,ym(1,:),'g.');
% plot(time.tspan,XP(1,:),'b','LineWidth',2);
plot(time.tspan,PFP(1,:)-ye(:,1)','k','LineWidth',2);
plot(time.tspan,PFminB(1,:)-ye(:,1)','b--','LineWidth',2);
plot(time.tspan,PFmaxB(1,:)-ye(:,1)','r--','LineWidth',2);
err_mean=norm(PFP(1,:)-ye(:,1)');
err_min=norm(PFminB(1,:)-ye(:,1)');
err_max=norm(PFminB(1,:)-ye(:,1)');
legend(['error_{mean},|error_{mean}|_2=',num2str(err_mean)],['error_{min},|error_{min}|_2=',num2str(err_min)],['error_{max},|error_{max}|_2=',num2str(err_max)],'Location','NorthEast');
xlabel('Time (sec)');ylabel('error_x');
ylim([-0.1 0.1])
% title('Estimation of x(t) Using Particle Filter');
saveas(fig4,'xerr.fig')
saveas(fig4,'xerr.jpg')


fig6=figure; hold on;box on
ym=cell2mat(y_meas);
% plot(time.tspan,ym(1,:),'g.');
% plot(time.tspan,XP(1,:),'b','LineWidth',2);
plot(time.tspan,PFP(2,:)-ye(:,2)','k','LineWidth',2);
plot(time.tspan,PFminB(2,:)-ye(:,2)','b--','LineWidth',2);
plot(time.tspan,PFmaxB(2,:)-ye(:,2)','r--','LineWidth',2);
err_mean=norm(PFP(2,:)-ye(:,2)');
err_min=norm(PFminB(2,:)-ye(:,2)');
err_max=norm(PFminB(2,:)-ye(:,2)');
legend(['error_{mean},|error_{mean}|_2=',num2str(err_mean)],['error_{min},|error_{min}|_2=',num2str(err_min)],['error_{max},|error_{max}|_2=',num2str(err_max)],'Location','NorthEast');
xlabel('Time (sec)');ylabel('error_{dx/dt}');
ylim([-0.2 0.2])
saveas(fig6,'xderr.fig')
saveas(fig6,'xderr.jpg')