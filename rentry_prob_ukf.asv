%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Kumar Vishwajeet
% contact: kumarvis@buffalo.edu
% date:    13rd April,2011
% problem:- Vehicle Re-entry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
%close all;
%% time
time.t0 = 0;
time.dt = 0.05;
time.tf = 200;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);

%% system model
model.fn = 5;               % state space dimensionality
model.fx = 'reEntryDynamics';
model.hn = 2;               % measurement dimensionality
model.hx = 'rangeAndBearing';
model.Q = blkdiag(0,0,2.4064e-5,2.4064e-5,0);
model.sQ = sqrtm(model.Q);
model.R = blkdiag(1e-3,0.017);

%% Measurement
measTs = 0.1;
measAvail = zeros(1,time.nSteps);
measAvail(1:measTs:time.nSteps) = 1;

%% Weights for UKF points (Refer: Dr. Crassidis Book)
alpha = 1; 
L = model.fn;
beta = 2;kappa = 0;
lambda = alpha*alpha*(L+kappa) - L;
xeta = sqrt(L + lambda);
Wm0 = lambda/(lambda+L);
Wc0 = lambda/(lambda+L)+ (1-alpha*alpha+beta);
Wm = [Wm0 1/(2*(L+lambda)) 1/(2*(L+lambda))];
Wc = [Wc0 1/(2*(L+lambda)) 1/(2*(L+lambda))];

%% Initial conditions for Monte carlo and Kalman filter
x_true = [6500.4 349.14 -1.8093 -6.7967 0.6932]';
P_true = blkdiag(1e-6,1e-6,1e-6,1e-6,0);
x_filter = [6500.4 349.14 -1.8093 -6.7967 0]';
P_filter = blkdiag(1e-6,1e-6,1e-6,1e-6,1);
meu = x_filter;
sig = P_filter;

%% Monte Carlo simulation
N = 1; % monte carlo runs
t0 = 0;
dt = 0.05; % 50 millisecond
tf = 200;  % final time
tt=t0:dt:tf; %time of integration

options=odeset('RelTol',1e-10,'AbsTol',1e-10);
Rmonte = zeros(tf/dt+1,model.fn,N);
saveUKFMeu = zeros(model.fn,time.nSteps,N);
saveUKFCov = zeros(model.fn,model.fn,time.nSteps,N);
for ct = 1: N
    x10 = repmat(x_true(1:4),1,N) + 0.001*randn(model.fn-1,N);
    x20 = repmat(x_true(5),1,N);
    x0 = [x10;x20];
    [~,Rosc]=ode45(@(t,x) reEntryDynamics(t,x),tt,x0(:,ct),options);
    Rmonte(:,:,ct) = Rosc;
    
    % Create Truth
    k=1;
    xt(:,1,ct) = x0(:,1);
    y_meas(:,k) = feval(model.hx,x0(:,1)) + chol(model.R)'* randn(model.hn,1);
    
    % get the trajectory of the sample over time
    saveUKFMeu(:,1,ct) = x_filter;
    saveUKFCov(:,:,k,ct) = P_filter;
    meu = x_filter;
    sig = P_filter;
    sig_up = P_filter;
    for k = 2 : time.nSteps
        tInt = [time.tspan(k-1) time.tspan(k)];
        [~,y]=ode45(@(t,x) reEntryDynamics(t,x),tInt,xt(:,k-1,ct),options);
        xt(:,k,ct) = y(end,:)' + model.sQ*eye(model.fn,1).* randn(model.fn,1);
        y_meas(:,k) = feval(model.hx, xt(:,k,ct)) + chol(model.R)'* randn(model.hn,1);
    
        % Unscented Kalman Filter      
        tInt = [time.tspan(k-1) time.tspan(k)];
        fprintf('MC Run = %d/%d \t Unscented Kalman Filter\t',ct,N);
        fprintf('Time step %d/%d\n',k,time.nSteps);
        
        % time update
        [meu, sig, sigPoints] = feval('UKF_time_update', model, meu, Wm,Wc,xeta,sig,tInt,options);    
        % measurement update
        if (measAvail(k) == 1)
            [meu_up, sig_up] = feval('UKF_measurement_update', model,meu,Wm,Wc,sig,y_meas(:,k),sigPoints);
            meu = meu_up;
            sig = sig_up;
        else
            meu_up = meu;
            sig_up = sig;                
        end
        saveUKFMeu(:,k,ct) = meu_up;
        saveUKFCov(:,:,k,ct) = sig_up;
    end
end

for i = 1 : time.nSteps
    for j = 1:model.fn
        mu(i,j) = sum(saveUKFMeu(j,i,:))./N;
        for k = 1:model.fn
            sigma(j,k,i) = sum(saveUKFCov(j,k,i,:))./N;
        end
    end
end
save 