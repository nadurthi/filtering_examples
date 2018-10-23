function [omu, osig, sigPoints] = UKF_time_update(model,imu,Wm,Wc,xeta,isig,T,opt)

L = numel(imu);
meu=[];
for i = 1:L
    meu = [meu;imu(i)];
end
%% Calculate mean
sig_sqrt = chol(isig)';
[~,x] = ode45(@UKF_propagation,T,meu(1:L),opt,model);
meu1{1} = x(end,:);

for i = 1:L
    [~,x] = ode45(@UKF_propagation,T,meu(1:L) + real(xeta.*sig_sqrt(1:L,i)),opt,model);
    meu2{i} = x(end,:);
    [~,x] = ode45(@UKF_propagation,T,meu(1:L) - real(xeta.*sig_sqrt(1:L,i)),opt,model);
    meu3{i} = x(end,:);
end

%% Calculate Covariance
sum2 = zeros(L,1);
sum3 = zeros(L,1);
for i = 1:L
    sum2 = sum2+[meu2{i}]';
    sum3 = sum3+[meu3{i}]';
end
omu = Wm(1).*[meu1{1}]'+Wm(2).*sum2+Wm(3).*sum3;

sum2 = zeros(L);
sum3 = zeros(L);
for i = 1:L
    sum2 = sum2 + ([meu2{i}]'-omu)*([meu2{i}]'-omu)';
    sum3 = sum3 + ([meu3{i}]'-omu)*([meu3{i}]'-omu)';
end
osig = Wc(1).*([meu1{1}]'-omu)*([meu1{1}]'-omu)' + Wc(2).*sum2 + Wc(3).*sum3+model.Q;
osig = osig(1:model.fn,1:model.fn);
%% generate new sigma points
sig_sqrt = chol(osig)';
meu2 = zeros(5,5);
meu3 = zeros(5,5);
for i = 1:L
    meu2(:,i) = omu(1:L) + real(xeta.*sig_sqrt(1:L,i));
    meu3(:,i) = omu(1:L) - real(xeta.*sig_sqrt(1:L,i));
end

sigPoints = [omu meu2 meu3];
