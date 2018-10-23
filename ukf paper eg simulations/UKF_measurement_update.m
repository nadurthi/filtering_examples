function [omu, osig] = UKF_measurement_update(model,imu,Wm,Wc,isig,y,sigPoints)
L = numel(imu);

%% Calculate mean
for i = 1:2*L+1
    Y(:,i) = feval(model.hx,sigPoints(:,i));
end

sum2 = zeros(model.hn,1);
sum3 = zeros(model.hn,1);
for i = 2:L+1
    sum2 = sum2+[Y(:,i)];
    sum3 = sum3+[Y(:,i+L)];
end

yMean = Wm(1).*Y(:,1)+Wm(2).*sum2+Wm(3).*sum3;

%% Calculate Covariance
sum2 = zeros(model.hn);
sum3 = zeros(model.hn);
for i = 2:L+1
    sum2 = sum2 + ([Y(:,i)]-yMean)*([Y(:,i)]-yMean)';
    sum3 = sum3 + ([Y(:,i+L)]-yMean)*([Y(:,i+L)]-yMean)';
end
Py = Wc(1).*([Y(:,1)]-yMean)*([Y(:,1)]-yMean)' + Wc(2).*sum2 + Wc(3).*sum3 + model.R;

sum2 = zeros(L,model.hn);
sum3 = zeros(L,model.hn);
for i = 2:L+1
    sum2 = sum2 + (sigPoints(:,i)-imu)*([Y(:,i)]-yMean)';
    sum3 = sum3 + (sigPoints(:,i+L)-imu)*([Y(:,i+L)]-yMean)';
end
Pxy = Wc(1).*(sigPoints(:,1)-imu)*([Y(:,1)]-yMean)' + Wc(2).*sum2 + Wc(3).*sum3;
Kk = Pxy*inv(Py);
omu = imu + Kk*(y - yMean);
osig = isig - Kk*Py*Kk';