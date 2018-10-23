%%%%%%%%%%%%%
%Understanding statisticla liner regression
%nonlinear regression
%EKF
%UKF-CKF-GH
%%%%%%%%%%%%%
clc
clear all
close all
%%%%%%%%
%% The measurement model
h=@(xk)(xk.^4);
hJ=@(xk)(4*xk.^3);
%% Initial condition P(xk)
muk=0.1;
stddk=1;

%% The MC true solution
xk_mc=normrnd(muk,stddk,1e5,1);
yk_mc=h(xk_mc);
plot(xk_mc,yk_mc,'*');
hold on
para_k1_mc=[mean(yk_mc),var(yk_mc)]

%% The EKF solution
para_k1_ekf=[h(muk),hJ(muk)^2*stddk^2]

%% Linear regression
xk_lr=normrnd(muk,stddk,10^5,1);
yk_lr=h(xk_lr);
A=[xk_lr,ones(size(xk_lr))];
B=yk_lr;
ab=A\B;
a=ab(1);
b=ab(2);
para_k1_lr=[a*muk+b,a^2*stddk^2]
plot(xk_lr,a*xk_lr+b,'r')
%% GH result
[xgh,wgh]=GH_points(muk,stddk^2,2);
mugh=sum(wgh.*h(xgh));
para_k1_gh=[mugh,sum(wgh.*(h(xgh)-mugh).^2)]

%% Linear regression with GH points
[xklrgh,wlrgh]=GH_points(muk,stddk^2,2);
yklrgh=h(xklrgh);
A=[xklrgh,ones(size(xklrgh))];
B=yklrgh;
ab=A\B;
a=ab(1);
b=ab(2);
para_k1_lrgh=[a*muk+b,a^2*stddk^2]

% %% Quadratic regression
% A=[xk_lr.^2,xk_lr,ones(size(xk_lr))];
% B=yk_lr;
% abc=A\B;
% a=abc(1);
% b=abc(2);
% c=abc(3);
% para_k1_lr=[a*(stddk^2+muk^2)+b*muk+c,a^2*stddk^2]