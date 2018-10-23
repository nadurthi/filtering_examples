clc
clear
% this codes test the new method sigma points with UT 2n+1,4n+1,6n+1, GH,
% montecarlo, when the dynamics is polynomial nonlinearity.
%% define your system dynamics F
F=@(t,x)[x(2)^3;-2*x(2)^2-3*x(1)^3-4*x(1)^2*x(2)^1];
%% ode 45 solution for mo
[t,x]=ode45(F,[0,15],[1,0]);
plot(x(:,1),x(:,2))

%% discrete system
x(1)=0
for i=1:1:100
    x(i)=