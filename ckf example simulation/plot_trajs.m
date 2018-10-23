function plot_trajs
t0=0;
 dt=1;
tF=100;
for i=-4:1:4
x0tr=[1000,300,1000,0,(-3+i)*pi/180]';
[t,x_mc]=ode45_discc(@CKF_eg_dyn_disc,t0,dt,tF,x0tr,1e-200);
hold on
plot(x_mc(:,1),x_mc(:,3))
pause(1)
end
end
    function xk1=CKF_eg_dyn_disc(xk)
% global T
T=1;
omg=xk(end,1);
xk1=[1,sin(omg*T)/omg,0,-(1-cos(omg*T))/omg,0;...
     0,cos(omg*T),0,-sin(omg*T),0;
     0,(1-cos(omg*T))/omg,1,sin(omg*T)/omg,0;
     0,sin(omg*T),0,cos(omg*T),0;
     0,0,0,0,1]*xk;

end