function y = fx_meas(x,tvec)
% keyboard
opt = odeset('reltol',1e-12,'abstol',1e-12);
% keyboard
[t y]=ode45(@fx,tvec,x,opt);
% y=y(end,:);


function dy=fx(t,y)

dy=zeros(4,1);
Uin=sin(3*t);
dy(1)=y(2);
dy(2)=Uin-y(3)*y(2)-y(4)*y(1)-2*y(1)^3;
dy(3)=0;
dy(4)=0;
% dy(5)=0;