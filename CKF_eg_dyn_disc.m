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