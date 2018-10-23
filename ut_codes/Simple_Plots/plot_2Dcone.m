function plot_2Dcone(xc,yc,alphaa,a,th,c,meth)
% th is the rotation matrix
% c is the colour 
% xc,yc is the vertex of cone
% equilateral triangle
% alpha is half angle, a is the center range 
% as the triag is plotted is upwar
th=th-pi/2;
R=[cos(th),-sin(th);sin(th),cos(th)];
x=[0,a*tan(alphaa),-a*tan(alphaa),0];
y=[0,a,a,0];
% keyboard
for i=1:1:length(x)
     xy=R*[x(i);y(i)]+[xc;yc];
x(i)=xy(1);
y(i)=xy(2);
end
plot(x,y,c,'linewidth',2)
if strcmpi(meth,'fill')
    fill(x,y,c)
    alpha 0.2;
end
