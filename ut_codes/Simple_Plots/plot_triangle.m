function plot_triangle(xc,yc,a,b,th,c,meth)
% th is the rotation matrix
% c is the colour 
% xc,yc is the center, a is the radius of cricle that circumscribes this
% equilateral triangle
% b is the skew parameter that stretches the triagle(non transformed) along the y axis
% b=1 it is equilateral  
R=[cos(th),-sin(th);sin(th),cos(th)];
x=[0,a*cosd(30),-a*cosd(30),0];
y=[a*b,-a*sind(30),-a*sind(30),a*b];
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
