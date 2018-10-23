function plot_ellipse(xc,yc,a,b,th,c,meth)
% R is the rotation matrix
% c is the colour 
% xc,yc is the center, a,b are the major and minor axis
i=1;
 R=[cos(th),-sin(th);sin(th),cos(th)];
for t=0:5:360
    xy=R*[a*cosd(t);b*sind(t)]+[xc;yc];
x(i)=xy(1);
y(i)=xy(2);
    i=i+1;
end
plot(x,y,c)
if strcmpi(meth,'fill')
    fill(x,y,c)
    alpha 0.2;
end