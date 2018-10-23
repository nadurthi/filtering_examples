function [x,y]=point_ellipse(xc,yc,a,b)
i=1;
for t=0:5:360
    x(i)=xc+a*cosd(t);
    y(i)=yc+b*sind(t);
i=i+1;
end
x(i)=xc;
y(i)=yc;