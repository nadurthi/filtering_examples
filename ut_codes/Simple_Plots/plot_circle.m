function plot_circle(xc,yc,r)
% R is the rotation matrix
% c is the colour 
% xc,yc is the center, a,b are the major and minor axis
i=1;
for t=0:5:360
    xy=[r*cosd(t);r*sind(t)]+[xc;yc];
x(i)=xy(1);
y(i)=xy(2);
    i=i+1;
end
plot(x,y)

end