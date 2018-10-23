function y = rangeAndBearing(x)

xr = 6374*cosd(0);
yr = 6374*sind(0);
r = sqrt((x(1) - xr).^2 + (x(2) - yr).^2);
theta = atan2((x(2) - yr),(x(1) - xr));
y(1,1) = r;
%  y(2,1) = theta;

%  y(1,1) = theta;
end
