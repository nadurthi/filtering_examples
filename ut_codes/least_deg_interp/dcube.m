function corners = dcube(d)
%DCUBE the 2^d corners of the unit d-cube
%
%        corners = dcube(d)
%
% returns the (d,2^d) array CORNERS containing all the vertices of the unit
% cube in R^d .

% cb: 26jul97

corners = [1 -1]; pone = [1];
for j=2:d
   pone = [pone pone];
   corners = [pone -pone; corners corners];
end
