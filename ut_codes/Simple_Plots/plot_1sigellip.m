function plot_1sigellip(mu,P,c,W)

switch nargin
    case 3
        w = 1;
    case 4
        w = W;
end

i=1;
for th=0:0.5:360
    XY(i,:)=(sqrtm(P)*[cosd(th);sind(th)]+mu(:))';
    i=i+1;
end
if isreal(XY)==0 
%     keyboard
end

    plot(XY(:,1),XY(:,2),c,'linewidth',w)

end