for n=2:8
[Xg2,Wg2]=Sigma_Quadrature_points(-1*ones(1,n),ones(1,n),4,'cut','uniform');
Nc(n-1)=size(Xg2,1);
end
Ng=3.^[2:8];
% N=[size(Xg2,1),size(Xg3,1),size(Xg4,1),size(Xg5,1);size(X2,1),size(X3,1),size(X4,1),size(X5,1)]

figure
bar3(2:1:8,[Ng',Nc'])
set(gca,'XTickLabel',{'GLgn3','CUT4-U'})
view(119,24)
xlabel('Dimension')
ylabel('Method')
zlabel('Number of points')
plot_prop_paper