% simulating the examples in UKF paper
function spherical_to_polar
clc
clear
global kappa
kappa=0;
%% Example 1 - (r,th)- (x,y) transformation
% 
% first by montecarlo method
sigr=(2)^2;
sigth=(15)^2;
sigphi=(15)^2;
mur=1;
muth=90;
muphi=45;
F=@(x)([x(:,1).*sind(x(:,2)).*cosd(x(:,3)),x(:,1).*sind(x(:,2)).*sind(x(:,3)),x(:,1).*cosd(x(:,2))]);
% N=3*10^6;
% x=mvnrnd([mur;muth],[sigr,0;0,sigth],N);
% w=(1/N)*ones(N,1);
[x,w]=monte_carlo_int_normal([mur;muth;muphi],[sigr,0,0;0,sigth,0;0,0,sigphi],10^6);
XY=F(x);
p_XY=XY;
mu_mc(1)=0;
mu_mc(2)=0;
mu_mc(3)=0;
for i=1:1:length(w)
    mu_mc(1)=mu_mc(1)+w(i)*XY(i,1);
    mu_mc(2)=mu_mc(2)+w(i)*XY(i,2);
    mu_mc(3)=mu_mc(3)+w(i)*XY(i,3);
end
P_mc=0;
for i=1:1:length(w)
    P_mc=P_mc+w(i)*(XY(i,:)'-mu_mc')*(XY(i,:)-mu_mc);
end
[xmc,ymc,zmc]=point_ellipsoid(mu_mc(1),mu_mc(2),mu_mc(3),sqrt(P_mc(1,1)),sqrt(P_mc(2,2)),sqrt(P_mc(3,3)));

save xmc
save ymc
save zmc

% UT 2n+1 points
[x,w]=UT_sigmapoints([mur;muth;muphi],[sigr,0,0;0,sigth,0;0,0,sigphi],2);
XY=F(x);
mu_ut(1)=0;
mu_ut(2)=0;
mu_ut(3)=0;
for i=1:1:length(w)
    mu_ut(1)=mu_ut(1)+w(i)*XY(i,1);
    mu_ut(2)=mu_ut(2)+w(i)*XY(i,2);
    mu_ut(3)=mu_ut(3)+w(i)*XY(i,3);
end
P_ut=0;
for i=1:1:length(w)
    P_ut=P_ut+w(i)*(XY(i,:)'-mu_ut')*(XY(i,:)-mu_ut);
end

[xut,yut,zut]=point_ellipsoid(mu_ut(1),mu_ut(2),mu_ut(3),sqrt(P_ut(1,1)),sqrt(P_ut(2,2)),sqrt(P_ut(3,3)));
save xut
save yut
save zut

% CKF points
[x,w]=cubature_KF_points([mur;muth;muphi],[sigr,0,0;0,sigth,0;0,0,sigphi]);
XY=F(x);
mu_ckf(1)=0;
mu_ckf(2)=0;
mu_ckf(3)=0;
for i=1:1:length(w)
    mu_ckf(1)=mu_ckf(1)+w(i)*XY(i,1);
    mu_ckf(2)=mu_ckf(2)+w(i)*XY(i,2);
    mu_ckf(3)=mu_ckf(3)+w(i)*XY(i,3);
end
P_ckf=0;
for i=1:1:length(w)
    P_ckf=P_ckf+w(i)*(XY(i,:)'-mu_ckf')*(XY(i,:)-mu_ckf);
end

[xckf,yckf,zckf]=point_ellipsoid(mu_ckf(1),mu_ckf(2),mu_ckf(3),sqrt(P_ckf(1,1)),sqrt(P_ckf(2,2)),sqrt(P_ckf(3,3)));
save xckf
save yckf
save zckf

% NM 4thmom points
[x,w]=conjugate_dir_gausspts([mur;muth;muphi],[sigr,0,0;0,sigth,0;0,0,sigphi]);
XY=F(x);
mu_nm4(1)=0;
mu_nm4(2)=0;
mu_nm4(3)=0;
for i=1:1:length(w)
    mu_nm4(1)=mu_nm4(1)+w(i)*XY(i,1);
    mu_nm4(2)=mu_nm4(2)+w(i)*XY(i,2);
    mu_nm4(3)=mu_nm4(3)+w(i)*XY(i,3);
end
P_nm4=0;
for i=1:1:length(w)
    P_nm4=P_nm4+w(i)*(XY(i,:)'-mu_nm4')*(XY(i,:)-mu_nm4);
end

[xnm4,ynm4,znm4]=point_ellipsoid(mu_nm4(1),mu_nm4(2),mu_nm4(3),sqrt(P_nm4(1,1)),sqrt(P_nm4(2,2)),sqrt(P_nm4(3,3)));
save xnm4
save ynm4
save znm4


% NM 6thmom points
[x,w]=conjugate_dir_gausspts_till_6moment_scheme2([mur;muth;muphi],[sigr,0,0;0,sigth,0;0,0,sigphi]);
XY=F(x);
mu_nm6(1)=0;
mu_nm6(2)=0;
mu_nm6(3)=0;
for i=1:1:length(w)
    mu_nm6(1)=mu_nm6(1)+w(i)*XY(i,1);
    mu_nm6(2)=mu_nm6(2)+w(i)*XY(i,2);
    mu_nm6(3)=mu_nm6(3)+w(i)*XY(i,3);
end
P_nm6=0;
for i=1:1:length(w)
    P_nm6=P_nm6+w(i)*(XY(i,:)'-mu_nm6')*(XY(i,:)-mu_nm6);
end


[xnm6,ynm6,znm6]=point_ellipsoid(mu_nm6(1),mu_nm6(2),mu_nm6(3),sqrt(P_nm6(1,1)),sqrt(P_nm6(2,2)),sqrt(P_nm6(3,3)));
save xnm6
save ynm6
save znm6 


load xmc
load ymc
load zmc

load xut
load yut
load zut

load xckf
load yckf
load zckf

load xnm4
load ynm4
load znm4

load xnm6
load ynm6
load znm6

% figure(1)
% [r,th,phi]=point_ellipsoid(mur,muth,muth,sqrt(sigr),sqrt(sigth),sqrt(sigphi));
% plot3(r,th,phi,'b.')
% legend('r-th-phi')
% xlabel('r')
% ylabel('th')
% zlabel('phi')
% set(gca,'FontSize',16)
% title('Polar (cm,rad)')
% h = get(gca, 'title');
% k = get(gca, 'xlabel');
% l = get(gca, 'ylabel');
% set(h, 'FontName', 'Helvetica', 'FontSize', 16)
% set(k, 'FontName', 'Helvetica', 'FontSize', 16)
% set(l, 'FontName', 'Helvetica', 'FontSize', 16)
figure(2)
plot3(xmc,ymc,zmc,'rd',xut,yut,zut,'ks',xckf,yckf,zckf,'b+',xnm4,ynm4,znm4,'go',xnm6,ynm6,znm6,'m*')
% n=length(0:10:360);
% hold on
% plot(xmc(1:1:n),ymc(1:1:n),'rd',xut(1:1:n),yut(1:1:n),'ks',xckf(1:1:n),yckf(1:1:n),'b+',xnm4(1:1:n),ynm4(1:1:n),'go',xnm6(1:1:n),ynm6(1:1:n),'m*')
% plot(xmc(n+1:1:2*n),ymc(n+1:1:2*n),'rd',xut(n+1:1:2*n),yut(n+1:1:2*n),'ks',xckf(n+1:1:2*n),yckf(n+1:1:2*n),'b+',xnm4(n+1:1:2*n),ynm4(n+1:1:2*n),'go',xnm6(n+1:1:2*n),ynm6(n+1:1:2*n),'m*')
% plot(xmc(2*n+1:1:3*n),ymc(2*n+1:1:3*n),'rd',xut(2*n+1:1:3*n),yut(2*n+1:1:3*n),'ks',xckf(2*n+1:1:3*n),yckf(2*n+1:1:3*n),'b+',xnm4(2*n+1:1:3*n),ynm4(2*n+1:1:3*n),'go',xnm6(2*n+1:1:3*n),ynm6(2*n+1:1:3*n),'m*')
% plot(xmc(3*n+1:1:4*n),ymc(3*n+1:1:4*n),'rd',xut(3*n+1:1:4*n),yut(3*n+1:1:4*n),'ks',xckf(3*n+1:1:4*n),yckf(3*n+1:1:4*n),'b+',xnm4(3*n+1:1:4*n),ynm4(3*n+1:1:4*n),'go',xnm6(3*n+1:1:4*n),ynm6(3*n+1:1:4*n),'m*')
% plot(xmc(4*n+1:1:5*n),ymc(4*n+1:1:5*n),'rd',xut(4*n+1:1:5*n),yut(4*n+1:1:5*n),'ks',xckf(4*n+1:1:5*n),yckf(4*n+1:1:5*n),'b+',xnm4(4*n+1:1:5*n),ynm4(4*n+1:1:5*n),'go',xnm6(4*n+1:1:5*n),ynm6(4*n+1:1:5*n),'m*')
% plot(xmc(5*n+1:1:6*n),ymc(5*n+1:1:6*n),'rd',xut(5*n+1:1:6*n),yut(5*n+1:1:6*n),'ks',xckf(5*n+1:1:6*n),yckf(5*n+1:1:6*n),'b+',xnm4(5*n+1:1:6*n),ynm4(5*n+1:1:6*n),'go',xnm6(5*n+1:1:6*n),ynm6(5*n+1:1:6*n),'m*')

legend('MC','UT','CKF','NM4','NM6')
xlabel('x')
ylabel('y')
% zlabel('z')
set(gca,'FontSize',16)
title('Cartesian (cm,cm)')
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 16)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
% figure(3)
% plot3(p_XY(:,1),p_XY(:,2),p_XY(:,3),'b+')
end
function [x,y,z]=point_ellipsoid(mux,muy,muz,sdx,sdy,sdz)
i=1;
for th=0:10:90
    for phi=0:10:360
       x(i)=mux+sdx*sind(th)*cosd(phi);
       y(i)=muy+sdy*sind(th)*sind(phi);
       z(i)=muz+sdz*cosd(th);
       i=i+1;
    end
end

end