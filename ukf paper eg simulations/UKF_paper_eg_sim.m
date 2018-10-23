% simulating the examples in UKF paper
clc
clear
%% Example 1 - (r,th)- (x,y) transformation
% 
% first by montecarlo method
sigr=(0.02)^2;
sigth=(30*pi/180)^2;
mur=50;
muth=0;
%% Analytical answer
%means
lam1=exp(-sigth/2);
lam2=0.5*(1+exp(-2*sigth));
lam3=0.5*(1-exp(-2*sigth));
x=mur*cos(muth);
y=mur*sin(muth);
amux=lam1*x
amuy=lam1*y
asigx=(lam1^2*x^2-2*lam1^2*x^2+x^2*lam2+(sigr*x^2*lam2)/(x^2+y^2)+y^2*lam3+(sigr*y^2*lam3)/(x^2+y^2))
asigy=(lam1^2*y^2-2*lam1^2*y^2+y^2*lam2+(sigr*y^2*lam2)/(x^2+y^2)+x^2*lam3+(sigr*x^2*lam3)/(x^2+y^2))
asigxy=(lam1^2*x*y-2*lam1^2*x*y+x*y*lam2+(sigr*x*y*lam2)/(x^2+y^2)-x*y*lam3-(sigr*x*y*lam3)/(x^2+y^2))


%%

F=@(x)([x(:,1).*cos(x(:,2)),x(:,1).*sin(x(:,2))]);
% N=3*10^6;
% x=mvnrnd([mur;muth],[sigr,0;0,sigth],N);
% w=(1/N)*ones(N,1);
no_mc_samps=3e4;
[x,w]=monte_carlo_int_normal([mur;muth],[sigr,0;0,sigth],no_mc_samps);
x_mcm=x;
XY=F(x);
p_XY=XY;
mu_mc(1)=0;
mu_mc(2)=0;
for i=1:1:length(w)
    mu_mc(1)=mu_mc(1)+w(i)*XY(i,1);
    mu_mc(2)=mu_mc(2)+w(i)*XY(i,2);
end
P_mc=0;
for i=1:1:length(w)
    P_mc=P_mc+w(i)*(XY(i,:)'-mu_mc')*(XY(i,:)-mu_mc);
end
[xmc,ymc]=point_ellipse(mu_mc(1),mu_mc(2),sqrt(P_mc(1,1)),sqrt(P_mc(2,2)));
XYmc=XY;
wmc=ones(no_mc_samps,1)./no_mc_samps;

% UT 2n+1 points
[x,w]=UT_sigmapoints([mur;muth],[sigr,0;0,sigth],2);
XY=F(x);
mu_ut(1)=0;
mu_ut(2)=0;
for i=1:1:length(w)
    mu_ut(1)=mu_ut(1)+w(i)*XY(i,1);
    mu_ut(2)=mu_ut(2)+w(i)*XY(i,2);
end
P_ut=0;
for i=1:1:length(w)
    P_ut=P_ut+w(i)*(XY(i,:)'-mu_ut')*(XY(i,:)-mu_ut);
end
XYut=XY;
wut=w;
[xut,yut]=point_ellipse(mu_ut(1),mu_ut(2),sqrt(P_ut(1,1)),sqrt(P_ut(2,2)));

% CKF points
[x,w]=cubature_KF_points([mur;muth],[sigr,0;0,sigth]);
XY=F(x);
mu_ckf(1)=0;
mu_ckf(2)=0;
for i=1:1:length(w)
    mu_ckf(1)=mu_ckf(1)+w(i)*XY(i,1);
    mu_ckf(2)=mu_ckf(2)+w(i)*XY(i,2);
end
P_ckf=0;
for i=1:1:length(w)
    P_ckf=P_ckf+w(i)*(XY(i,:)'-mu_ckf')*(XY(i,:)-mu_ckf);
end

[xckf,yckf]=point_ellipse(mu_ckf(1),mu_ckf(2),sqrt(P_ckf(1,1)),sqrt(P_ckf(2,2)));

% NM 4thmom points
[x,w]=conjugate_dir_gausspts([mur;muth],[sigr,0;0,sigth]);
XY=F(x);
mu_nm4(1)=0;
mu_nm4(2)=0;
for i=1:1:length(w)
    mu_nm4(1)=mu_nm4(1)+w(i)*XY(i,1);
    mu_nm4(2)=mu_nm4(2)+w(i)*XY(i,2);
end
P_nm4=0;
for i=1:1:length(w)
    P_nm4=P_nm4+w(i)*(XY(i,:)'-mu_nm4')*(XY(i,:)-mu_nm4);
end
XYcut4=XY;
wcut4=w;
[xnm4,ynm4]=point_ellipse(mu_nm4(1),mu_nm4(2),sqrt(P_nm4(1,1)),sqrt(P_nm4(2,2)));


% NM 6thmom points
[x,w]=conjugate_dir_gausspts_till_6moment_scheme2([mur;muth],[sigr,0;0,sigth]);
XY=F(x);
mu_nm6(1)=0;
mu_nm6(2)=0;
for i=1:1:length(w)
    mu_nm6(1)=mu_nm6(1)+w(i)*XY(i,1);
    mu_nm6(2)=mu_nm6(2)+w(i)*XY(i,2);
end
P_nm6=0;
for i=1:1:length(w)
    P_nm6=P_nm6+w(i)*(XY(i,:)'-mu_nm6')*(XY(i,:)-mu_nm6);
end
XYcut6=XY;
wcut6=w;
[xnm6,ynm6]=point_ellipse(mu_nm6(1),mu_nm6(2),sqrt(P_nm6(1,1)),sqrt(P_nm6(2,2)));

% %GH4 thmom points
% [x,w]=GH_points([mur;muth],[sigr,0;0,sigth],4);
% XY=F(x);
% mu_gh4(1)=0;
% mu_gh4(2)=0;
% for i=1:1:length(w)
%     mu_gh4(1)=mu_gh4(1)+w(i)*XY(i,1);
%     mu_gh4(2)=mu_gh4(2)+w(i)*XY(i,2);
% end
% P_gh4=0;
% for i=1:1:length(w)
%     P_gh4=P_gh4+w(i)*(XY(i,:)'-mu_gh4')*(XY(i,:)-mu_gh4);
% end
% wcut8=w;
% XYcut8=XY;
% [xgh4,ygh4]=point_ellipse(mu_gh4(1),mu_gh4(2),sqrt(P_gh4(1,1)),sqrt(P_gh4(2,2)));


%GH3 thmom points
[x,w]=GH_points([mur;muth],[sigr,0;0,sigth],3);
XY=F(x);
mu_gh3(1)=0;
mu_gh3(2)=0;
for i=1:1:length(w)
    mu_gh3(1)=mu_gh3(1)+w(i)*XY(i,1);
    mu_gh3(2)=mu_gh3(2)+w(i)*XY(i,2);
end
P_gh3=0;
for i=1:1:length(w)
    P_gh3=P_gh3+w(i)*(XY(i,:)'-mu_gh3')*(XY(i,:)-mu_gh3);
end
wgh3=w;
XYgh3=XY;
[xgh3,ygh3]=point_ellipse(mu_gh3(1),mu_gh3(2),sqrt(P_gh3(1,1)),sqrt(P_gh3(2,2)));

% 
%GH4 thmom points
[x,w]=GH_points([mur;muth],[sigr,0;0,sigth],4);
XY=F(x);
mu_gh4(1)=0;
mu_gh4(2)=0;
for i=1:1:length(w)
    mu_gh4(1)=mu_gh4(1)+w(i)*XY(i,1);
    mu_gh4(2)=mu_gh4(2)+w(i)*XY(i,2);
end
P_gh4=0;
for i=1:1:length(w)
    P_gh4=P_gh4+w(i)*(XY(i,:)'-mu_gh4')*(XY(i,:)-mu_gh4);
end
wgh4=w;
XYgh4=XY;
[xgh4,ygh4]=point_ellipse(mu_gh4(1),mu_gh4(2),sqrt(P_gh4(1,1)),sqrt(P_gh4(2,2)));

%GH5 thmom points
[x,w]=GH_points([mur;muth],[sigr,0;0,sigth],5);
XY=F(x);
mu_gh5(1)=0;
mu_gh5(2)=0;
for i=1:1:length(w)
    mu_gh5(1)=mu_gh5(1)+w(i)*XY(i,1);
    mu_gh5(2)=mu_gh5(2)+w(i)*XY(i,2);
end
P_gh5=0;
for i=1:1:length(w)
    P_gh5=P_gh5+w(i)*(XY(i,:)'-mu_gh5')*(XY(i,:)-mu_gh5);
end
wgh5=w;
XYgh5=XY;
[xgh5,ygh5]=point_ellipse(mu_gh5(1),mu_gh5(2),sqrt(P_gh5(1,1)),sqrt(P_gh5(2,2)));


figure(1)
[r,th]=point_ellipse(mur,muth,sqrt(sigr),sqrt(sigth));
plot(x_mcm(1:10^4,1),x_mcm(1:10^4,2),'ro',r,th,'b.','LineWidth',2)
legend('r-th')
xlabel('r')
ylabel('th')
set(gca,'FontSize',16)
title('Polar (cm,rad)')
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 16)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

 figure(2)
plot(xmc,ymc,'kd',xut,yut,'ks',xckf,yckf,'k^',xgh3,ygh3,'k+',xgh4,ygh4,'kp',xnm4,ynm4,'ko',xnm6,ynm6,'k*','LineWidth',2,'MarkerSize',15)
legend('MC','UT','CKF','GH3','GH4','CUT4','CUT6')
% axis square
xlabel('x')
ylabel('y')
%  axis([-0.35 0.35 0.9 1.04])
plot_prop_paper
% 
% 
% plot(xmc,ymc,'kd',xut,yut,'ks','LineWidth',2,'MarkerSize',10)
% legend('MC','UT')
% xlabel('x')
% ylabel('y')
% set(gca,'FontSize',16)
% % title('Cartesian (cm,cm), mur=50, muth=0, sigr=0.02, sigth=30')
% h = get(gca, 'title');
% k = get(gca, 'xlabel');
% l = get(gca, 'ylabel');
% set(h, 'FontName', 'Helvetica', 'FontSize', 16)
% set(k, 'FontName', 'Helvetica', 'FontSize', 16)
% set(l, 'FontName', 'Helvetica', 'FontSize', 16)
% 
% %figure(3)
% %plot(p_XY(1:10^4,1),p_XY(1:10^4,2),'b+','LineWidth',2)
% %xlabel('x')
% %ylabel('y')
