%% simulating the helicopter hovering problem

%% True plant parameters
% x=[r,v,gam,p0,Bc,g,LD]
mu0 =[3451,2.4,-9*pi/180];
paras={'p0','Bc','g','LD'};

bnds=[3445,3457;
      2.16,2.64;
      -0.1745,-0.1396];
bnds_lb=bnds(:,1);
bnds_ub=bnds(:,2);
n=length(bnds_lb);

Tt=0:1:50;

%% Monti-carlo simulations of the moments
% N=100000;
% Tt=0:0.1:15;
% Y=zeros(length(Tt),4,N);
% X=repmat(pl,N,1) + repmat((pu-pl),N,1).*rand(N,4);
% for i=1:1:N
%     p=X(i,:);
%     [t,x]=ode45(@(t,x)heli_model(t,x,[p,p0(5),p0(6)]),Tt,x0);
%     Y(:,:,i)=x;
% end
% w=(1/N)*ones(N,1);
% M1=Evol_moments_samples(Y,w,1,'central');
% M2=Evol_moments_samples(Y,w,2,'central');
% M3=Evol_moments_samples(Y,w,3,'central');
% save('MCruns','M1','M2','M3')


%% Gauss legendre simulations of the moments
Ngl_tru=6;

[Xgl,wgl]=GLeg_pts(Ngl_tru*ones(1,n),bnds_lb,bnds_ub);
Ygl=zeros(length(Tt),n,length(wgl));
for i=1:1:length(wgl)
    xx=Xgl(i,:);
    [t,y]=ode45(@(t,x)reentryvehicle(t,x,paras),Tt,xx);
  i
    Ygl(:,:,i)=y;

end

[y1,M1gl_tru6]=Evol_moments_samples(Ygl,wgl,1,'central');
[y2,M2gl_tru6]=Evol_moments_samples(Ygl,wgl,2,'central');
[y3,M3gl_tru6]=Evol_moments_samples(Ygl,wgl,3,'central');

%% Gauss legendre simulations of the moments
Ngl=5;
Tt=0:0.1:15;
[Xgl,wgl]=GLeg_pts(Ngl*ones(1,4),bnds_lb,bnds_ub);
Ygl=zeros(length(Tt),4,length(wgl));
for i=1:1:length(wgl)
    xx=Xgl(i,:);
    [t,y]=ode45(@(t,x)heli_model(t,x,[p,p0(5),p0(6)]),Tt,xx);
  
    Ygl(:,:,i)=y;

end

[y1,M1gl]=Evol_moments_samples(Ygl,wgl,1,'central');
[y2,M2gl]=Evol_moments_samples(Ygl,wgl,2,'central');
[y3,M3gl]=Evol_moments_samples(Ygl,wgl,3,'central');
%% Gauss legendre simulations of the moments
Ncut=8;
Tt=0:0.1:15;
[Xcut,wcut]=uniform_sigma_pts(bnds_lb,bnds_ub,Ncut);

Ycut=zeros(length(Tt),4,length(wcut));
for i=1:1:length(wcut)
    xx=Xcut(i,:);
    [t,y]=ode45(@(t,x)heli_model(t,x,[p,p0(5),p0(6)]),Tt,xx);
  
    Ycut(:,:,i)=y;

end

[y,M1cut]=Evol_moments_samples(Ycut,wcut,1,'central');
[y,M2cut]=Evol_moments_samples(Ycut,wcut,2,'central');
[y,M3cut]=Evol_moments_samples(Ycut,wcut,3,'central');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
MCnorms=zeros(length(t),3);
Mglnorms=zeros(length(t),3);
Mcutnorms=zeros(length(t),3);
for i=1:1:length(t)
MCnorms(i,:)=[sqrt(sum(M1(i,:).^2)),sqrt(sum(M2(i,:).^2)),sqrt(sum(M3(i,:).^2))];
Mglnorms(i,:)=[sqrt(sum(M1gl(i,:).^2)),sqrt(sum(M2gl(i,:).^2)),sqrt(sum(M3gl(i,:).^2))];
Mcutnorms(i,:)=[sqrt(sum(M1cut(i,:).^2)),sqrt(sum(M2cut(i,:).^2)),sqrt(sum(M3cut(i,:).^2))];
end

figure(4)
plot(t,MCnorms(:,1),'k.-',t,Mglnorms(:,1),'ro-',t,Mcutnorms(:,1),'b^-','linewidth',2)
legend('MC','GL3','CUT4')
title('Mean-2norm')
plot_prop_paper

figure(5)
plot(t,MCnorms(:,2),'k.-',t,Mglnorms(:,2),'ro-',t,Mcutnorms(:,2),'b^-','linewidth',2)
legend('MC','GL3','CUT4')
title('2^{nd} Moment-2norm')
plot_prop_paper

figure(6)
plot(t,MCnorms(:,3),'k.-',t,Mglnorms(:,3),'ro-',t,Mcutnorms(:,3),'b^-','linewidth',2)
legend('MC','GL3','CUT4')
title('3^{rd} Mom-2norm')
plot_prop_paper
% keyboard
%%%%%%%%%%%%%%%%%%%%%



%with respect to gh7
[sqrt(sum(sqrt(sum((M1gl-M1gl_tru).^2,2)/size(y1,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M1cut-M1gl_tru).^2,2)/size(y1,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M1-M1gl_tru).^2,2)/size(y1,1)).^2)/length(t));
 sqrt(sum(sqrt(sum((M2gl-M2gl_tru).^2,2)/size(y2,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M2cut-M2gl_tru).^2,2)/size(y2,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M2-M2gl_tru).^2,2)/size(y2,1)).^2)/length(t));
 sqrt(sum(sqrt(sum((M3gl-M3gl_tru).^2,2)/size(y3,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M3cut-M3gl_tru).^2,2)/size(y3,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M3-M3gl_tru).^2,2)/size(y3,1)).^2)/length(t))]

%with respect to gh8
[sqrt(sum(sqrt(sum((M1gl-M1gl_tru8).^2,2)/size(y1,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M1cut-M1gl_tru8).^2,2)/size(y1,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M1-M1gl_tru8).^2,2)/size(y1,1)).^2)/length(t));
 sqrt(sum(sqrt(sum((M2gl-M2gl_tru8).^2,2)/size(y2,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M2cut-M2gl_tru8).^2,2)/size(y2,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M2-M2gl_tru8).^2,2)/size(y2,1)).^2)/length(t));
 sqrt(sum(sqrt(sum((M3gl-M3gl_tru8).^2,2)/size(y3,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M3cut-M3gl_tru8).^2,2)/size(y3,1)).^2)/length(t)),sqrt(sum(sqrt(sum((M3-M3gl_tru8).^2,2)/size(y3,1)).^2)/length(t))]


%gl-8
[sqrt(sum(sqrt(sum((M1gl-M1gl_tru8).^2,2)/size(y1,1)).^2)/length(t));
 sqrt(sum(sqrt(sum((M2gl-M2gl_tru8).^2,2)/size(y2,1)).^2)/length(t));
 sqrt(sum(sqrt(sum((M3gl-M3gl_tru8).^2,2)/size(y3,1)).^2)/length(t))]

%  6-8
[sqrt(sum(sqrt(sum((M1gl_tru6-M1gl_tru8).^2,2)/size(y1,1)).^2)/length(t));
 sqrt(sum(sqrt(sum((M2gl_tru6-M2gl_tru8).^2,2)/size(y2,1)).^2)/length(t));
 sqrt(sum(sqrt(sum((M3gl_tru6-M3gl_tru8).^2,2)/size(y3,1)).^2)/length(t))]
 

% str=strcat('CUT',num2str(Ncut));
% 
% figure(1)
% subplot(2,2,1)
% plot(Tt,M1(:,1),Tt,M1gl(:,1),Tt,M1cut(:,1))
% legend('MC','GL',str)
% subplot(2,2,2)
% plot(Tt,M1(:,2),Tt,M1gl(:,2),Tt,M1cut(:,2))
% legend('MC','GL',str)
% subplot(2,2,3)
% plot(Tt,M1(:,3),Tt,M1gl(:,3),Tt,M1cut(:,3))
% legend('MC','GL',str)
% subplot(2,2,4)
% plot(Tt,M1(:,4),Tt,M1gl(:,4),Tt,M1cut(:,4))
% legend('MC','GL',str)
% 
% figure(2)
% subplot(2,2,1)
% plot(Tt,M2(:,1),Tt,M2gl(:,1),Tt,M2cut(:,1))
% legend('MC','GL',str)
% subplot(2,2,2)
% plot(Tt,M2(:,2),Tt,M2gl(:,2),Tt,M2cut(:,2))
% legend('MC','GL',str)
% subplot(2,2,3)
% plot(Tt,M2(:,3),Tt,M2gl(:,3),Tt,M2cut(:,3))
% legend('MC','GL',str)
% subplot(2,2,4)
% plot(Tt,M2(:,4),Tt,M2gl(:,4),Tt,M2cut(:,4))
% legend('MC','GL',str)
% 
% figure(3)
% subplot(2,2,1)
% plot(Tt,M3(:,1),Tt,M3gl(:,1),Tt,M3cut(:,1))
% legend('MC','GL',str)
% subplot(2,2,2)
% plot(Tt,M3(:,2),Tt,M3gl(:,2),Tt,M3cut(:,2))
% legend('MC','GL',str)
% subplot(2,2,3)
% plot(Tt,M3(:,3),Tt,M3gl(:,3),Tt,M3cut(:,3))
% legend('MC','GL',str)
% subplot(2,2,4)
% plot(Tt,M3(:,4),Tt,M3gl(:,4),Tt,M3cut(:,4))
% legend('MC','GL',str)


figure
hold on
for i=1:1:length(wgl)
    plot(t,Ygl(:,2,i))
end