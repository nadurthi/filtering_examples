% simulating the example in CKF paper

 function CKF_paper_eg()
%  matlabpool open 5
clc
clear
global T
t0=0;
% dt=0.1;
tF=50;
    % Dynamics and Measurement fns required by system and filters
    %dynamics are continuous
    %Measuremnt are discrete
    
    Dyn_disc=@CKF_eg_dyn_disc;
    Meas_disc=@CKF_eg_meas_disc;
    Dyn_disc_Jac=@CKF_eg_dyn_disc_Jac;
    Meas_disc_Jac=@CKF_eg_meas_disc_Jac;
    para_set=0;
    for f=100000:100000:100000
           for rr=1:5:2
para_set=para_set+1;
    T=1;
    dt=T;
    omg=-3*pi/180;
    q1=0.1;
    q2=1.75*10^(-4);
    sigr=10;
    sigth=0.1833*(pi/180)*rr;
%    sigr=10;
% sigth=0.3*pi/180;
    M=[T^3/3,T^2/2;T^2/2,T];
    Q=blkdiag(q1*M,q1*M,q2*T);
    
    

    R=diag([sigr^2,sigth^2]);
    
    % frq of update : every f steps
%     f=1;
    
       % MC truth conditions
    x0tr=[1000,300,1000,0,omg]';
%     P0tr=diag([100,10,100,10,200*10^(-4)]);
  P0tr=diag([100,10,100,10,100e-6]);

  CKFeg_simulation_para.x0tr=x0tr;
  CKFeg_simulation_para.P0tr=P0tr;
  CKFeg_simulation_para.freq=f;
  CKFeg_simulation_para.R=R;
  CKFeg_simulation_para.Q=Q;
  CKFeg_simulation_para.dt=dt;
  n=length(x0tr);
    
%% 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %monte carlo truth generation using cont dynamics       
    [t,x_mc]=ode45_discc(Dyn_disc,t0,dt,tF,x0tr,1e-200);
    %generating the measurement
    for i=1:1:length(t)
        ym(i,:)=Meas_disc(x_mc(i,:))+(sqrtm(R)*randn(2,1));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure(1)
% plot(x_mc(:,1),x_mc(:,3),ym(:,1).*cos(ym(:,2)),ym(:,1).*sin(ym(:,2)),'k*')
% legend('MC','sensor')
NNN=1;
for j=1:1:NNN
    
    


    x_0f=mvnrnd(x0tr,P0tr)';

    P_0f=P0tr;
    
    %filter initial conditions
     % EKF filter
%     xu=x_0f;
%     Pu=P_0f;
%     x_ekf=xu';
%     Cov_ekf=reshape(Pu,1,n^2);
    
    % UKF filter
    mu_ut=x_0f;
    P_ut=P_0f;
    x_ukf=mu_ut';
    Cov_ukf=reshape(P_ut,1,n^2);
    
      % CKF filter
    mu_ckf=x_0f;
    P_ckf=P_0f;
    x_ckf=mu_ckf';
    Cov_ckf=reshape(P_ckf,1,n^2);
 
      % CUT4 filter
    mu_cut4=x_0f;
    P_cut4=P_0f;
    x_cut4=mu_cut4';
    Cov_cut4=reshape(P_cut4,1,n^2);
    
          % CUT6 filter
    mu_cut6=x_0f;
    P_cut6=P_0f;
    x_cut6=mu_cut6';
    Cov_cut6=reshape(P_cut6,1,n^2);
    
          % CUT8 filter
    mu_cut8=x_0f;
    P_cut8=P_0f;
    x_cut8=mu_cut8';
    Cov_cut8=reshape(P_cut8,1,n^2);
    
       % ghKF filter
    mu_gh=x_0f;
    P_gh=P_0f;
    x_gh=mu_gh';
    Cov_gh=reshape(P_gh,1,n^2);
    tic
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %Filtering
    for tt=t0:dt:tF-dt
                disp(['frq=',num2str(f),'  sig_th=',num2str(rr),'  CKF eg: Time step no. ',num2str(tt),' of MC Iteration no. ',num2str(j)]);
            
     pp=find(abs(t-(tt+dt))<=1e-10);
        if rem(pp,f)==0
            zm=ym(pp,:)';
        else
            zm=-1234
        end

      [mu_ut,P_ut]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_ut,P_ut,zm,Q,R,'ut',1);
      [mu_ckf,P_ckf]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_ckf,P_ckf,zm,Q,R,'ckf',0);
      [mu_cut4,P_cut4]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_cut4,P_cut4,zm,Q,R,'cut4',0);
      [mu_cut6,P_cut6]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_cut6,P_cut6,zm,Q,R,'cut6',0);
      [mu_cut8,P_cut8]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_cut8,P_cut8,zm,Q,R,'cut8',0);
      [mu_gh,P_gh]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_gh,P_gh,zm,Q,R,'gh',7); 

   %    [xu,Pu]=EKF_disc(Dyn_disc,Dyn_disc_Jac,Meas_disc,Meas_disc_Jac,xu,Pu,zm,Q,R);  
        
    x_ukf=vertcat(x_ukf,mu_ut');
    Cov_ukf=vertcat(Cov_ukf,reshape(P_ut,1,n*n));
%       x_ukf=vertcat(x_ukf,mu_ut);
%     Cov_ukf=vertcat(Cov_ukf,P_ut);
    
    x_ckf=vertcat(x_ckf,mu_ckf');
    Cov_ckf=vertcat(Cov_ckf,reshape(P_ckf,1,n*n));
%     x_ckf=vertcat(x_ckf,mu_ckf);
%     Cov_ckf=vertcat(Cov_ckf,P_ckf);
     

    x_cut4=vertcat(x_cut4,mu_cut4');
    Cov_cut4=vertcat(Cov_cut4,reshape(P_cut4,1,n*n));
%     x_cut4=vertcat(x_cut4,mu_cut4);
%     Cov_cut4=vertcat(Cov_cut4,P_cut4);
    
    x_cut6=vertcat(x_cut6,mu_cut6');
    Cov_cut6=vertcat(Cov_cut6,reshape(P_cut6,1,n*n));
%     x_cut6=vertcat(x_cut6,mu_cut6);
%     Cov_cut6=vertcat(Cov_cut6,P_cut6);
    
    x_cut8=vertcat(x_cut8,mu_cut8');
    Cov_cut8=vertcat(Cov_cut8,reshape(P_cut8,1,n*n));
%     x_cut8=vertcat(x_cut8,mu_cut8);
%     Cov_cut8=vertcat(Cov_cut8,P_cut8);

    x_gh=vertcat(x_gh,mu_gh');
    Cov_gh=vertcat(Cov_gh,reshape(P_gh,1,n*n));
%     x_gh=vertcat(x_gh,mu_gh);
%     Cov_gh=vertcat(Cov_gh,P_gh);    

%     x_ekf=vertcat(x_ekf,xu');
%     Cov_ekf=vertcat(Cov_ekf,reshape(Pu,1,n^2));
%     mu_ut=mu_ut(end,:)';
%     mu_ckf=mu_ckf(end,:)';
%     mu_cut4=mu_cut4(end,:)';
%     mu_cut6=mu_cut6(end,:)';
%     mu_cut8=mu_cut8(end,:)';
%     mu_gh=mu_gh(end,:)';
%     P_ut=reshape(P_ut(end,:),n,n);
%     P_ckf=reshape(P_ckf(end,:),n,n);
%     P_cut4=reshape(P_cut4(end,:),n,n);
%     P_cut6=reshape(P_cut6(end,:),n,n);
%     P_cut8=reshape(P_cut8(end,:),n,n);
%     P_gh=reshape(P_gh(end,:),n,n);
    
    end
    toc
%      save(strcat('ckf_eg_run_no_',num2str(j)),'x_mc','x_ukf','x_ckf','x_cut4','x_cut6','x_cut8','x_gh','ym','Cov_ukf','Cov_fin_ckf','Cov_fin_cut4','Cov_fin_cut6','Cov_fin_cut8','Cov_fin_gh')
%     x100_mc(:,:,j)=xc;
% x100_ekf(:,:,j)=x_ekf;
x100_ukf(:,:,j)=x_ukf;
x100_ckf(:,:,j)=x_ckf;
x100_cut4(:,:,j)=x_cut4;
x100_cut6(:,:,j)=x_cut6;
x100_cut8(:,:,j)=x_cut8;
x100_gh(:,:,j)=x_gh;

% P100_ekf(:,:,j)=Cov_ekf;
P100_ukf(:,:,j)=Cov_ukf;
P100_ckf(:,:,j)=Cov_ckf;
P100_cut4(:,:,j)=Cov_cut4;
P100_cut6(:,:,j)=Cov_cut6;
P100_cut8(:,:,j)=Cov_cut8;
P100_gh(:,:,j)=Cov_gh;
end


%%%%%%%%%%%%%%%%%---Averaging the distance estimate over runs---%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for r=1:1:length(t)
        for c=0:1:2
%            x_fin_mc(r,c)= mean(x100_mc(r,c,:));
%            x_fin_ekf(r,c)=mean(x100_ekf(r,c,:));
           if c<=1
           est_fin_ukf(r,c+1)= sqrt(mean((x100_ukf(r,c+1,:)-x_mc(r,c+1)).^2+(x100_ukf(r,c+3,:)-x_mc(r,c+3)).^2));
           est_fin_ckf(r,c+1) =sqrt(mean((x100_ckf(r,c+1,:)-x_mc(r,c+1)).^2+(x100_ckf(r,c+3,:)-x_mc(r,c+3)).^2));
           est_fin_cut4(r,c+1)=sqrt(mean((x100_cut4(r,c+1,:)-x_mc(r,c+1)).^2+(x100_cut4(r,c+3,:)-x_mc(r,c+3)).^2));
           est_fin_cut6(r,c+1)=sqrt(mean((x100_cut6(r,c+1,:)-x_mc(r,c+1)).^2+(x100_cut6(r,c+3,:)-x_mc(r,c+3)).^2));
           est_fin_cut8(r,c+1)=sqrt(mean((x100_cut8(r,c+1,:)-x_mc(r,c+1)).^2+(x100_cut8(r,c+3,:)-x_mc(r,c+3)).^2));
           est_fin_gh(r,c+1)=  sqrt(mean((x100_gh(r,c+1,:)-x_mc(r,c+1)).^2+(x100_gh(r,c+3,:)-x_mc(r,c+3)).^2));
           else
           est_fin_ukf(r,c+1)= sqrt(mean((x100_ukf(r,5,:)-x_mc(r,5)).^2));
           est_fin_ckf(r,c+1) =sqrt(mean((x100_ckf(r,5,:)-x_mc(r,5)).^2));
           est_fin_cut4(r,c+1)=sqrt(mean((x100_cut4(r,5,:)-x_mc(r,5)).^2));
           est_fin_cut6(r,c+1)=sqrt(mean((x100_cut6(r,5,:)-x_mc(r,5)).^2));
           est_fin_cut8(r,c+1)=sqrt(mean((x100_cut8(r,5,:)-x_mc(r,5)).^2));
           est_fin_gh(r,c+1)=  sqrt(mean((x100_gh(r,5,:)-x_mc(r,5)).^2)); 
           end
        end
   end
save(strcat('ckf_eg_run_nos_',num2str(NNN),'paraset_',num2str(para_set)),'CKFeg_simulation_para','t','est_fin_ukf','est_fin_ckf','est_fin_cut4','est_fin_cut6','est_fin_cut8','est_fin_gh','x_mc','x100_ukf','x100_ckf','x100_cut4','x100_cut6','x100_cut8','x100_gh','ym','P100_ukf','P100_ckf','P100_cut4','P100_cut6','P100_cut8','P100_gh')
           end
    end
figure(1)
plot(t,est_fin_ckf(:,1),t,est_fin_ukf(:,1),t,est_fin_cut4(:,1),t,est_fin_cut6(:,1),t,est_fin_cut8(:,1),t,est_fin_gh(:,1))
legend('ckf','ukf','cut4','cut6','cut8','gh')
figure(2)
plot(t,est_fin_ckf(:,2),t,est_fin_ukf(:,2),t,est_fin_cut4(:,2),t,est_fin_cut6(:,2),t,est_fin_cut8(:,2),t,est_fin_gh(:,2))
legend('ckf','ukf','cut4','cut6','cut8','gh')
figure(3)
plot(t,(180/pi)*est_fin_ckf(:,3),t,(180/pi)*est_fin_ukf(:,3),t,(180/pi)*est_fin_cut4(:,3),t,(180/pi)*est_fin_cut6(:,3),t,(180/pi)*est_fin_cut8(:,3),t,(180/pi)*est_fin_gh(:,3))
legend('ckf','ukf','cut4','cut6','cut8','gh')
figure(4)
plot(x100_ckf(:,1,end),x100_ckf(:,3,end),x100_ukf(:,1,end),x100_ukf(:,3,end),x100_cut4(:,1,end),x100_cut4(:,3,end),x100_cut6(:,1,end),x100_cut6(:,3,end),x100_cut8(:,1,end),x100_cut8(:,3,end),x100_gh(:,1,end),x100_gh(:,3,end),x_mc(:,1),x_mc(:,3),ym(:,1).*cos(ym(:,2)),ym(:,1).*sin(ym(:,2)),'k*','LineWidth',2)
legend('ckf','ukf','cut4','cut6','cut8','gh','MC','sensor')
%    for r=1:1:length(t)
%        for c=1:1:n^2
% %            P_fin_ekf(r,c)=mean(P100_ekf(r,c,:));
%            P_fin_ukf(r,c)=mean(P100_ukf(r,c,:));
%            P_fin_ckf(r,c)=mean(P100_ckf(r,c,:));
%            P_fin_cut4(r,c)=mean(P100_cut4(r,c,:));
%            P_fin_cut6(r,c)=mean(P100_cut6(r,c,:));
%            P_fin_cut8(r,c)=mean(P100_cut8(r,c,:));
%            P_fin_gh(r,c)=mean(P100_gh(r,c,:));
%        end
%    end

% save(strcat('ckf_eg_run_no_',num2str(j)),'x_fin_mc','x_fin_ukf','x_fin_ckf','x_fin_cut4','x_fin_cut6','x_fin_cut8','x_fin_gh','ym','P_fin_ukf','P_fin_ckf','P_fin_cut4','P_fin_cut6','P_fin_cut8','P_fin_gh')
%  save(strcat('Cov_ckf_eg_run_no_',num2str(ppmax)),x_fin_mc,x_fin_ukf,x_fin_ckf,x_fin_cut4,x_fin_cut6,x_fin_cut8,x_fin_gh,ym)
%    save x_fin_mc
% % save x_fin_ekf
% save x_fin_ukf
% save x_fin_ckf
% save x_fin_cut4
% save x_fin_cut6
% save x_fin_cut8
% save x_fin_gh
% 
% % save P_fin_ekf
% save P_fin_ukf
% save P_fin_ckf
% save P_fin_cut4
% save P_fin_cut6
% save P_fin_cut8
% save P_fin_gh


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load x_fin_mc
% % load x_fin_ekf
% load x_fin_ukf
% load x_fin_ckf
% load x_fin_cut4
% load x_fin_cut6
% load x_fin_cut8
% load x_fin_gh
% 
% % plot(t,sqrt((x_fin_mc(:,1)-x_fin_ukf(:,1)).^2+(x_fin_mc(:,3)-x_fin_ukf(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_ekf(:,1)).^2+(x_fin_mc(:,3)-x_fin_ekf(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_ckf(:,1)).^2+(x_fin_mc(:,3)-x_fin_ckf(:,3)).^2)...
% %     ,t,sqrt((x_fin_mc(:,1)-x_fin_cut4(:,1)).^2+(x_fin_mc(:,3)-x_fin_cut4(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_gh(:,1)).^2+(x_fin_mc(:,3)-x_fin_gh(:,3)).^2))
% % legend('ukf','ekf','ckf','cut4','gh')
% % axis([70,100,10,30])
% % plot(t,sqrt((x_fin_mc(:,1)-x_fin_ckf(:,1)).^2+(x_fin_mc(:,3)-x_fin_ckf(:,3)).^2)...
% %     ,t,sqrt((x_fin_mc(:,1)-x_fin_gh(:,1)).^2+(x_fin_mc(:,3)-x_fin_gh(:,3)).^2))
% figure(2)
% subplot(2,3,1)
% plot(t,x_fin_mc(:,1),t,x_fin_ckf(:,1))
% legend('MC','ckf')
% subplot(2,3,2)
% plot(t,x_fin_mc(:,1),t,x_fin_ukf(:,1))
% legend('MC','ukf')
% subplot(2,3,3)
% plot(t,x_fin_mc(:,1),t,x_fin_cut4(:,1))
% legend('MC','cut4')
% subplot(2,3,4)
% plot(t,x_fin_mc(:,1),t,x_fin_cut6(:,1))
% legend('MC','cut6')
% subplot(2,3,5)
% plot(t,x_fin_mc(:,1),t,x_fin_cut8(:,1))
% legend('MC','cut8')
% subplot(2,3,6)
% plot(t,x_fin_mc(:,1),t,x_fin_gh(:,1))
% legend('MC','GH')
% 
% figure(3)
% subplot(2,3,1)
% plot(t,x_fin_mc(:,3),t,x_fin_ckf(:,3))
% legend('MC','ckf')
% subplot(2,3,2)
% plot(t,x_fin_mc(:,3),t,x_fin_ukf(:,3))
% legend('MC','ukf')
% subplot(2,3,3)
% plot(t,x_fin_mc(:,3),t,x_fin_cut4(:,3))
% legend('MC','cut4')
% subplot(2,3,4)
% plot(t,x_fin_mc(:,3),t,x_fin_cut6(:,3))
% legend('MC','cut6')
% subplot(2,3,5)
% plot(t,x_fin_mc(:,3),t,x_fin_cut8(:,3))
% legend('MC','cut8')
% subplot(2,3,6)
% plot(t,x_fin_mc(:,3),t,x_fin_gh(:,3))
% legend('MC','GH')
% 
% figure(4)
% plot(t,sqrt((x_fin_mc(:,1)-x_fin_ukf(:,1)).^2+(x_fin_mc(:,3)-x_fin_ukf(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_ckf(:,1)).^2+(x_fin_mc(:,3)-x_fin_ckf(:,3)).^2)...
%      ,t,sqrt((x_fin_mc(:,1)-x_fin_cut4(:,1)).^2+(x_fin_mc(:,3)-x_fin_cut4(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_cut6(:,1)).^2+(x_fin_mc(:,3)-x_fin_cut6(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_cut8(:,1)).^2+(x_fin_mc(:,3)-x_fin_cut8(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_gh(:,1)).^2+(x_fin_mc(:,3)-x_fin_gh(:,3)).^2))
%  legend('ukf','ckf','cut4','cut6','cut8','gh')
% 
%  
%  matlabpool close
%  
%  

 end
%%discrete dynamic equations for MC
function xk1=CKF_eg_dyn_disc(xk)
% global T
T=1;
omg=xk(end,1);
xk1=[1,sin(omg*T)/omg,0,-(1-cos(omg*T))/omg,0;...
     0,cos(omg*T),0,-sin(omg*T),0;
     0,(1-cos(omg*T))/omg,1,sin(omg*T)/omg,0;
     0,sin(omg*T),0,cos(omg*T),0;
     0,0,0,0,1]*xk;

end

%% discrete mesuremtn function
    function yk=CKF_eg_meas_disc(x)
    zi=x(1)-0;
    ni=x(3)+0;
    
    yk(1,1)=sqrt((zi)^2+(ni)^2);
    yk(2,1)=atan2(ni,zi);
    end

%% discrete Jacobian dynamic equations for EKF filters
function J=CKF_eg_dyn_disc_Jac(x)
global T
a=5*10^(-5);
J=[1,sin(T*x(5))/x(5),0,(-1+cos(T*x(5)))/x(5),-((x(4)*(-1 + cos(T*x(5))))/x(5)^2)+(T*x(2)*cos(T*x(5)))/x(5)-(x(2)*sin(T*x(5)))/x(5)^2-(T*x(4)*sin(T*x(5)))/x(5);...
    0, cos(T*x(5)),0,-sin(T*x(5)),-T*x(4)*cos(T*x(5))-T*x(2)*sin(T*x(5));...
   1,(1-cos(T*x(5)))/x(5),1,(sin(T*x(5)))/x(5),-((x(2)*(1 - cos(T*x(5))))/x(5)^2)+(T*x(4)*cos(T*x(5)))/x(5)-(x(4)*sin(T*x(5)))/x(5)^2+(T*x(2)*sin(T*x(5)))/x(5);...
    0,sin(T*x(5)),0,cos(T*x(5)), T*x(2)*cos(T*x(5))- T*x(4)*sin(T*x(5));...
   0,0,0,0,1];
end

%% discrete Jacobian measurement equations for filters
function H=CKF_eg_meas_disc_Jac(x)
 zi=x(1);
 ni=x(3);
H=[x(1)/sqrt(x(1)^2+x(3)^2),0,x(3)/sqrt(x(1)^2+x(3)^2),0,0;-(x(3)/x(1)^2)/(1+(x(3)/x(1))^2),0,(1/x(1))/(1+(x(3)/x(1))^2),0,0];
end