% simulating the example in CKF paper

 function reentry_lol_paper_eg_2()
% matlabpool open 3
warning off
clc
clear
global dt
t0=0;
dt=50*10^(-3);
tF=20;
    % Dynamics and Measurement fns required by system and filters
    %dynamics are continuous
    %Measuremnt are discrete
    
    Dyn_disc=@reEntryDynamics;
    Meas_disc=@rangeAndBearing;
%     Dyn_disc_Jac=@UKF_eg_dyn_disc_Jac;
%     Meas_disc_Jac=@UKF_eg_meas_disc_Jac;
    
    % MC truth conditions
    x0tr=[10*1e3,100*1e3,90*1e3,-900,-1750,-400]';
%     P0tr=diag([1e-6,1e-6,1e-6,1e-6,0]);
%  P0tr=diag([200,10,200,10,1e-2]);
    n=length(x0tr);
    
Qt=diag([dt*1,dt*1,dt*1,dt*(10000),dt*(10000),dt*(10000)]);

% Qt=diag([0,0,2.4064*1e-5,2.4064*1e-5,0]);   

% R=diag([(1*1e-3)^2,(17*1e-3)^2]);
% paraset=0;
%     % frq of update : every f steps
% for f=80:15:82
%     for rr=2:1:2
% paraset=paraset+1;
f=1;
% R=diag([(1*1e-3)^2,(1*17*1e-3)^2]);
% 
%     x_0f=[6502.4,349.14,-1.2093,-6.1967,0]';
%     P_0f=diag([3e-1,1e-6,1e-6,1e-6,1.2]);
%     
% lol_eg_para_set.freq=f;
% lol_eg_para_set.R=R;
% lol_eg_para_set.x0tr=x0tr;
% lol_eg_para_set.P0tr=P0tr;
% lol_eg_para_set.Qf=Qf;
% lol_eg_para_set.Qt=Qt;
% lol_eg_para_set.x_0f=x_0f;
% lol_eg_para_set.P_0f=P_0f;

    
%% 

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %monte carlo truth generation using cont dynamics       
%     [t,x_mc]=ode45_discc(Dyn_disc,t0,dt,tF,x0tr,1e-200);
%     %generating the measurement
%     for i=1:1:length(t)
%         ym(i,:)=Meas_disc(x_mc(i,:))+(sqrtm(R)*randn(2,1));
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure(1)
% plot(xc(:,1),xc(:,3),ym(:,1).*cos(ym(:,2)),ym(:,1).*sin(ym(:,2)),'k*')
% legend('MC','sensor')
NNN=1;
for j=1:1:NNN
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %monte carlo truth generation using cont dynamics    
%     x0trr = mvnrnd(x0tr', P0tr);
    [t,X]=ode45_discc(Dyn_disc,t0,dt,tF,x0tr,Qt);
    %generating the measurement
%     for i=1:1:length(t)
%         ym(i,:)=Meas_disc(X(i,:))+(sqrtm(R)*randn(2,1));
%     end
    x_mc(:,:,j)=X;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
%     ,ym(:,1).*cos(ym(:,2)),ym(:,1).*sin(ym(:,2)),'k*'
% for i=1:1:length(t)
plot3(X(:,1),X(:,2),X(:,3),0,0,0,'ro')
hold on
% pause(0.01)
% % legend('MC','sensor')
% end

d


    
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
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %Filtering
    for tt=t0:dt:tF-dt
        disp(['freq = ',num2str(f),' rr = ',num2str(rr),' Time step no. ',num2str(tt),' of MC Iteration no. ',num2str(j)]);
            
     pp=find(abs(t-(tt+dt))<=1e-10);
        if rem(pp,f)==0
            zm=ym(pp,:)';
        else
            zm=-1234;
        end

      [mu_ut,P_ut]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_ut,P_ut,zm,Qf,R,'ut',1);
      [mu_ckf,P_ckf]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_ckf,P_ckf,zm,Qf,R,'ckf',0);
      [mu_cut4,P_cut4]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_cut4,P_cut4,zm,Qf,R,'cut4',0);
      [mu_cut6,P_cut6]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_cut6,P_cut6,zm,Qf,R,'cut6',0);
      [mu_cut8,P_cut8]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_cut8,P_cut8,zm,Qf,R,'cut8',0);
      [mu_gh,P_gh]=QUADpts_disc_disc_modified(Dyn_disc,Meas_disc,mu_gh,P_gh,zm,Qf,R,'gh',2); 

   %    [xu,Pu]=EKF_disc(Dyn_disc,Dyn_disc_Jac,Meas_disc,Meas_disc_Jac,xu,Pu,zm,Q,R);  
        
    x_ukf=vertcat(x_ukf,mu_ut');
    Cov_ukf=vertcat(Cov_ukf,reshape(P_ut,1,n*n));
%       x_ukf=vertcat(x_ukf,mu_ut);
%     Cov_ukf=vertcat(Cov_ukf,P_ut);
    
    x_ckf=vertcat(x_ckf,mu_ckf');
    Cov_ckf=vertcat(Cov_ckf,reshape(P_ckf,1,n*n));
%     x_ckf=vertcat(x_ckf,mu_ckf);
%     Cov_ckf=vertcat(Cov_ckf,P_ckf);
     
% 
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
%     
    end
    
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
        for c=1:1:5
           est_fin_ukf(r,c)= sqrt(mean((x100_ukf(r,c,:)-x_mc(r,c,:)).^2))^2;
           est_fin_ckf(r,c)= sqrt(mean((x100_ckf(r,c,:)-x_mc(r,c,:)).^2))^2;
           est_fin_cut4(r,c)= sqrt(mean((x100_cut4(r,c,:)-x_mc(r,c,:)).^2))^2;
           est_fin_cut6(r,c)= sqrt(mean((x100_cut6(r,c,:)-x_mc(r,c,:)).^2))^2;
           est_fin_cut8(r,c)= sqrt(mean((x100_cut8(r,c,:)-x_mc(r,c,:)).^2))^2;
           est_fin_gh(r,c)= sqrt(mean((x100_gh(r,c,:)-x_mc(r,c,:)).^2))^2;
        end
        for c=1:1:n*n
             Pest_fin_ukf(r,c)= sqrt(mean(P100_ukf(r,c,:)));
             Pest_fin_ckf(r,c)= sqrt(mean(P100_ckf(r,c,:)));
             Pest_fin_cut4(r,c)= sqrt(mean(P100_cut4(r,c,:)));
             Pest_fin_cut6(r,c)= sqrt(mean(P100_cut6(r,c,:)));
             Pest_fin_cut8(r,c)= sqrt(mean(P100_cut8(r,c,:)));
             Pest_fin_gh(r,c)= sqrt(mean(P100_gh(r,c,:)));
        end
   end

save(strcat('ukf_set2_eg_run_nos_',num2str(NNN),'_paraset_no_',num2str(paraset)),'ukf_eg_para_set','t','Pest_fin_ukf','Pest_fin_ckf','Pest_fin_cut4','Pest_fin_cut6','Pest_fin_cut8','Pest_fin_gh','est_fin_ukf','est_fin_ckf','est_fin_cut4','est_fin_cut6','est_fin_cut8','est_fin_gh','x_mc','x100_ukf','x100_ckf','x100_cut4','x100_cut6','x100_cut8','x100_gh','P100_ukf','P100_ckf','P100_cut4','P100_cut6','P100_cut8','P100_gh')


figure(1)
subplot(2,1,1)
plot(t,est_fin_ckf(:,1),t,est_fin_ukf(:,1),t,est_fin_cut4(:,1),t,est_fin_cut6(:,1),t,est_fin_cut8(:,1),t,est_fin_gh(:,1))
legend('ckf','ukf','cut4','cut6','cut8','gh')
subplot(2,1,2)
plot(t,Pest_fin_ckf(:,1),t,Pest_fin_ukf(:,1),t,Pest_fin_cut4(:,1),t,Pest_fin_cut6(:,1),t,Pest_fin_cut8(:,1),t,Pest_fin_gh(:,1))
legend('ckf','ukf','cut4','cut6','cut8','gh')

figure(2)
subplot(2,1,1)
plot(t,est_fin_ckf(:,2),t,est_fin_ukf(:,2),t,est_fin_cut4(:,2),t,est_fin_cut6(:,2),t,est_fin_cut8(:,2),t,est_fin_gh(:,2))
legend('ckf','ukf','cut4','cut6','cut8','gh')
subplot(2,1,2)
plot(t,Pest_fin_ckf(:,7),t,Pest_fin_ukf(:,7),t,Pest_fin_cut4(:,7),t,Pest_fin_cut6(:,7),t,Pest_fin_cut8(:,7),t,Pest_fin_gh(:,7))
legend('ckf','ukf','cut4','cut6','cut8','gh')

figure(3)
subplot(2,1,1)
plot(t,est_fin_ckf(:,3),t,est_fin_ukf(:,3),t,est_fin_cut4(:,3),t,est_fin_cut6(:,3),t,est_fin_cut8(:,3),t,est_fin_gh(:,3))
legend('ckf','ukf','cut4','cut6','cut8','gh')
subplot(2,1,2)
plot(t,Pest_fin_ckf(:,13),t,Pest_fin_ukf(:,13),t,Pest_fin_cut4(:,13),t,Pest_fin_cut6(:,13),t,Pest_fin_cut8(:,13),t,Pest_fin_gh(:,13))
legend('ckf','ukf','cut4','cut6','cut8','gh')

figure(4)
subplot(2,1,1)
plot(t,est_fin_ckf(:,4),t,est_fin_ukf(:,4),t,est_fin_cut4(:,4),t,est_fin_cut6(:,4),t,est_fin_cut8(:,4),t,est_fin_gh(:,4))
legend('ckf','ukf','cut4','cut6','cut8','gh')
subplot(2,1,2)
plot(t,Pest_fin_ckf(:,19),t,Pest_fin_ukf(:,19),t,Pest_fin_cut4(:,19),t,Pest_fin_cut6(:,19),t,Pest_fin_cut8(:,19),t,Pest_fin_gh(:,19))
legend('ckf','ukf','cut4','cut6','cut8','gh')

figure(5)
subplot(2,1,1)
plot(t,est_fin_ckf(:,5),t,est_fin_ukf(:,5),t,est_fin_cut4(:,5),t,est_fin_cut6(:,5),t,est_fin_cut8(:,5),t,est_fin_gh(:,5))
legend('ckf','ukf','cut4','cut6','cut8','gh')
subplot(2,1,2)
plot(t,Pest_fin_ckf(:,25),t,Pest_fin_ukf(:,25),t,Pest_fin_cut4(:,25),t,Pest_fin_cut6(:,25),t,Pest_fin_cut8(:,25),t,Pest_fin_gh(:,25))
legend('ckf','ukf','cut4','cut6','cut8','gh')

figure(6)
plot(x100_ckf(:,1,1),x100_ckf(:,2,1),x100_ukf(:,1,1),x100_ukf(:,2,1),x100_cut4(:,1,1),x100_cut4(:,2,1),x100_cut6(:,1,1),x100_cut6(:,2,1),x100_cut8(:,1,1),x100_cut8(:,2,1),x100_gh(:,1,1),x100_gh(:,2,1))
legend('ckf','ukf','cut4','cut6','cut8','gh')
%  matlabpool close
%  
%  

 end
