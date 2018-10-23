% simulating the example in Crassidis book
%% poly nonlinear function
 function poly_nonlin_fn()
clc
clear
global dt
t0=0;
dt=1;
tF=60;
    % Dynamics and Measurement fns required by system and filters
    %dynamics are continuous
    %Measuremnt are discrete
    
    Dyn_cont=@crass_eg_dyn_cont;
    Meas_dis=@crass_eg_meas_disc;
    Dyn_cont_Jac=@crass_eg_dyn_disc_Jac;
    Meas_dis_Jac=@crass_eg_meas_disc_Jac;
    
    % MC truth conditions
    x0tr=[3*10^5;2*10^4;1*10^(-3)];
    n=length(x0tr);
    
    % Common Filter  conditions
    x_0f=[3*10^5;2*10^4;3*10^(-5)];
    P_0f=diag([10^6,4*10^6,10^(-4)]);
    Q=0;
    R=10^4;
    % frq of update : every f steps
    f=1;
    
   

    
%% 
for j=1:1:10
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %monte carlo truth generation using cont dynamics       
    [t,xc]=ode45(@crass_eg_dyn_cont,t0:dt:tF,x0tr);
    %generating the measurement
    for i=1:1:length(t)
        ym(i,:)=crass_eg_meas_disc(xc(i,:))+sqrt(10^4)*randn;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %filter initial conditions
     % EKF filter
    xu=x_0f;
    Pu=P_0f;
    x_ekf=xu';
    Cov_ekf=reshape(Pu,1,n^2);
    
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
 
      % NM4KF filter
    mu_nm4=x_0f;
    P_nm4=P_0f;
    x_nm4=mu_nm4';
    Cov_nm4=reshape(P_nm4,1,n^2);
    
       % GHKF filter
    mu_GH=x_0f;
    P_GH=P_0f;
    x_GH=mu_GH';
    Cov_GH=reshape(P_GH,1,n^2);
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %Filtering
    for tt=t0:dt:tF-dt
            
     pp=find(abs(t-(tt+dt))<=1e-10);
        if rem(pp,f)==0
            zm=ym(pp,:)';
        else
            zm=-1;
        end
  [mu_ut,P_ut]=UKF_cont_disc(Dyn_cont,Meas_dis,mu_ut,P_ut,zm,Q,R,tt,dt,tt+dt);
  [mu_ckf,P_ckf]=CKF_cont_disc(Dyn_cont,Meas_dis,mu_ckf,P_ckf,zm,Q,R,tt,dt,tt+dt);
  [mu_nm4,P_nm4]=NM4_cont_disc(Dyn_cont,Meas_dis,mu_nm4,P_nm4,zm,Q,R,tt,dt,tt+dt);
  [mu_GH,P_GH]=GHKF_cont_disc(Dyn_cont,Meas_dis,mu_GH,P_GH,zm,Q,R,tt,dt,tt+dt,3); 
  [xu,Pu]=EKF(Dyn_cont,Dyn_cont_Jac,Meas_dis,Meas_dis_Jac,xu,Pu,zm,Q,R,tt,dt,tt+dt);  
        
    x_ukf=vertcat(x_ukf,mu_ut');
    Cov_ukf=vertcat(Cov_ukf,reshape(P_ut,1,n*n));
    
    x_ckf=vertcat(x_ckf,mu_ckf');
    Cov_ckf=vertcat(Cov_ckf,reshape(P_ckf,1,n*n));
    
    x_nm4=vertcat(x_nm4,mu_nm4');
    Cov_nm4=vertcat(Cov_nm4,reshape(P_nm4,1,n*n));
    
    x_GH=vertcat(x_GH,mu_GH');
    Cov_GH=vertcat(Cov_GH,reshape(P_GH,1,n*n));
    
    x_ekf=vertcat(x_ekf,xu');
    Cov_ekf=vertcat(Cov_ekf,reshape(Pu,1,n^2));
    
    end
    
x100_mc(:,:,j)=xc;
x100_ekf(:,:,j)=x_ekf;
x100_ukf(:,:,j)=x_ukf;
x100_ckf(:,:,j)=x_ckf;
x100_nm4(:,:,j)=x_nm4;
x100_GH(:,:,j)=x_GH;

P100_ekf(:,:,j)=Cov_ekf;
P100_ukf(:,:,j)=Cov_ukf;
P100_ckf(:,:,j)=Cov_ckf;
P100_nm4(:,:,j)=Cov_nm4;
P100_GH(:,:,j)=Cov_GH;

j
end
%%%%%%%%%%%%%%%%%---Averaging over runs---%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for r=1:1:length(t)
       for c=1:1:n
           x_fin_mc(r,c)= mean(x100_mc(r,c,:));
           x_fin_ekf(r,c)=mean(x100_ekf(r,c,:));
           x_fin_ukf(r,c)=mean(x100_ukf(r,c,:));
           x_fin_ckf(r,c)=mean(x100_ckf(r,c,:));
           x_fin_nm4(r,c)=mean(x100_nm4(r,c,:));
           x_fin_GH(r,c)=mean(x100_GH(r,c,:));
       end
   end
   
   for r=1:1:length(t)
       for c=1:1:n^2
           P_fin_ekf(r,c)=mean(P100_ekf(r,c,:));
           P_fin_ukf(r,c)=mean(P100_ukf(r,c,:));
           P_fin_ckf(r,c)=mean(P100_ckf(r,c,:));
           P_fin_nm4(r,c)=mean(P100_nm4(r,c,:));
           P_fin_GH(r,c)=mean(P100_GH(r,c,:));
       end
   end
      
save x_fin_mc
save x_fin_ekf
save x_fin_ukf
save x_fin_ckf
save x_fin_nm4
save x_fin_GH

save P_fin_ekf
save P_fin_ukf
save P_fin_ckf
save P_fin_nm4
save P_fin_GH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plot(t,abs(x_fin_mc(:,1)-x_fin_ukf(:,1)),t,abs(x_fin_mc(:,1)-x_fin_ekf(:,1)),t,abs(x_fin_mc(:,1)-x_fin_ckf(:,1))...
    ,t,abs(x_fin_mc(:,1)-x_fin_nm4(:,1)),t,abs(x_fin_mc(:,1)-x_fin_GH(:,1)))
legend('ukf','ekf','ckf','nm4','GH')



 end
%% continuous dynamic equations for MC
function dx=crass_eg_dyn_cont(t,x)
a=5*10^(-5);
dx=[-x(2);-exp(-a*x(1))*x(2)^2*x(3);0];
end

%% discrete mesuremtn function
    function yk=crass_eg_meas_disc(x)
    M=10^5;
    Z=10^5;
    yk=sqrt(M^2+(x(1)-Z)^2);
    end

%% discrete Jacobian dynamic equations for EKF filters
function J=crass_eg_dyn_disc_Jac(x)
global dt
a=5*10^(-5);
J=exp(-a*x(1))*[0,-exp(a*x(1)),0;a*x(2)^2*x(3),-2*x(2)*x(3),-x(2)^2;0,0,0];
end

%% discrete Jacobian measurement equations for filters
function H=crass_eg_meas_disc_Jac(x)
global dt
 M=10^5;
 Z=10^5;
H=[(x(1)-Z)/sqrt(M^2+(x(1)-Z)^2),0,0];
end