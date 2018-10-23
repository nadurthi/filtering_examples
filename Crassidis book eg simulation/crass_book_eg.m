% simulating the example in Crassidis book
%% free fall body tracking
function body_fall_track_system()
clc
clear
global dt
t0=0;
dt=0.05;
tF=60;
%% monte carlo truth generation using cont dynamics
for j=1:1:10
x0tr=[3*10^5;2*10^4;1*10^(-3)];
N=length(x0tr);
% for i=1:1:1
    [t,xc]=ode45(@crass_eg_dyn_cont,t0:dt:tF,x0tr);
    %generating the measurement
    for i=1:1:length(t)
        ym(i,:)=crass_eg_meas_disc(xc(i,:))+sqrt(10^4)*randn;
    end
    
    % Filter  conditions
    x_0f=[3*10^5;2*10^4;3*10^(-5)];
    P_0f=diag([10^6,4*10^6,10^(-4)]);
    Q=0;
    R=10^4;
    % frq of update : every f steps
    f=100;
    
    n=length(x_0f);
    
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
    for tt=t0:dt:tF-dt
        
        %ukf filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %state forecast step
     pp=find(abs(t-(tt+dt))<=1e-10);
        if rem(pp,f)==0
%              pp=find(abs(t-tt)<=1e-10);
            zm=ym(pp,:)';
        
        else
            zm=-1;
        end
   [mu_ut,P_ut]=UKF_cont_disc(@crass_eg_dyn_cont,@crass_eg_meas_disc,mu_ut,P_ut,zm,Q/dt,R,tt,dt,tt+dt);
  [mu_ckf,P_ckf]=CKF_cont_disc(@crass_eg_dyn_cont,@crass_eg_meas_disc,mu_ckf,P_ckf,zm,Q/dt,R,tt,dt,tt+dt);
   [mu_nm4,P_nm4]=NM4_cont_disc(@crass_eg_dyn_cont,@crass_eg_meas_disc,mu_nm4,P_nm4,zm,Q/dt,R,tt,dt,tt+dt); 
    
    x_ukf=vertcat(x_ukf,mu_ut');
    Cov_ukf=vertcat(Cov_ukf,reshape(P_ut,1,n*n));
    
    x_ckf=vertcat(x_ckf,mu_ckf');
    Cov_ckf=vertcat(Cov_ckf,reshape(P_ckf,1,n*n));
    
    x_nm4=vertcat(x_nm4,mu_nm4');
    Cov_nm4=vertcat(Cov_nm4,reshape(P_nm4,1,n*n));
 %Ekf filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     pp=find(abs(t-(tt+dt))<=1e-10);
        if rem(pp,f)==0
%              pp=find(abs(t-tt)<=1e-10);
            zm=ym(pp,:)';
        
        else
            zm=-1;
        end
        [xu,Pu]=EKF(@crass_eg_dyn_cont,@crass_eg_dyn_disc_Jac,@crass_eg_meas_disc,@crass_eg_meas_disc_Jac...
            ,xu,Pu,zm,Q,R,tt,dt,tt+dt);
        x_ekf=vertcat(x_ekf,xu');
        Cov_ekf=vertcat(Cov_ekf,reshape(Pu,1,n^2));
%         pause
    end
x100_mc(:,:,j)=xc;
x100_ekf(:,:,j)=x_ekf;
x100_ukf(:,:,j)=x_ukf;
x100_ckf(:,:,j)=x_ckf;
x100_nm4(:,:,j)=x_nm4;

P100_ekf(:,:,j)=Cov_ekf;
P100_ukf(:,:,j)=Cov_ukf;
P100_ckf(:,:,j)=Cov_ckf;
P100_nm4(:,:,j)=Cov_nm4;
j
end
save x100_mc
save x100_ekf
save x100_ukf
save x100_ckf
save x100_nm4

   for r=1:1:length(t)
       for c=1:1:N
           x_fin_mc(r,c)=mean(x100_mc(r,c,:));
           x_fin_ekf(r,c)=mean(x100_ekf(r,c,:));
           x_fin_ukf(r,c)=mean(x100_ukf(r,c,:));
           x_fin_ckf(r,c)=mean(x100_ckf(r,c,:));
           x_fin_nm4(r,c)=mean(x100_nm4(r,c,:));
       end
   end
      for r=1:1:length(t)
       for c=1:1:N^2
           
           P_fin_ekf(r,c)=mean(P100_ekf(r,c,:));
           P_fin_ukf(r,c)=mean(P100_ukf(r,c,:));
           P_fin_ckf(r,c)=mean(P100_ckf(r,c,:));
           P_fin_nm4(r,c)=mean(P100_nm4(r,c,:));
       end
      end
   figure(1)
plot(t,abs(x_fin_mc(:,1)-x_fin_ukf(:,1)),t,abs(x_fin_mc(:,1)-x_fin_ekf(:,1)),t,abs(x_fin_mc(:,1)-x_fin_ckf(:,1))...
    ,t,abs(x_fin_mc(:,1)-x_fin_nm4(:,1)))
legend('ukf','ekf','ckf','nm4')
figure(2)
plot(t,P_fin_ekf(:,1),t,P_fin_ukf(:,1),t,P_fin_ckf(:,1)...
    ,t,P_fin_nm4(:,1))
legend('ukf','ekf','ckf','nm4')














end
%% continuous dynamic equations for MC
function dx=crass_eg_dyn_cont(t,x)
a=5*10^(-5);
dx=[-x(2);-exp(-a*x(1))*x(2)^2*x(3);0];
end
%% discrete equivalent dynamic equations for filters
function xk1=crass_eg_dyn_disc(x)
global dt
a=5*10^(-5);
xk1=[x(1)-x(2)*dt;x(2)-exp(-a*x(1))*x(2)^2*x(3)*dt;x(3)+0*dt];
end

%% discrete mesuremtn function
    function yk=crass_eg_meas_disc(x)
    M=10^5;
    Z=10^5;
    yk=sqrt(M^2+(x(1)-Z)^2);
    end

%% discrete Jacobian dynamic equations for filters
function J=crass_eg_dyn_disc_Jac(x)
global dt
a=5*10^(-5);
% xk1=[x(1)-x(2)*dt;x(2)-exp(-a*x(1))*x(2)^2*x(3)*dt;x(3)+0*dt];
% J=[1,-dt,0;-exp(-a*x(1))*(-a)*x(2)^2*x(3)*dt,1-exp(-a*x(1))*2*x(2)*x(3)*dt,-exp(-a*x(1))*x(2)^2*dt;0,0,1];
J=exp(-a*x(1))*[0,-exp(a*x(1)),0;a*x(2)^2*x(3),-2*x(2)^1*x(3),-x(2)^2;0,0,0];
end

%% discrete Jacobian measurement equations for filters
function H=crass_eg_meas_disc_Jac(x)
global dt
 M=10^5;
 Z=10^5;
H=[(x(1)-Z)/sqrt(M^2+(x(1)-Z)^2),0,0];
end