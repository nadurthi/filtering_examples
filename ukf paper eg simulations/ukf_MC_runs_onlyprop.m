
function ukf_MC_runs_onlyprop()
%% ------------------------------------------------------------------------
% time
global kappa
%--------------------------------------------------------------------------
time.t0 = 0;
time.dt = 0.05;
time.tf = 200;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);
%__________________________________________________________________________

%% ------------------------------------------------------------------------
% model
T=time.dt;
    
Qf=diag([0,0,1.1104*1e-3,2.4064*1e-3,1e-4]);

Qt=diag([0,0,2.4064*1e-5,2.4064*1e-5,0]);   

R=diag([(20*1e-3)^2,(1*17*1e-3)^2]);
     
model.fn = 5;               % state space dimensionality
model.fx = @reEntryDynamics;
model.fx_jac=@KIRB_CT_eg_dyn_jac_disc;

model.hn = 2;               % measurement dimensionality
model.hx =@rangeAndBearing;
model.hx_jac=@KIRB_eg_meas_jac_disc;

model.Q = zeros(model.fn,model.fn);
model.sQ =zeros(model.fn,model.fn);
model.R = R;
model.sR=sqrtm(R);
model.Qtruth=zeros(model.fn,model.fn);
model.Qt_sq=zeros(model.fn,model.fn);


model.x0tr=[6500.4,349.14,-1.8093,-6.7967,0.6932]';
%%%%%%%%%%%%% Run conditions 1 %%%%%%%%%%%%%%%%%%%%%%%%%
model.P0tr=diag([1e-1,1e-1,1e-1,1e-1,1]);
% for N=1:1:10
% X_mc=zeros(time.nSteps,model.fn,5000);
% x0trr = mvnrnd(model.x0tr', model.P0tr,5000);
% for i=1:1:5000
% 
%     [t,x_mc]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0trr(i,:)',model.Qt_sq);
%     X_mc(:,:,i)=x_mc;
% %       plot(X_mc(:,1),X_mc(:,2))
% %        axis([6360,6520,0,360])
% end
% save(strcat('COND_1_MC_runs_no_',num2str(N),'__',num2str(5000)),'model','X_mc')

% %UT runs
% kappa=1;
% [x0,w] = UT_sigmapoints(model.x0tr,model.P0tr,2);
% X=zeros(time.nSteps,model.fn,size(x0,1));
% for i=1:1:size(x0,1)
%     [t,x]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0(i,:)',model.Qt_sq);
%     X(:,:,i)=x;
% end
%  save(strcat('COND_1_UT_runs_no_',num2str(1),'__',num2str(size(x0,1))),'model','X','w')
% 
% % ckf
% [x0,w] = cubature_KF_points(model.x0tr,model.P0tr);
% X=zeros(time.nSteps,model.fn,size(x0,1));
% for i=1:1:size(x0,1)
%     [t,x]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0(i,:)',model.Qt_sq);
%     X(:,:,i)=x;
% end
%  save(strcat('COND_1_ckf_runs_no_',num2str(1),'__',num2str(size(x0,1))),'model','X','w')
% 
%  % cut4
% [x0,w]  = conjugate_dir_gausspts(model.x0tr,model.P0tr);
% X=zeros(time.nSteps,model.fn,size(x0,1));
% for i=1:1:size(x0,1)
%     [t,x]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0(i,:)',model.Qt_sq);
%     X(:,:,i)=x;
% end
%  save(strcat('COND_1_cut4_runs_no_',num2str(1),'__',num2str(size(x0,1))),'model','X','w')
% 
%  % cut6
% [x0,w]  = conjugate_dir_gausspts_till_6moment_scheme2(model.x0tr,model.P0tr);
% X=zeros(time.nSteps,model.fn,size(x0,1));
% for i=1:1:size(x0,1)
%     [t,x]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0(i,:)',model.Qt_sq);
%     X(:,:,i)=x;
% end
%  save(strcat('COND_1_cut6_runs_no_',num2str(1),'__',num2str(size(x0,1))),'model','X','w')
% 
%  % cut8
% [x0,w]  = conjugate_dir_gausspts_till_8moment(model.x0tr,model.P0tr);
% X=zeros(time.nSteps,model.fn,size(x0,1));
% for i=1:1:size(x0,1)
%     [t,x]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0(i,:)',model.Qt_sq);
%     X(:,:,i)=x;
% end
%  save(strcat('COND_1_cut8_runs_no_',num2str(1),'__',num2str(size(x0,1))),'model','X','w')
% 
%  
% gh
[x0,w]  = GH_points(model.x0tr,model.P0tr,5);
X=zeros(time.nSteps,model.fn,size(x0,1));
for i=1:1:size(x0,1)
    [t,x]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0(i,:)',model.Qt_sq);
    X(:,:,i)=x;
end
 save(strcat('COND_1_gh_runs_no_',num2str(1),'__',num2str(size(x0,1))),'model','X','w')


%%%%%%%%%%%%% Run conditions 2 %%%%%%%%%%%%%%%%%%%%%%%%%
% model.P0tr=diag([1,1,1e-2,1e-2,1]);
% for N=1:1:10
% X_mc=zeros(time.nSteps,model.fn,5000);
% x0trr = mvnrnd(model.x0tr', model.P0tr,5000);
% for i=1:1:5000
% 
%     [t,x_mc]=ode45_discc(model.fx,time.t0,time.dt,time.tf,x0trr(i,:)',model.Qt_sq);
%     X_mc(:,:,i)=x_mc;
%        
% end
% save(strcat('COND_2_MC_runs_no_',num2str(N),'__',num2str(5000)),'model','X_mc')
% end
%%%%%%%%%%% %  %%%%%%%%%%%%%%%% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(MU(:,1,9),MU(:,2,9),'k')
end