function [xu,Pu]=EKF_disc(model,xu,Pu,ym)
%this is discrete dynamics
% global Js xs n q
% Js=F;
% xs=xu;
% q=Q;

% n=length(xu);
% forecast

%given xu and Pu at time k ....this func gives xu and Pu at time k+1

%ym is meas at this time step
% the struct model should have the following
%model.fx is the model dyn function....as model.fx=@Dynamic_model etc...
%model.Q is procss noise cov 
%model.fx_jac is jacobian of dyns .... model.fx_jac=@Dynamic_model_jacobian etc...
%model.hx is meas model.... model.hx=@Meas_model etc...
%model.hx_jac is meas model.... model.hx_jac=@Meas_model_jac etc...
%model.R is meas noise cov
%model.fn  is the dimension of dyn system
%model.hn is dim of meas model
% the model.fx takes in a column vector and gives out a column vector...same with all function meas etc
% jac takes a column vector and gives out a square matrix

%if ym=-1234 that means dont do meas update

xf=model.fx(xu,model.para_dt);
BB=model.fx_jac(xu,model.para_dt);
Pf=BB*Pu*BB'+model.Q;

%update
if ym==-1234
    xu=xf;
    Pu=Pf;
else
 
    AA=model.hx_jac(xf,model.para_dt);
    K=Pf*AA'*inv(AA*Pf*AA'+model.R);
    xu=xf+K*(ym-model.hx(xf,model.para_dt));
    Pu=(eye(model.fn)-K*AA)*Pf;
end
end

% function dp=cov_prop(t,p)
% global Js xs n q
% F=Js;
% xu=xs;
% Q=q;
% P=reshape(p,n,n);
% dpp=F(xu)*P+P*F(xu)'+Q;
% dp=reshape(dpp,1,n*n);
% dp=dp';
% end