clc
clear
global err_tol
%% define the cubature method for which u wanna
% the error
% [X,w]=GLeg_pts(2*ones(1,2),-1*ones(1,2),1*ones(1,2));
[X,w]=uniform_sigma_pts(-1*ones(1,2),1*ones(1,2),2);
%% error_anals
f=@(x,y)10*cos(x+y);
err_tol=10e-4;
%% define the function space
m=4;
p=2;
q=2;

%% calculate the 1-norm of the kernelspq
Ip=1000000;
for N=2:2:200
[Xt,wt]=GLeg_pts(N*ones(1,2),-1*ones(1,2),1*ones(1,2));
epq=0;
for i=1:1:length(wt)
    epq=epq+wt(i)*abs(kernelpq2D(p,q,Xt(i,1),Xt(i,2),X,w));
end
if abs(epq-Ip)<=err_tol
    break
end
Ip=epq;
end

%% calculate the 1-norm of the kernels mj
for j=0:1:q-1
Ip=1000000;
for N=2:2:200
[Xt,wt]=GLeg_pts(N,-1,1);
emj(j+1)=0;
for i=1:1:length(wt)
    emj(j+1)=emj(j+1)+wt(i)*abs(kernelmj2D(m,j,Xt(i,1),X,w));
end
if abs(emj(j+1)-Ip)<=err_tol
    break
end
Ip=emj(j+1);
end

end

%% calculate the 1-norm of the kernels
for i=0:1:p-1
Ip=1000000;
for N=2:2:200
[Xt,wt]=GLeg_pts(N,-1,1);
emi(i+1)=0;
for k=1:1:length(wt)
    emi(i+1)=emi(i+1)+wt(k)*abs(kernelmi2D(m,i,Xt(k,1),X,w));
end
if abs(emi(i+1)-Ip)<=err_tol
    break
end
Ip=emi(i+1);
% keyboard
end

end

%% the error bound is
Ib=sum(emj)+sum(emi)+epq
%% computing the actual error bound
Ip=1000000;
for N=2:2:200
[Xt,wt]=GLeg_pts(N*ones(1,2),-1*ones(1,2),1*ones(1,2));
It=0;
for i=1:1:length(wt)
    It=It+wt(i)*f(Xt(i,1),Xt(i,2));
end
if abs(It-Ip)<=err_tol
    break
end
Ip=It;
end

Ia=0;
for i=1:1:length(w)
        Ia=Ia+w(i)*f(X(i,1),X(i,2));
end

TrI=abs(It-Ip)
epq
emj
emi