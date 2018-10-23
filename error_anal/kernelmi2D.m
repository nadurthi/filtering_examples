function K=kernelmi2D(m,i,v,X,w)
%% paeno error kernels for 2D.
% the boundary is b/n -1 and 1
% so u would have to shift the function to theis region
% the taylor series is about origin 0,0
global err_tol
a0=0;
c0=0;
facmi=factorial(m-i-1);
faci=factorial(i);
% err_tol=10e-6;
f=@(x,y,m,i,v)(((x-a0)^(i)/faci)*(y-v)^(m-i-1)/facmi*step_fn(c0,v,y));
Ip=1000000;
for N=2:2:200
[Xt,wt]=GLeg_pts(N*ones(1,2),-1*ones(1,2),1*ones(1,2));
It=0;
for k=1:1:length(wt)
    It=It+wt(k)*f(Xt(k,1),Xt(k,2),m,i,v);
end
if abs(It-Ip)<=err_tol
    break
end
Ip=It;
end

Im=0;
for i=k:1:length(w)
    Im=Im+w(k)*f(X(k,1),X(k,2),m,i,v);
end

K=It-Im;