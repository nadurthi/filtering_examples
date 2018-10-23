function K=kernelpq2D(p,q,u,v,X,w)
%% paeno error kernels for 2D.
global err_tol
% the boundary is b/n -1 and 1
% so u would have to shift the function to theis region
% the taylor series is about origin 0,0
a0=0;
c0=0;
facp=factorial(p-1);
facq=factorial(q-1);
% err_tol=10e-6;
f=@(x,y,p,q,u,v)(((x-u)^(p-1)/facp)*step_fn(a0,u,x)*(y-v)^(q-1)/facq*step_fn(c0,v,y));
Ip=1000000;
for N=2:2:200
[Xt,wt]=GLeg_pts(N*ones(1,2),-1*ones(1,2),1*ones(1,2));
It=0;
for i=1:1:length(wt)
    It=It+wt(i)*f(Xt(i,1),Xt(i,2),p,q,u,v);
end
if abs(It-Ip)<=err_tol
    break
end
Ip=It;
end

Im=0;
for i=1:1:length(w)
    Im=Im+w(i)*f(X(i,1),X(i,2),p,q,u,v);
end

K=It-Im;