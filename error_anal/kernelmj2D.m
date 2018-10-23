function K=kernelmj2D(m,j,u,X,w)
%% paeno error kernels for 2D.
% the boundary is b/n -1 and 1
% so u would have to shift the function to theis region
% the taylor series is about origin 0,0
global err_tol
a0=0;
c0=0;
facmj=factorial(m-j-1);
facj=factorial(j);
% err_tol=10e-6;
f=@(x,y,m,j,u)(((x-u)^(m-j-1)/facmj)*step_fn(a0,u,x)*(y-c0)^(j)/facj);
Ip=1000000;
for N=2:2:200
[Xt,wt]=GLeg_pts(N*ones(1,2),-1*ones(1,2),1*ones(1,2));
It=0;
for i=1:1:length(wt)
    It=It+wt(i)*f(Xt(i,1),Xt(i,2),m,j,u);
end
if abs(It-Ip)<=err_tol
    break
end
Ip=It;
end

Im=0;
for i=1:1:length(w)
    Im=Im+w(i)*f(X(i,1),X(i,2),m,j,u);
end

K=It-Im;