clear
clc
F=@(x1,x2,x3)(x1^2*x2^2+x2^2*x3^2+x1^2*x3^2);
w=[0.3583,0.0045,0.0045,0.0045,0.0045,0.0045,0.0045,0.0512,0.0512,0.0512,0.0512,0.0512,0.0512,0.0512,0.0512,0.0512,0.0512,0.0512,0.0512];
w=w';
s1=3.2530;
s2=1.4862;
x=[0,0,0;s1,0,0;-s1,0,0;0,s1,0;0,-s1,0;0,0,s1;0,0,-s1;s2,s2,0;-s2,-s2,0;s2,-s2,0;-s2,s2,0;s2,0,s2;-s2,0,-s2;s2,0,-s2;-s2,0,s2;0,s2,s2;0,-s2,-s2;0,-s2,s2;0,s2,-s2];

P=[4,2,1;2,5,3;1,3,6];

[u,s,v]=svd(inv(P));
u=v';
A=sqrtm(s)*u;
y1=[];
for i=1:1:length(w);
    r=inv(A)*x(i,:)';
    y1=[y1;r'];
end
sumjpts=0;
xj=y1;
for i=1:1:length(w)
    sumjpts=sumjpts+w(i)*F(xj(i,1),xj(i,2),xj(i,3));
end


[xgh,wgh]=conjugate_dir_gausspts([0;0;0],eye(3));
y2=[];
for i=1:1:length(wgh);
    r=inv(A)*xgh(i,:)';
    y2=[y2;r'];
end
sumgh=0;
xgh=y2;
for i=1:1:length(wgh)
    sumgh=sumgh+wgh(i)*F(xgh(i,1),xgh(i,2),xgh(i,3));
end
[sumjpts,sumgh]