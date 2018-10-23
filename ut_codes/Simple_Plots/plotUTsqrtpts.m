clc
 x1 = -10:.2:10; x2 = -10:.2:10;
[X1,X2] = meshgrid(x1,x2);
P=[5,2;2,3];
F = mvnpdf([X1(:) X2(:)],[0,0],P);
F = reshape(F,length(x2),length(x1));
% figure(1)
% surf(x1,x2,F);
% caxis([min(F(:))-0.5*range(F(:)),max(F(:))]);
% % axis([-10 10 -10 10])
% xlabel('x1'); 
% ylabel('x2');
% zlabel('Probability Density');
figure(2)
contour(X1,X2,F,[0.03])
[u,s,v]=svd(P);
hold on
plot(u(1,1).*x1,u(2,1).*x1)
plot(u(1,2).*x1,u(2,2).*x1)
w0=1-2/3;
w1=(1-w0)/4;
PP=P;
A=chol(PP);
% A=sqrtm((2/(1-w0))*P);
A=A';

[u,s,v]=svd(PP);
AA=chol(s)*u;
AA=AA'
plot(AA(1,1),AA(2,1),'rx',-AA(1,1),-AA(2,1),'rx')
plot(AA(1,2),AA(2,2),'ro',-AA(1,2),-AA(2,2),'ro')
plot(A(1,1),A(2,1),'kx',-A(1,1),-A(2,1),'kx')
plot(A(1,2),A(2,2),'ko',-A(1,2),-A(2,2),'ko')
v1=[AA(1,1);AA(2,1)];
v2=-[AA(1,1);AA(2,1)];
v3=[AA(1,2);AA(2,2)];
v4=-[AA(1,2);AA(2,2)];
acosd((v1-v2)'*(v3-v4)/(norm((v1-v2))*norm(v3-v4)))

w1*(v1*v1'+v2*v2'+v3*v3'+v4*v4')
