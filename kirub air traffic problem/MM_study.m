x0=[10,10,0,10,-15*pi/180]';
X=x0';
p=1;
for t=1:1:10
    xk=X(end,:)';
%    keyboard
xk1=(1-p)*KIRB_CT_eg_dyn_disc(xk,1)+p*[KIRB_UM_eg_dyn_disc(xk(1:4,1),1);0];
X=vertcat(X,xk1');
end
plot(X(:,1),X(:,2),'m')
hold on
axis([0,100,0,100])