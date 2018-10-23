function dx=reentryvehicle(t,x,paras)
% if the para is in the cell sof trings paras then it is considered as
% parameter and not a state
% x=[r,v,gam,p0,Bc,g,LD]


r=x(1);
v=x(2);
gam=x(3);

k=4;

if  cell2mat(  strfind(paras,'p0')  )
  p0=0.0019;
else
    p0=x(k);
    k=k+1;
end

if  cell2mat(  strfind(paras,'Bc')  )
  Bc=72.8;
else
    Bc=x(k);
    k=k+1;
end

if  cell2mat(  strfind(paras,'g')  )
  g=3.71;
else
    g=x(k);
    k=k+1;
end

if  cell2mat(  strfind(paras,'LD')  )
  LD=0.3;
else
    LD=x(k);
    k=k+1;
end

Rm=3397;
h1=9.8;
h2=20;
h=r-Rm;
p=p0*exp((h2-h)/h1);

dx(1,1)=v*sin(gam);
dx(2,1)=-p*v^2/(2*Bc)-g*sin(gam);
dx(3,1)=(v/r-g/v)*cos(gam)+1/(2*Bc)*LD*v;







