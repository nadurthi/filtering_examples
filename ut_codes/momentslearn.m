clc
clear
w1=0.7;
w2=1-0.7;
mu1=-5;
mu2=0.5;
p1=1;
p2=1;
m1=w1*mu1+w2*mu2;
m2=w1*(mu1^2+p1)+w2*(mu2^2+p2);
m3=w1*(mu1^3+3*mu1*p1)+w2*(mu2^3+3*mu2*p2);
m4=w1*(mu1^4+6*mu1^2*p1+3*p1^2)+w2*(mu2^4+6*mu2^2*p2+3*p2^2);
m5=w1*(mu1^5+10*mu1^3*p1+15*mu1*p1^2)+w2*(mu2^5+10*mu2^3*p2+15*mu2*p2^2);
[m1,m2,m3,m4,m5]
mc1=m1;
pc=m2-m1^2;
mc2=mc1^2+pc;
mc3=mc1^3+3*mc1*pc;
mc4=mc1^4+6*mc1^2*pc+3*pc^2;
mc5=mc1^5+10*mc1^3*pc+15*mc1*pc^2;
[mc1,mc2,mc3,mc4,mc5]

X=-50:0.1:50;


    
Y1 = normpdf(X,mu1,sqrt(p1));
Y2 = normpdf(X,mu2,sqrt(p2));
Y3 = normpdf(X,mc1,sqrt(pc));
Y=w1*Y1+w2*Y2;

plot(X,Y,X,Y3)
axis([-20,20,0,1])
