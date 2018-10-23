I=0;
for i=2:20
[xint,wint] = GH_points(0,1,i);
I(i-1)=sum(wint.*fun(xint));
end
plot(2:20,I)
[xint,wint] = GH_points(0,1,3);
sum(wint.*fun(xint))

x=-10:0.05:10;
y=1- 0.5*x.^2 + 0.375*x.^4 - 0.3125*x.^6 + 0.273438*x.^8 -... 
 0.246094*x.^10 + 0.225586*x.^12 - 0.209473*x.^14 + 0.196381*x.^16 -... 
 0.185471*x.^18 + 0.176197*x.^20;

% figure
% plot(x,1./sqrt(1+x.^2),'r',x,y,'b','linewidth',2)
% legend('f(x)','Taylor Series of f upto 20 terms')
% hold on
% plot([0,0],[-1,1],'k--',[-10,10],[0,0],'k--','linewidth',2)  
%   axis([-10,10,0,1])
%   plot_prop_paper
  
  figure
  plot(x,fun(x),'r',x,pn_bpa(x),'linewidth',2)
legend('f(x)','Best Polynomial approximation')
hold on
plot([0,0],[0,8],'k--',[-10,10],[0,0],'k--','linewidth',2)  
plot([5,5],[0,8],'k-.',[-5,-5],[0,8],'k-.','linewidth',1)  
axis([-7,7,0,8])
plot_prop_paper