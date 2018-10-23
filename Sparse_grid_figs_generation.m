close all

figure
[x,w]=GH_points(0,1,1);
plot(x,0,'bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)
axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('y')
title('Q_1^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,3);
plot([0,0,0],x,'bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)
axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('y')
title('Q_3^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,3);
plot([0,0,0],x,'bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)
          
axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('y')
title('Q_1^1 \otimes Q_3^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,3);
plot3([0,0,0],x,w,'bo',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)
           hold on
            plot_dotlines_pts([[0;0;0],x],w,'b--')
 fill3([-2,-2,2,2],[-2,2,2,-2],[0,0,0,0],[0.9,0.9,0.9])
 alpha(0.5)
axis([-2,2,-2,2,-1,1])
xlabel('x')
grid
ylabel('y')
zlabel('w')
title('Q_1^1 \otimes Q_3^1')
plot_prop_paper


%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
[x,w]=GH_points(0,1,2);
[x1,w1]=GH_points(0,1,1);
plot(x,w,'bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)      
axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('w')
title('Q_2^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,2);
[x1,w1]=GH_points(0,1,1);
plot(-x1,-w1,'ro',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',12)            
axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('w')
title('-Q_1^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,2);
[x1,w1]=GH_points(0,1,1);
plot(x,[0,0],'bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)   
            hold on
plot(-x1,0,'ro',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)            
axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('y')
title('Q_2^1-Q_1^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,2);

plot([0,0],x,'bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)   
        
axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('y')
title(' Q_2^1')
plot_prop_paper



figure
[x,w]=GH_points(0,1,2);
[x1,w1]=GH_points(0,1,1);

[xd,wd]=tens_prod_vec([x;0],x,[w;-1],w);
plot(xd(:,1),xd(:,2),'bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)   
 grid
axis([-2,2,-2,2])
xlabel('x')

ylabel('y')
title('(Q_2^1-Q_1^1)\otimes Q_2^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,2);
[x1,w1]=GH_points(0,1,1);

[xd,wd]=tens_prod_vec([x;0],x,[w;-1],w);
plot3(xd(1:4,1),xd(1:4,2),wd(1:4),'bo',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)   
 grid
 hold on
 plot_dotlines_pts(xd(1:4,:),wd(1:4),'b--')
 plot3(xd(5:6,1),xd(5:6,2),wd(5:6),'ro',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',12)  
            plot_dotlines_pts(xd(5:6,:),wd(5:6),'r--')
 fill3([-2,-2,2,2],[-2,2,2,-2],[0,0,0,0],[0.9,0.9,0.9])
 alpha(0.5)
axis([-2,2,-2,2,-1,1])
xlabel('x')

ylabel('y')
zlabel('w')
title('(Q_2^1-Q_1^1)\otimes Q_2^1')
plot_prop_paper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
[x,w]=GH_points(0,1,3);
plot(x,w,'bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)
            hold on
            plot_dotlines_pts(x,w,'b--')
axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('w')
title('Q_3^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,2);
plot(x,-w,'bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',12)
             hold on
            plot_dotlines_pts(x,-w,'r--')
axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('w')
title('Q_2^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,2);
[x1,w1]=GH_points(0,1,3);
plot(x(:,1),[0,0]','bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)
             hold on
plot(x1(:,1),[0,0,0]','bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)
axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('y')
title('Q_3^1-Q_2^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,2);
[x1,w1]=GH_points(0,1,3);
plot([0,0,0,0,0],[x;x1],'bo',[-2,2],[0,0],'k--',[0,0],[-2,2],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)

axis([-2,2,-2,2])
xlabel('x')
grid
ylabel('y')
title('(Q_3^1-Q_2^1)\otimes Q_1^1')
plot_prop_paper

figure
[x,w]=GH_points(0,1,3);
[x1,w1]=GH_points(0,1,2);
[x2,w2]=GH_points(0,1,1);
[xd,wd]=tens_prod_vec([x;-x1],x2,[w;-w1],w2);
plot3(xd(1:3,1),xd(1:3,2),wd(1:3),'bo',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)   
            hold on
            plot3(xd(4:5,1),xd(4:5,2),wd(4:5),'bo',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',12)   
            
 grid
  plot_dotlines_pts(xd(1:3,:),wd(1:3),'b--')
  plot_dotlines_pts(xd(4:5,:),wd(4:5),'r--')
 fill3([-2,-2,2,2],[-2,2,2,-2],[0,0,0,0],[0.9,0.9,0.9])
 alpha(0.5)
axis([-2,2,-2,2,-1,1])
xlabel('x')
zlabel('w')
ylabel('y')
title('(Q_3^1-Q_2^1)\otimes Q_1^1')
plot_prop_paper

%%%%%%%%%%%%%%%%
%combining them all together
figure
[x,w]=GH_points(0,1,2);
[x1,w1]=GH_points(0,1,1);

[xd,wd]=tens_prod_vec([x;0],x,[w;-1],w);
plot3(xd(1:4,1),xd(1:4,2),wd(1:4),'bo',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)   
 
 hold on
 plot_dotlines_pts(xd(1:4,:),wd(1:4),'b--')
 plot3(xd(5:6,1),xd(5:6,2),wd(5:6),'ro',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',12)  
            plot_dotlines_pts(xd(5:6,:),wd(5:6),'r--')
%  fill3([-2,-2,2,2],[-2,2,2,-2],[0,0,0,0],[0.9,0.9,0.9])
%  alpha(0.5)
axis([-2,2,-2,2,-1,1])
xlabel('x')

ylabel('y')
zlabel('w')
% title('(Q_2^1-Q_1^1)\otimes Q_2^1')

[x,w]=GH_points(0,1,3);
[x1,w1]=GH_points(0,1,2);
[x2,w2]=GH_points(0,1,1);
[xd,wd]=tens_prod_vec([x;-x1],x2,[w;-w1],w2);
plot3(xd(1:3,1),xd(1:3,2),wd(1:3),'bo',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)   
%             hold on
            plot3(xd(4:5,1),xd(4:5,2),wd(4:5),'bo',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',12)   
            
 grid
  plot_dotlines_pts(xd(1:3,:),wd(1:3),'b--')
  plot_dotlines_pts(xd(4:5,:),wd(4:5),'r--')
 fill3([-2,-2,2,2],[-2,2,2,-2],[0,0,0,0],[0.9,0.9,0.9])
 alpha(0.35)
axis([-2,2,-2,2,-1,1])
xlabel('x')
zlabel('w')
ylabel('y')
title('Q_1^1 \otimes Q_3^1+(Q_2^1-Q_1^1)\otimes Q_2^1+(Q_3^1-Q_2^1)\otimes Q_1^1')

[x,w]=GH_points(0,1,3);
plot3([0,0,0],x,w,'bo',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)
%            hold on
            plot_dotlines_pts([[0;0;0],x],w,'b--')
%  fill3([-2,-2,2,2],[-2,2,2,-2],[0,0,0,0],[0.9,0.9,0.9])
%  alpha(0.5)
% axis([-2,2,-2,2,-1,1])
% xlabel('x')
% grid
% ylabel('y')
% zlabel('w')
% title('')

plot_prop_paper


figure
[x,w]=smolyak_sparse_grid(2,3,'GH');
plot3(x([1,2,3,4,5,6,7,10,11],1),x([1,2,3,4,5,6,7,10,11],2),w([1,2,3,4,5,6,7,10,11]),'bo',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',12)   
            hold on
            plot3(x([8,9,12,13],1),x([8,9,12,13],2),w([8,9,12,13]),'bo',[-2,2],[0,0],[0,0],'k--',[0,0],[-2,2],[0,0],'k--','MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',12)   
            
 grid
  plot_dotlines_pts(x([1,2,3,4,5,6,7,10,11],:),w([1,2,3,4,5,6,7,10,11]),'b--')
   plot_dotlines_pts(x([8,9,12,13],:),w([8,9,12,13]),'r--')
 fill3([-2,-2,2,2],[-2,2,2,-2],[0,0,0,0],[0.9,0.9,0.9])
 alpha(0.35)
axis([-2,2,-2,2,-1,1])
xlabel('x')
zlabel('w')
ylabel('y')
title('Q_3^2')
plot_prop_paper