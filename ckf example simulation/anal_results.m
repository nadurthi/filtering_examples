clear
clc
p=31;
t=0:1:100;
scrsz = get(0,'ScreenSize');
for para_set=1:1:100

    load(strcat('ckf_eg_run_nos_',num2str(10),'paraset_',num2str(para_set)));
 figure(p)
set(gca,'Position',scrsz)
p=p+1;
subplot(2,3,1)
plot(t,x_fin_mc(:,1),t,x_fin_ckf(:,1),'LineWidth',2)
legend('MC','ckf')
xlabel('t')
ylabel('x')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)


subplot(2,3,2)
plot(t,x_fin_mc(:,1),t,x_fin_ukf(:,1),'LineWidth',2)
legend('MC','ukf')
xlabel('t')
ylabel('x')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,3,3)
plot(t,x_fin_mc(:,1),t,x_fin_cut4(:,1),'LineWidth',2)
legend('MC','cut4')
xlabel('t')
ylabel('x')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,3,4)
plot(t,x_fin_mc(:,1),t,x_fin_cut6(:,1),'LineWidth',2)
legend('MC','cut6')
xlabel('t')
ylabel('x')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,3,5)
plot(t,x_fin_mc(:,1),t,x_fin_cut8(:,1),'LineWidth',2)
legend('MC','cut8')
xlabel('t')
ylabel('x')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,3,6)
plot(t,x_fin_mc(:,1),t,x_fin_gh(:,1),'LineWidth',2)
legend('MC','GH')
xlabel('t')
ylabel('x')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
saveas(gcf,strcat('ckf eg Run no ',num2str(i),'  pic1','.jpg'))
close

figure(p)
set(gca,'Position',scrsz)
p=p+1;
subplot(2,3,1)
plot(t,x_fin_mc(:,3),t,x_fin_ckf(:,3),'LineWidth',2)
legend('MC','ckf')
xlabel('t')
ylabel('y')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,3,2)
plot(t,x_fin_mc(:,3),t,x_fin_ukf(:,3),'LineWidth',2)
legend('MC','ukf')
xlabel('t')
ylabel('y')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,3,3)
plot(t,x_fin_mc(:,3),t,x_fin_cut4(:,3),'LineWidth',2)
legend('MC','cut4')
xlabel('t')
ylabel('y')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,3,4)
plot(t,x_fin_mc(:,3),t,x_fin_cut6(:,3),'LineWidth',2)
legend('MC','cut6')
xlabel('t')
ylabel('y')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,3,5)
plot(t,x_fin_mc(:,3),t,x_fin_cut8(:,3),'LineWidth',2)
legend('MC','cut8')
xlabel('t')
ylabel('y')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,3,6)
plot(t,x_fin_mc(:,3),t,x_fin_gh(:,3),'LineWidth',2)
legend('MC','GH')
xlabel('t')
ylabel('y')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
saveas(gcf,strcat('ckf eg Run no ',num2str(i),'  pic2','.jpg'))
close

figure(p)
set(gca,'Position',scrsz)
p=p+1;
subplot(2,1,1)
plot(t,sqrt((x_fin_mc(:,1)-x_fin_ukf(:,1)).^2+(x_fin_mc(:,3)-x_fin_ukf(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_ckf(:,1)).^2+(x_fin_mc(:,3)-x_fin_ckf(:,3)).^2)...
     ,t,sqrt((x_fin_mc(:,1)-x_fin_cut4(:,1)).^2+(x_fin_mc(:,3)-x_fin_cut4(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_cut6(:,1)).^2+(x_fin_mc(:,3)-x_fin_cut6(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_cut8(:,1)).^2+(x_fin_mc(:,3)-x_fin_cut8(:,3)).^2),t,sqrt((x_fin_mc(:,1)-x_fin_gh(:,1)).^2+(x_fin_mc(:,3)-x_fin_gh(:,3)).^2),'LineWidth',2)
 legend('ukf','ckf','cut4','cut6','cut8','gh')
set(gca,'FontSize',7)
title(strcat('ckf eg run no',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,1,2)
plot(x_fin_mc(:,1),x_fin_mc(:,3),ym(:,1).*cos(ym(:,2)),ym(:,1).*sin(ym(:,2)),'k*','LineWidth',2)
legend('MC','sensor')
xlabel('x')
ylabel('y')
set(gca,'FontSize',7)
title(strcat('ckf eg Run no ',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
% set(gcf, 'PaperPositionMode', 'auto');
% set(gcf, 'PaperSize', [28,10]);
saveas(gcf,strcat('ckf eg Run no ',num2str(i),'  pic3','.jpg'))
close

 figure(p)
set(gca,'Position',scrsz)
p=p+1;
for jj=1:1:length(t)
Pukf=reshape(P_fin_ukf(jj,:),5,5);
[u,s,v]=svd(Pukf);
s_ukf(jj,:)=sqrt(sum(s,1));

Pckf=reshape(P_fin_ckf(jj,:),5,5);
[u,s,v]=svd(Pckf);
s_ckf(jj,:)=sqrt(sum(s,1));

Pcut4=reshape(P_fin_cut4(jj,:),5,5);
[u,s,v]=svd(Pcut4);
s_cut4(jj,:)=sqrt(sum(s,1));

Pcut6=reshape(P_fin_cut6(jj,:),5,5);
[u,s,v]=svd(Pcut6);
s_cut6(jj,:)=sqrt(sum(s,1));

Pcut8=reshape(P_fin_cut8(jj,:),5,5);
[u,s,v]=svd(Pcut8);
s_cut8(jj,:)=sqrt(sum(s,1));

Pgh=reshape(P_fin_gh(jj,:),5,5);
[u,s,v]=svd(Pgh);
s_gh(jj,:)=sqrt(sum(s,1));

end
subplot(2,1,1)
plot(t,s_ukf(:,1),t,s_ckf(:,1),t,s_cut4(:,1),t,s_cut6(:,1),t,s_cut8(:,1),t,s_gh(:,1),'LineWidth',2)
 legend('ukf','ckf','cut4','cut6','cut8','gh')
xlabel('t')
ylabel('x')
set(gca,'FontSize',7)
title(strcat('ckf eg Run no ',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)

subplot(2,1,2)
plot(t,s_ukf(:,3),t,s_ckf(:,3),t,s_cut4(:,3),t,s_cut6(:,3),t,s_cut8(:,3),t,s_gh(:,3),'LineWidth',2)
 legend('ukf','ckf','cut4','cut6','cut8','gh')
xlabel('t')
ylabel('y')
set(gca,'FontSize',7)
title(strcat('ckf eg Run no ',num2str(i)))
h = get(gca, 'title');
k = get(gca, 'xlabel');
l = get(gca, 'ylabel');
set(h, 'FontName', 'Helvetica', 'FontSize', 14)
set(k, 'FontName', 'Helvetica', 'FontSize', 16)
set(l, 'FontName', 'Helvetica', 'FontSize', 16)
saveas(gcf,strcat('ckf eg Run no ',num2str(i),'  pic4','.jpg'))
close

     clear(strcat('ckf_eg_run_no_',num2str(i)))
end