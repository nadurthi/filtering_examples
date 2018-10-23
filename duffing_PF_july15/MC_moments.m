clear all;close all;clc
load duff_mcdata_1e5
[mcr tl]=size(y1');
y1=y1'; y2=y2';
y1m=zeros(1,mcr); y1v=zeros(1,mcr); y1s=zeros(1,mcr);
y2m=zeros(1,mcr); y2v=zeros(1,mcr); y2s=zeros(1,mcr);

for ind=1:mcr
    y1m(ind)=sum(y1(1:ind,201))/ind;
    y2m(ind)=sum(y2(1:ind,201))/ind;
    y1v(ind)=sum((y1(1:ind,201)-y1m(ind)).^2)/ind;
    y2v(ind)=sum((y2(1:ind,201)-y2m(ind)).^2)/ind;
    y1s(ind)=sum((y1(1:ind,201)-y1m(ind)).^3)/ind;
    y2s(ind)=sum((y2(1:ind,201)-y2m(ind)).^3)/ind;
    ind
end
figure;box on; hold on
set(gca,'Fontsize',16);
plot(y1m,'LineWidth',1.5);
xlabel('Number of Monte Carlo Runs');
ylabel('E[x]');

figure;box on; hold on
set(gca,'Fontsize',16);
plot(y2m,'LineWidth',1.5);
xlabel('Number of Monte Carlo Runs');
ylabel('E[dx/dt]');


figure;box on; hold on
set(gca,'Fontsize',16);
plot(y1v,'LineWidth',1.5);
xlabel('Number of Monte Carlo Runs');
ylabel('E[(x-E[x])^2]');

figure;box on; hold on
set(gca,'Fontsize',16);
plot(y2v,'LineWidth',1.5);
xlabel('Number of Monte Carlo Runs');
ylabel('E[(dx/dt-E[dx/dt])^2]');   


figure;box on; hold on
set(gca,'Fontsize',16);
plot(y1s,'LineWidth',1.5);
xlabel('Number of Monte Carlo Runs');
ylabel('E[(x-E[x])^3]');

figure;box on; hold on
set(gca,'Fontsize',16);
plot(y2s,'LineWidth',1.5);
xlabel('Number of Monte Carlo Runs');
ylabel('E[(dx/dt-E[dx/dt])^3]');   
