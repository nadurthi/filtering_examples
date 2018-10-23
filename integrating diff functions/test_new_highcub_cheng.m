n=2;

% A=randn(n);
% P=(A'*A).^2/20;
% eig(P)
P=eye(n);
mu=zeros(n,1);

[xut,wut]=UT_sigmapoints(mu,P,2);

[xcheng5,wcheng5]=high_cub_chengpaer(mu,P);
[xcheng7,wcheng7]=high_cub_cheng_paper(mu,P,7);
[xcheng9,wcheng9]=high_cub_cheng_paper(mu,P,9);

[xcut4,wcut4]=conjugate_dir_gausspts(mu,P);
[xcut6,wcut6]=conjugate_dir_gausspts_till_6moment_scheme2(mu,P);
[xcut8,wcut8]=conjugate_dir_gausspts_till_8moment(mu,P);

[xgh3,wgh3]=GH_points(mu,P,3);
[xgh4,wgh4]=GH_points(mu,P,4);
[xgh5,wgh5]=GH_points(mu,P,5);
[xgh6,wgh6]=GH_points(mu,P,6);
[xgh7,wgh7]=GH_points(mu,P,7);
[xgh8,wgh8]=GH_points(mu,P,8);

[xsp2,ws]=smolyak_sparse_grid_modf(mu,P,n,3,'gh');
min(ws)
sum(ws)
length(ws)

[xng wng] = nwSpGr('GQN', n, 3);
length(wng)
min(wng)
sum(abs(ws))

[xsp3,wsp3]=smolyak_sparse_grid_modf(mu,P,n,3,'GH');
[xsp4,wsp4]=smolyak_sparse_grid_modf(mu,P,n,4,'GH');
[xsp5,wsp5]=smolyak_sparse_grid_modf(mu,P,n,5,'GH');
[xsp6,wsp6]=smolyak_sparse_grid_modf(mu,P,n,6,'GH');
[xsp7,wsp7]=smolyak_sparse_grid_modf(mu,P,n,7,'GH');
[xsp8,wsp8]=smolyak_sparse_grid_modf(mu,P,n,8,'GH');

mut=cal_moments_wrt_pts(xut,wut,8);
mcheng5=cal_moments_wrt_pts(xcheng5,wcheng5,8);
mcheng7=cal_moments_wrt_pts(xcheng7,wcheng7,8);
mcheng9=cal_moments_wrt_pts(xcheng9,wcheng9,8);

mcut4=cal_moments_wrt_pts(xcut4,wcut4,8);
mcut6=cal_moments_wrt_pts(xcut6,wcut6,8);
mcut8=cal_moments_wrt_pts(xcut8,wcut8,8);

mgh3=cal_moments_wrt_pts(xgh3,wgh3,8);
msp3=cal_moments_wrt_pts(xsp3,wsp3,8);

mgh7=cal_moments_wrt_pts(xgh6,wgh6,8);

SS=[mut(:,end),mcheng5(:,end),mcheng7(:,end),mcheng9(:,end),mcut4(:,end),mcut6(:,end),mgh3(:,end),msp3(:,end),mgh7(:,end)];

ind=find(sum(abs(SS),2)>=1);
SS(ind,:)
%  F=@(x)(1-sqrt(1+sum(abs(x).^2,2)).^(-3));
%  F=@(x)ones(size(x,1),1);
% F=@(x)sqrt(sum(abs(x).^1,2)); 
% F=@(x)(1-100*mvnpdf(x,ones(n,1)',0.5*eye(n)));
 F=@(x)cos(sqrt(sum(x.^2,2))); 
%   F=@(x)(sum(0.001*x.^10,2));
 
 Iut=sum(wut.*F(xut));
 
 Icheng5=sum(wcheng5.*F(xcheng5));
  Icheng7=sum(wcheng7.*F(xcheng7));
   Icheng9=sum(wcheng9.*F(xcheng9));
   
 Icut4=sum(wcut4.*F(xcut4));
  Icut6=sum(wcut6.*F(xcut6));
%    Icut8=sum(wcut8.*F(xcut8));
   
   Igh3=sum(wgh3.*F(xgh3));
   Igh4=sum(wgh4.*F(xgh4));
   Igh5=sum(wgh5.*F(xgh5));
   Igh6=sum(wgh6.*F(xgh6));
%    Igh7=sum(wgh7.*F(xgh7));
%    Igh8=sum(wgh8.*F(xgh8));
    
    Isp3=sum(wsp3.*F(xsp3));
   Isp4=sum(wsp4.*F(xsp4));
   Isp5=sum(wsp5.*F(xsp5));
%    Isp6=sum(wsp6.*F(xsp6));
%    Isp7=sum(wsp7.*F(xsp7));
%    Isp8=sum(wsp8.*F(xsp8));
   [Iut,Icheng5,Icheng7,Icut4,Icut6,Igh3,Igh4,Igh5,Igh6]
   
 KJ=[Iut,Icheng5,Icheng7,Icut4,Icut6,Icut8,Igh3,Igh4,Igh5,Igh6,Igh7,Igh8,Isp3,Isp4,Isp5,Isp6,Isp7,Isp8];
 NJ=[length(wut),length(wcheng),length(wcut4),length(wcut6),length(wcut8),length(wgh3),length(wgh4),length(wgh5),length(wgh6),length(wgh7),length(wgh8),length(wsp3),length(wsp4),length(wsp5),length(wsp6),length(wsp7),length(wsp8)];
 WJ=[min(wut),min(wcheng),min(wcut4),min(wcut6),min(wcut8),min(wgh3),min(wgh4),min(wgh5),min(wgh6),min(wgh7),min(wgh8),min(wsp3),min(wsp4),min(wsp5),min(wsp6),min(wsp7),min(wsp8)];
 SJ=[sum(abs(wut)),sum(abs(wcheng)),sum(abs(wcut4)),sum(abs(wcut6)),sum(abs(wcut8)),sum(abs(wgh3)),sum(abs(wgh4)),sum(abs(wgh5)),sum(abs(wgh6)),sum(abs(wgh7)),sum(abs(wgh8)),sum(abs(wsp3)),sum(abs(wsp4)),sum(abs(wsp5)),sum(abs(wsp6)),sum(abs(wsp7)),sum(abs(wsp8))];
 Itr=-0.543583844;
 MM=cell(17,5);
 MM{1,1}='Iut';MM{2,1}='Icheng';MM{3,1}='Icut4';MM{4,1}='Icut6';MM{5,1}='Icut8';MM{6,1}='Igh3';MM{7,1}='Igh4';MM{8,1}='Igh5';MM{9,1}='Igh6';MM{10,1}='Igh7';MM{11,1}='Igh8';MM{12,1}='Isp3';MM{13,1}='Isp4';MM{14,1}='Isp5';MM{15,1}='Isp6';MM{16,1}='Isp7';MM{17,1}='Isp8';
 for i=1:1:17
     MM{i,2}=KJ(i);
     %MM{i,3}=[abs(Itr-KJ(i))/abs(Itr)]*100
     MM{i,3}=WJ(i);
     MM{i,4}=SJ(i);
     MM{i,5}=NJ(i);
 end
 MM
 
ind=find(wcheng<0)
sqrt(sum(xcheng(ind,:).^2,2))
 
ind=find(wsp3<0)
sqrt(sum(xsp3(ind,:).^2,2))
length(ind)
length(wsp3)
length(ind)/length(wsp3)*100
nsp3=[wsp3(ind),F(xsp3(ind,:))];

ind=find(wsp4<0)
sqrt(sum(xsp4(ind,:).^2,2))
length(ind)
length(wsp4)
length(ind)/length(wsp4)*100
nsp4=[wsp4(ind),F(xsp4(ind,:))];

ind=find(wcheng<0)
sqrt(sum(xcheng(ind,:).^2,2))
length(ind)
length(wcheng)
length(ind)/length(wcheng)*100
ncheng=[wcheng(ind),F(xcheng(ind,:))];

[nsp3]
[nsp4]
[ncheng]


n=6;
[xcheng7,wcheng7,xs,ws]=high_cub_cheng_paper(zeros(n,1),eye(n),9);
length(wcheng7)
length(ws)

n=8;
[xcut6,wcut6]=conjugate_dir_gausspts_till_6moment_scheme2(zeros(n,1),eye(n));
length(wcut6)


[xsp3,wsp3]=smolyak_sparse_grid_modf(zeros(n,1),eye(n),n,m,'GH');
[nodes, weights] = nwspgr('GQN', n, m);
length(weights)
length(wsp3)

n=6;
m=5;
[nodes, weights] = nwspgr('GQN', n, m);
min(weights)


%%
[xgh4,wgh4]=GH_points(mu,P,4);
mgh4=cal_moments_wrt_pts(xgh4,wgh4,12);

[xgh8,wgh8]=GH_points(mu,P,8);
mgh8=cal_moments_wrt_pts(xgh8,wgh8,12);

round([mgh4,mgh8(:,end)])