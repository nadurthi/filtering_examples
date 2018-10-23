 [ xp, wp ] =patterson_rule ( 15, -1, 1 );
 
 [x,w]=GLeg_pts(4, -1, 1);
 w=2*w;
 
 [sum(wp.*xp.^12),sum(w.*x.^10)]
 
 n=4
 [xgh,wgh]=GH_points(zeros(5,1),eye(5),5);
 [xnew, wnew ] = gqn ( n );
 
 
 n=9;
P=eye(n);
mu=zeros(n,1);
[xcut4,wcut4]=conjugate_dir_gausspts(mu,P);
length(wcut4) 
 
clc
 n=9;
P=eye(n);
mu=zeros(n,1);
[xcut6,wcut6]=conjugate_dir_gausspts_till_6moment_scheme2(mu,P);
length(wcut6) 

clc
 n=5;
P=eye(n);
mu=zeros(n,1);
[xcut8,wcut8]=conjugate_dir_gausspts_till_8moment(mu,P);
length(wcut8) 

 clc
 n=5;
[xngkpn wngkpn] = nwSpGr('KPN', n, 5);  % KPN   GQN
[sum(abs(wng)),length(wng)]

 clc
 n=5;
[xngqn wngqn] = nwSpGr('GQN', n, 5);  % KPN   GQN
[sum(abs(wng)),length(wng)]


min(wng)
sum(abs(wng))
max(wng)
mgh3=cal_moments_wrt_pts(xng, wng,8);
mgh3(:,end)

[xsp3,ws]=smolyak_sparse_grid_modf(mu,P,6,3,'gh');
length(ws)
min(ws)
sum(abs(ws))
max(ws)


 [xghT,wghT]=GH_points(zeros(5,1),eye(5),7);

mgh5=cal_moments_wrt_pts(xgh,wgh,12);
mkpn=cal_moments_wrt_pts(xngkpn,wngkpn,12);
mgqn=cal_moments_wrt_pts(xngqn,wngqn,12);
mcut8=cal_moments_wrt_pts(xcut8,wcut8,12);

mghT=cal_moments_wrt_pts(xghT,wghT,12);


max(round(abs([mcut8(:,end)-mghT(:,end),mgh5(:,end)-mghT(:,end),mkpn(:,end)-mghT(:,end),mgqn(:,end)-mghT(:,end)])),[],1)

%%
n=9;
% [xcheng5,wcheng5]=high_cub_chengpaer(zeros(n,1),eye(n));
[xcheng7,wcheng7]=high_cub_cheng_paper(zeros(n,1),eye(n),7);
[sum(abs(wcheng7)),length(wcheng7)]