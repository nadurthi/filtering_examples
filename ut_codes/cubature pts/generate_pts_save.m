%% gaussian CUT
d=6;
N=8;
mu=zeros(d,1);
P=eye(d);

if N==4
[X,w]=conjugate_dir_gausspts(mu,P);
elseif N==6
    [X,w]=conjugate_dir_gausspts_till_6moment_scheme2(mu,P);
elseif N==8
    [X,w]=conjugate_dir_gausspts_till_8moment(mu,P);
else
    error('lol')
end
    
nme=strcat('cut',num2str(N),'_',num2str(d),'D_gaussian');
save(nme,'X','w')
csvwrite(strcat(nme,'.txt'),[w,X])







%% Uniform CUT

d=2;
N=8;
bdd_low=-1*ones(1,d);
bdd_up=ones(1,d);

[X,w]=uniform_sigma_pts(bdd_low,bdd_up,N);
nme=strcat('cut',num2str(N),'_',num2str(d),'D_uniform');
save(nme,'X','w')
csvwrite(strcat(nme,'.txt'),[w,X])

dlmwrite(strcat(nme,'.txt'), [w,X], 'precision','%.14f')


%% 

%% gaussian CUT
for d=1:7;
N=7;
mu=zeros(d,1);
P=eye(d);


[X,w]=GH_points(mu,P,N);

nme=strcat('GH',num2str(N),'_',num2str(d),'D_gaussian');
save(nme,'X','w')
csvwrite(strcat(nme,'.txt'),[w,X])
end





%% Uniform Gauss legendre


for d=7:8;
N=8;
bdd_low=-1*ones(1,d);
bdd_up=ones(1,d);

[X,w] = GLeg_pts(N*ones(1,d), bdd_low, bdd_up);

nme=strcat('GLgn',num2str(N),'_',num2str(d),'D_uniform');
save(nme,'X','w')
csvwrite(strcat(nme,'.txt'),[w,X])
end




