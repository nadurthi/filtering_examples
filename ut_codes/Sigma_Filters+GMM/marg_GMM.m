function nwGMM=marg_GMM(GMM,sts)
sts=sort(sts);
% sts are the states to be kept others are removed
ng=length(GMM.w);
nx=size(GMM.mu,2);
% define a new Gaussian mixture with only sts as first two states in order
nwGMM.mu=zeros(ng,length(sts));
nwGMM.P=zeros(ng,(length(sts))^2);
nwGMM.w=zeros(ng,1);
nwGMM.w=GMM.w;
% keyboard
nwGMM.mu=GMM.mu(:,[sts]);
for gg=1:1:ng
P=reshape(GMM.P(gg,:),nx,nx);
nwGMM.P(gg,:)=reshape(P([sts],[sts]),1,length(sts)^2);    
end
