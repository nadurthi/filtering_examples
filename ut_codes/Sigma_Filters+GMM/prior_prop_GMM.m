function GMM=prior_prop_GMM(prior_GMM,F,para,Q)
% F is the dynamics

ng=length(prior_GMM.w);
nx=size(prior_GMM.mu,2);
GMM.mu=zeros(size(prior_GMM.mu));
GMM.P=zeros(size(prior_GMM.P));
GMM.w=zeros(size(prior_GMM.w));
[X,W]=GH_pts(zeros(nx,1),eye(nx),2);
for i=1:1:ng
    mu=prior_GMM.mu(i,:)';
    sqP=sqrtm(reshape(prior_GMM.P(i,:),nx,nx));
    
    for j=1:1:length(W)
       Xj=sqP*X(j,:)'+mu;
       fx=F(Xj,para);
       GMM.mu(i,:)=GMM.mu(i,:)+W(j)*fx';
       GMM.P(i,:)=GMM.P(i,:)+W(j)*reshape(fx*fx',1,nx^2);
    end
    GMM.P(i,:)=reshape(reshape(GMM.P(i,:),nx,nx)-GMM.mu(i,:)'*GMM.mu(i,:)+Q,1,nx^2);
end
GMM.w=prior_GMM.w;
end