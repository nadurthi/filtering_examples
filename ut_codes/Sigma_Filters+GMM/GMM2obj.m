function obj=GMM2obj(GMM)
obj.mu=GMM.mu;
ng=size(GMM.mu,1);
nx=size(GMM.mu,2);

obj.Sigma=[];
for i=1:1:ng
    obj.Sigma=cat(3,obj.Sigma,reshape(GMM.P(i,:),nx,nx));
end
obj.PComponents=GMM.w';

obj = gmdistribution(obj.mu,obj.Sigma,obj.PComponents);
end
