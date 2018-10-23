function GMM=obj2GMM(obj)
GMM.mu=obj.mu;
ng=size(obj.mu,1);
nx=size(obj.mu,2);

GMM.P=zeros(ng,nx^2);
for i=1:1:ng
    GMM.P(i,:)=reshape(obj.Sigma(:,:,i),1,nx^2);
end
GMM.w=obj.PComponents';
end
    
    
    