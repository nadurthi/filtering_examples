function ef=Expt_GMM(GMM,f)
ng=length(GMM.w);
nx=size(GMM.mu,2);
% calculates the expectation of any function w.r.t a gaussian mixture
[X,W]=GH_pts(zeros(nx,1),eye(nx),2);
H=zeros(ng,1);
for i=1:1:ng
    mu=GMM.mu(i,:)';
    sqP=sqrtm(reshape(GMM.P(i,:),nx,nx));
    
    for j=1:1:length(W)
          Xj=sqP*X(j,:)'+mu;
          H(i)=H(i)+W(j)*real(f(Xj));       
    end
end
ef=H'*GMM.w;
end