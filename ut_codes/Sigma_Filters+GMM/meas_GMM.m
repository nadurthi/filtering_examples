function GMMz=meas_GMM(prior_GMM,nz,h,g,R,s)

GMMz.w=prior_GMM.w;

ng=length(prior_GMM.w);
nx=size(prior_GMM.mu,2);

GMM.mu=zeros(size(prior_GMM.mu));
GMM.P=zeros(size(prior_GMM.P));
GMMz.mu=zeros(ng,nz);
GMMz.P=zeros(ng,nz^2);
Pxz=zeros(ng,nx*nz);

[X,W]=GH_pts(zeros(nx,1),eye(nx),2);

for i=1:1:ng
    mu=prior_GMM.mu(i,:)';
    sqP=sqrtm(reshape(prior_GMM.P(i,:),nx,nx));
    
    for j=1:1:length(W)
       Xj=sqP*X(j,:)'+mu;
       hx=h(Xj,s);
       GG=g(Xj,s)*R*g(Xj,s)';

       GMMz.mu(i,:)=GMMz.mu(i,:)+W(j)*hx';
       GMMz.P(i,:)=GMMz.P(i,:)+W(j)*reshape(GG+hx*hx',1,nz^2);
       
    end
    GMMz.P(i,:)=reshape(reshape(GMMz.P(i,:),nz,nz)-GMMz.mu(i,:)'*GMMz.mu(i,:),1,nz^2);

end
