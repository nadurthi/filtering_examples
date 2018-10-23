function [mu,P]=GMM_mean_cov(GMM)
nx=size(GMM.mu,2);
ng=size(GMM.mu,1);
mu=zeros(nx,1);
P=zeros(nx,nx);
 
    for i=1:1:ng
        mu=mu+GMM.w(i)*GMM.mu(i,:)';
        P=P+GMM.w(i)*(reshape(GMM.P(i,:),nx,nx)+GMM.mu(i,:)'*GMM.mu(i,:));
    end
    P=P-mu*mu';
end