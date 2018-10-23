function [mu,P]=cal_moms_wrt_allstates(X,w,N,type)
% given X(time,states,samples) and constant w-weight calculate the moments
% with respect to time for all states
% N is the the order of required moment
nt=size(X,1);%no. of time steps
nx=size(X,2);%no. of states
ns=size(X,3);%no. of samples
if size(w,1)==1
    w=w';
end
% calculating the mean

% switch between central and mean by substracting the mean from the samples
% if strcmp(type,'central')==1
%    X=X-mu; 
% end
P=zeros(nt,nx^2);
mu=zeros(nt,nx);
for i=1:1:nt
    for j=1:1:ns
        x(j,:)=X(i,:,j);
    end
    W=repmat(w,1,nx);
    mk=sum(W.*x,1)';
    MU=repmat(mk',ns,1);
    xx=x-MU; 
    Pk=xx'*(W.*xx);
    mu(i,:)=mk;
    P(i,:)=reshape(Pk,1,nx^2);
end
end      
    