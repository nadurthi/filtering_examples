function [y,M]=Cal_moments_samples(X,w,N,type)
% given X(samples,states) and constant w-weight calculate the moments

% N is the the order of required moment
nx=size(X,2);%no. of states
ns=size(X,1);%no. of samples
if size(w,1)==1
    w=w';
end
W=repmat(w,1,nx);
mu=sum(W.*X,1);
if N==1
    M=mu';
    y=eye(nx);
    return;
end
if strcmp(type,'central')==1
  disp('central moms')
    X=X-repmat(mu,ns,1);

elseif strcmp(type,'raw')==1
     disp('raw moms')    
else
    disp('err')
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
combos = GenerateIndex(nx,(N+1)*ones(1,nx));
combos(find(combos==(N+1)))=0;

y=[];
for i=1:1:length(combos)
    if sum(combos(i,:))==N
     y=vertcat(y,combos(i,:));
%      x=vertcat(x,wrev(combos(i,:)));
    end
end
y=sortrows(y,-1);

[yy,yyy]=size(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=zeros(yy,1);
         
        for i=1:1:yy
        M(i)=sum(w.*prod(X.^repmat(y(i,:),ns,1),2));
        end
