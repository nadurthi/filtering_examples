function m=cal_all_moms(X,w,N,type)
% type='central' or 'raw'
% calculate all the moments of order N given samples X and corr. weight w
[r,n]=size(X);
% r is the number of samples
% n is the dimension of the system

if N==1
    W=repmat(w,1,n);
    m=sum(W.*X,1);
    return
end

combos = GenerateIndex(n,(N+1)*ones(1,n));
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


m=zeros(1,yy);
switch lower(type)
    case 'raw'
        for i=1:1:yy
        p=w;
           for j=1:1:n
            p=p.*X(:,j).^y(i,j);
           end
          m(i)=sum(p);
        end
    case 'central'
        %first calculate the means
        W=repmat(w,1,n);
        mu=sum(W.*X,1);
        MU=repmat(mu,r,1);
        X=X-MU;
        for i=1:1:yy
        p=w;
        for j=1:1:n
        p=p.*X(:,j).^y(i,j);
        end
         m(i)=sum(p);
        end
    otherwise
        error('are u nuts !!!!only raw or central are allowed')
        
        
end


end