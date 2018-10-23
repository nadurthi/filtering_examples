function P=pdf_MaxEnt_extd(x,lam,y,xl,xu,A)
P=zeros(size(x,1),1);
ns=size(x,2);
nm=length(lam);
for i=1:1:size(x,1)
    if sum(sign(x(i,:)-xl))==length(xl) && sum(sign(xu-x(i,:)))==length(xl)
        P(i)=exp(lam'*A*prod(repmat(x,nm,1).^y,2));
    else
        P(i)=0;
        
    end
    
end