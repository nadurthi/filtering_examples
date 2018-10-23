function p=pdf_MaxEnt(x,lam,y)
if size(x,1)>1
    x=x';
end
if size(lam,1)==1
    lam=lam';
end

ns=size(x,2);
nm=length(lam);

p=exp(sum(prod(repmat(x,nm,1).^y,2).*lam));
end