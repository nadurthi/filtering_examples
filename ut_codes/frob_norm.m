function f=frob_norm(P)
if size(P,1)==1
    n=sqrt(length(P));
    P=reshape(P,n,n);
end
f=sqrt(trace(P*P'));
end
