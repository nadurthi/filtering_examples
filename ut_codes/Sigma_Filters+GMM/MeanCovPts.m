function [mu,P]=MeanCovPts(w,X)

[y1,mu]=Cal_moments_samples(X,w,1,'central');
[y2,M2]=Cal_moments_samples(X,w,2,'central');

mu=mu(:);
P=zeros(length(mu),length(mu));
for i=1:1:size(y2,1)
    if max(y2(i,:))==2
        ind=find(y2(i,:)==2);
        P(ind,ind)=M2(i);
    else
        ind=find(y2(i,:)==1);
        P(ind(1),ind(2))=M2(i);
        P(ind(2),ind(1))=M2(i);
    end
    
end

end