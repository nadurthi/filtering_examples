function D=dist_GM(g,h)
N=length(g.w);
K=length(h.w);
d=zeros(N,K);
for i=1:1:N
    for k=1:1:K
        Pi=reshape(g.P(i,:),sqrt(length(g.P(i,:))),sqrt(length(g.P(i,:))));
        mi=g.mu(i,:)';
%         wi=g.w(i);
        
        Pk=reshape(h.P(k,:),sqrt(length(h.P(k,:))),sqrt(length(h.P(k,:))));
        mk=h.mu(k,:)';
%         wk=h.w(k);
        
        inPi=inv(Pi);
        inPk=inv(Pk);
        C=1/2*(mi'*inPi+mk'*inPk)*inv(inPi+inPk)*(inPi*mi+inPk*mk)-1/2*(mi'*inPi*mi+mk'*inPk*mk);
        mc=inv(inPi+inPk)*(inPi*mi+inPk*mk);
        inPc=(inPi+inPk)/2;
        Pc=inv(inPc);
        BC=sqrt(1/sqrt(det(2*pi*Pi))*1/sqrt(det(2*pi*Pk)))*exp(C/2)*sqrt(det(2*pi*Pc));
        d(i,k)=abs(sqrt(1-BC));
    end
end
A=[];
B=[];
for i=1:1:K
    A=vertcat(A,[zeros(1,(i-1)*N),ones(1,N),zeros(1,N*K-N-(i-1)*N)]);
    B=horzcat(B,[1,zeros(1,N-1)]);
end
A=vertcat(A,B);
for i=1:1:N-1
    B=circshift(B',1)';
    A=vertcat(A,B);   
end
options=optimset('Display','off');
[x,D] = linprog(reshape(d,N*K,1),[],[],A,[h.w;g.w],zeros(N*K,1),[],[],options);
end





