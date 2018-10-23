function [x,g]=permute_moments(P,N)
n=length(P);
if N==2
   g=[];
   for i=1:1:n
       g=horzcat(g,P(i,i:1:n));
   end
  g=g'; 
   return
end

% N is the moment(always keep it even)
% n is dimension of system
n=length(P);
combos = combntns(1:N,2);
ind=combntns(1:length(combos),N/2);
A=[];
a=[1:1:N];
for i=1:1:length(ind)
    r=[];
    for j=1:1:N/2
        r=horzcat(r,combos(ind(i,j),:));
    end
    if sum(abs(sort(r)-a))==0
    A=vertcat(A,r) ;
    end   
end

% combos = combntns(0:N,n)
combos = GenerateIndex(n,(N+1)*ones(1,n));
combos(find(combos==(N+1)))=0;

x=[];
for i=1:1:length(combos)
    if sum(combos(i,:))==N
     x=vertcat(x,combos(i,:));
%      x=vertcat(x,wrev(combos(i,:)));
    end
end
size(x);
% x=vertcat(x,N/2*ones(1,n));
x=sortrows(x,-1);
nn=size(x);
% m=zeros(nn(1),1);
g=[];
m=0;
for i=1:1:nn(1)
    
    p=[];
    for h=1:1:n
        p=horzcat(p,h*ones(1,x(i,h)));
    end
    C=A;
    for h=1:1:N
        C(find(A==h))=p(h);
    end
    nn=size(A);
    for j=1:1:nn(1);
        cc=1;
        for k=1:2:N-1
            cc=cc*P(C(j,k),C(j,k+1));
        end
        m=m+cc;
    end
    g=vertcat(g,m);
    m=0;
end
%  g=[x,g];   
% if N==4
%     m=[3*P(1,1)^2;3*P(1,1)*P(1,2);P(1,1)*P(2,2)+2*P(1,2)^2;3*P(2,2)*P(2,1);3*P(2,2)^2];
% end
% EE=@(x1,x2,x3,x4,x5,x6)(P(x1,x2)*P(x3,x4)*P(x5,x6)+P(x1,x2)*P(x3,x5)*P(x4,x6)+P(x1,x2)*P(x3,x6)*P(x4,x5)...
%     +P(x1,x3)*P(x2,x4)*P(x5,x6)+P(x1,x3)*P(x2,x5)*P(x4,x6)+P(x1,x3)*P(x2,x6)*P(x4,x5)...
%     +P(x1,x4)*P(x2,x3)*P(x5,x6)+P(x1,x4)*P(x2,x5)*P(x3,x6)+P(x1,x4)*P(x2,x6)*P(x3,x5)...
%     +P(x1,x5)*P(x2,x3)*P(x4,x6)+P(x1,x5)*P(x2,x4)*P(x3,x6)+P(x1,x5)*P(x2,x6)*P(x3,x4)...
%     +P(x1,x6)*P(x2,x3)*P(x4,x5)+P(x1,x6)*P(x2,x4)*P(x3,x5)+P(x1,x6)*P(x2,x5)*P(x3,x4));
% if N==6
%     m=[EE(1,1,1,1,1,1),EE(1,1,1,1,1,2),EE(1,1,1,1,2,2),EE(1,1,1,2,2,2),EE(1,1,2,2,2,2)...
%         ,EE(1,2,2,2,2,2),EE(2,2,2,2,2,2)];
%     m=m';
% end
end