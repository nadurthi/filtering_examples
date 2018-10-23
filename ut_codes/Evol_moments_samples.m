function [y,M]=Evol_moments_samples(X,w,N,type)
% global Xdat
% X=Xdat;
% X=XXmc;
% w=wmc;
% given X(time,states,samples) and constant w-weight calculate the moments
% with respect to time for all states
% N is the the order of required moment
nt=size(X,1);%no. of time steps
nx=size(X,2);%no. of states
ns=size(X,3);%no. of samples
if size(w,1)==1
    w=w';
end
mu=zeros(nt,nx);
for i=1:1:nt
    for j=1:1:nx
        D(1:1:ns,1)=X(i,j,:);
mu(i,j)=sum(w.*D);
    end
end
if N==1
    M=mu;
    y=eye(nx);
    return;
end
if strcmp(type,'central')==1
  disp('central moms')
    for i=1:1:nt
        for j=1:1:nx
            X(i,j,:)=X(i,j,:)-mu(i,j);
        end
    end
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
% if ns<10000
M=zeros(nt,yy);
   for t=1:1:nt
      D=zeros(ns,nx);
         for hh=1:1:nx
        D(1:1:ns,hh)=X(t,hh,1:1:ns);
         end
         
        for i=1:1:yy
        M(t,i)=sum(w.*prod(D.^repmat(y(i,:),ns,1),2));
        end
   end
% else
%    %%%%%%%%%%%%%%%%%%%%%%%%
%    M=zeros(nt,yy);
%    for t=1:1:nt
% %       for i=1:1:yy
% %          for hh=1:1:nx
% %         D(1:1:ns,hh)=X(t,hh,1:1:ns);
% %          end
% XX=0;
%             for ss=0:2000:ns
% %                 for k=1:1:nx
%                 XX(1:2000,1:nx)=X(t,:,1+ss:1:2000+ss);
% %                 end
%             M(t,:)=M(t,:)+reshape(w(ss).*prod(repmat(XX,yy,1).^y,2),1,yy);
%             end
% %       end
%    end
%    
% end