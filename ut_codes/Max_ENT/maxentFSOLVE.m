function [ceq,jac]=maxentFSOLVE(lam,y,M,X,W)



nm=size(y,1);
nq=length(W);
ceq=zeros(nm,nq);
jac=zeros(nm^2,nq);
for i=1:length(W)
    ceq(:,i)=W(i)*prod(repmat(X(i,:),nm,1).^y,2)*pdf_MaxEnt(X(i,:),lam,y);
    jac(:,i)=reshape(W(i)*(prod(repmat(X(i,:),nm,1).^y,2)*prod(repmat(X(i,:),nm,1).^y,2)')*pdf_MaxEnt(X(i,:),lam,y),nm^2,1);
end
ceq=sum(ceq,2)-M;
jac=reshape(sum(jac,2),nm,nm);
% [xx,zz]=meshgrid(linspace(-1,1,100),linspace(-1,1,100));
% pent=zeros(size(xx));
% for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%  pent(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam,y);
%     end
%  end
%    contour(xx,zz,pent,linspace(0.0001,max(max(pent,[],1))/2,15))
%  axis([-1,1,-1,1])
%  saveas(gcf,num2str(pic_cnt),'jpg')
%  pic_cnt=pic_cnt+1;
%  pause(0.2)
%  close
end