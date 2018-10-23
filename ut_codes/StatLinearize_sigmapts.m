function H=StatLinearize_sigmapts(y,x,w)
[N,nx]=size(x);
[N,ny]=size(y);

% A=zeros(nx);
% B=zeros(ny,nx);
% for i=1:1:N
% A=A+w(i)*x(i,:)'*x(i,:);
% B=B+w(i)*y(i,:)'*x(i,:);
% end
% H=B/A;

a=fminunc(@(a)myfun(a,w,y,x,nx,ny,N),reshape(ones(ny,nx),1,ny*nx));
H=reshape(a,ny,nx);
end


function f=myfun(a,w,y,x,nx,ny,N)
A=reshape(a,ny,nx);
f=0;
for i=1:1:N
    f=f+w(i)*norm(y(i,:)'-A*x(i,:)')^2;
end
end