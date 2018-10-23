function y=recSIGfun(x,a,b)
i=find(x>=a & x<=b);
y=zeros(size(x));
y(i)=1;
end
