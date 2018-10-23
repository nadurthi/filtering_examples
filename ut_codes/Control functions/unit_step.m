function g=unit_step(t)
dim=size(t);

if size(t,1)==1
    t=t';
end
g=zeros(size(t));
for i=1:1:length(t)
if t(i)>=0
    g(i)=1;
else
    g(i)=0;
end
end
g=reshape(g,dim);
end