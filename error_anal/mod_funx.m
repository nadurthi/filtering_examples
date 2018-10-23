function y=mod_funx(u,x,a,p)
if x>=a
for i=1:1:length(u)
if u(i)>x
    y(i)=0;
elseif u(i)<a
    y(i)=0;
else
    y(i)=(x-u(i))^p;
end
end
end

if x<a
for i=1:1:length(u)
if u(i)<x
    y(i)=0;
elseif u(i)>a
    y(i)=0;
else
    y(i)=-(x-u(i))^p;
end
end
end

end