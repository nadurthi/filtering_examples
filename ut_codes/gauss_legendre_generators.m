
[xint,wint] = get_colocation(5*ones(1,4),-1*ones(1,4),1*ones(1,4));
xint=abs(xint);
for i=1:1:length(xint)
if max(xint(i,:))~=0
    xint(i,:)=xint(i,:)/max(xint(i,:));
end
end
xint=sort(xint,2);
x=xint(1,:);
for i=2:1:length(xint)
    [r,c]=size(x);
    cnt=0;
    for j=1:1:r
       
        if sum(abs(xint(i,:)-x(j,:)))==0
          cnt=cnt+1;
        end
        
    end
    if cnt==0
    x=vertcat(x,xint(i,:));
    end
end
xint
xlswrite('8thmom4D',x)