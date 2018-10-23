m=4;
p=2;
q=2;
A=[]
for j=0:1:q-1
    for i=0:1:m-j-1
        A=vertcat(A,[i,j]);
    end
end
for i=0:1:p-1
    for j=q:1:m-i-1
        A=vertcat(A,[i,j]);
    end
end
A
B=[];
for i=0:1:m
    for j=0:1:m
        if i+j<m
            B=vertcat(B,[i,j]);
        end
    end
end
B
