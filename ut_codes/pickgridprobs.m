clc
n=21;
x=-n:1:n;
A=[];
N=2*n;
for i=0:1:N
    A=vertcat(A,0.1^i*x.^i);
end
A
inv(A)