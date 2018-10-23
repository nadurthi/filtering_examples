function f=phi_bistep(a,b,c)
if a<=b && b<c
    f=1;
elseif a>b && b>=c
    f=-1;
else
    f=0;
end