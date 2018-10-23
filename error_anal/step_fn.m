function f=step_fn(a,b,g)
if a<=b && b<g
    f=1;
elseif g<=b && b<a
    f=-1;
else
    f=0;
end
