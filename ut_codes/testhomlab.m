
global nvar degrees FFUN
nvar=6; degrees=[3 5 5 7 7 7];
FFUN='simfn';
HomStruct=ones(1,7);
LPDstruct=ones(sum(degrees),nvar+1);
lpdsolve
disp('The solutions are:'); disp(dehomog(xsoln,1e-8,HomStruct))
