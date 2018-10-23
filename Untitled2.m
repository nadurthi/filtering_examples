     nvar=2;
     LPDdegrees=[2 2];
eop=[0 -1 0];
tableau=[
    1 2 0
    -1 1 0
    -2 0 0
    eop
    1 1 1
    -1 0 0
    eop
    ];
xw=[1 0 1]; yw=[0 1 1];
LPDstruct=[
    1 0 1
    1 0 1
    1 0 1
    0 1 1];
HomStruct=[];
lpdtab
 disp('the sols are');disp(dehomog(xsoln,1e-8))