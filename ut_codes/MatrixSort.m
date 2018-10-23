function Asorted=MatrixSort(A,dim,ind)
% sort the matrix A according to the dim nad index ind
% i.e. dim=1 .. rows.... ind is the row number
% dime =2 ... col.. ind is column number
if dim==1
    if ind>size(A,2)
        error('index greater than matrix')
    end
v=A(:,ind);

[Vsorted,Indsorted]=sort(v);
Asorted=A(Indsorted,:);

elseif dim==2
    
    if ind>size(A,1)
        error('index greater than matrix')
    end
v=A(ind,:);

[Vsorted,Indsorted]=sort(v);
Asorted=A(:,Indsorted);


end
