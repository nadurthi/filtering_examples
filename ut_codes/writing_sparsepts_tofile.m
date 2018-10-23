fid = fopen('sparse_points_uniform.m', 'wt');
fprintf(fid, '%s\n', 'function XW=sparse_points_uniform(n,k)');
fprintf(fid, '%s\n', '%file with all the points');
%% writing sparse points for uniform distribution
fprintf(fid, '%s\n', '%sparse dim 2 to 10; k 2 to 10');

for n=2:1:9
    fprintf(fid, '%s\n', strcat('if n==',num2str(n)));
    for k=2:1:8
        [n,k]
    fprintf(fid, '%s\n', strcat('if k==',num2str(k)));
    [X,w]=smolyak_sparse_grid(n,k,'GLgn');
    
    fprintf(fid, '%s', 'XW=[');
        for i=1:1:length(w)
        fprintf(fid, '% 5.15f ',[X(i,:),w(i)]);
        fprintf(fid, '%s\n', ';');
        end
    fprintf(fid, '%s\n', '];');
    fprintf(fid, '\n%s\n', 'end');
    end
fprintf(fid, '%s\n', 'end');%%% ending cut4
end
fprintf(fid, '%s\n', 'end');%%% ending program
fclose(fid)
