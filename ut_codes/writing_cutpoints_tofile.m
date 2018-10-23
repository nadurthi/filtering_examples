fid = fopen('cut_points_uniform.m', 'wt');
fprintf(fid, '%s\n', 'function XW=cut_points_uniform(n,N)');
fprintf(fid, '%s\n', '%file with all the points');
%% writing cut4 points for gaussian distribution
fprintf(fid, '%s\n', '%CUT4 dim 2 to 10');
fprintf(fid, '%s\n', strcat('if N==',num2str(4)));
for n=2:1:8
[X,w]=uniform_sigma_pts(-1*ones(1,n),ones(1,n),4);
fprintf(fid, '%s\n', strcat('if n==',num2str(n)));
fprintf(fid, '%s', 'XW=[');
for i=1:1:length(w)
fprintf(fid, '% 5.15f ',[X(i,:),w(i)]);
fprintf(fid, '%s\n', ';');
end
fprintf(fid, '%s\n', '];');
fprintf(fid, '\n%s\n', 'end');
end
fprintf(fid, '%s\n', 'end');%%% ending cut4


%%  cut 6 points
%writing cut6 points for gaussian distribution
fprintf(fid, '%s\n', '%CUT6 dim 2 to 9');
fprintf(fid, '%s\n', strcat('if N==',num2str(6)));
for n=2:1:9
[X,w]=uniform_sigma_pts(-1*ones(1,n),ones(1,n),6);
fprintf(fid, '%s\n', strcat('if n==',num2str(n)));
fprintf(fid, '%s', 'XW=[');
for i=1:1:length(w)
fprintf(fid, '% 5.15f ',[X(i,:),w(i)]);
fprintf(fid, '%s\n', ';');
end
fprintf(fid, '%s\n', '];');
fprintf(fid, '\n%s\n', 'end');
end
fprintf(fid, '%s\n', 'end');%%% ending cut6

%% cut8 points
fprintf(fid, '%s\n', '%CUT8 dim 2 to 6');
fprintf(fid, '%s\n', strcat('if N==',num2str(8)));
for n=2:1:5
[X,w]=uniform_sigma_pts(-1*ones(1,n),ones(1,n),8);
fprintf(fid, '%s\n', strcat('if n==',num2str(n)));
fprintf(fid, '%s', 'XW=[');
for i=1:1:length(w)
fprintf(fid, '% 5.15f ',[X(i,:),w(i)]);
fprintf(fid, '%s\n', ';');
end
fprintf(fid, '%s\n', '];');
fprintf(fid, '\n%s\n', 'end');
end
fprintf(fid, '%s\n', 'end');%%% ending cut8
fprintf(fid, '\n%s\n', 'end');% ending the function
fclose(fid);
