function doit()
X=[];
for nt=1:1:4001
for i=1:1:9
    load(strcat('COND_1_MC_runs_no_',num2str(i),'__5000'))
    X=vertcat(X,X_mc(nt,:,:));

%     MU(:,:,i)=mu;
%     PP(:,:,i)=P;
end
 w=ones(size(X_mc,3),1)/size(X_mc,3);
 [mu,P]=cal_moms_wrt_allstates(X_mc,w,1,'raw');
save('MOMS','MU','PP')
end
