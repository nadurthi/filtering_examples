function f=Matlab_prog_submit_Script()
%% Example submisstion script for using the MDCS at The Center for Computational Research, University at Buffalo
%%N. Barlow 2/2010

clear
clc

%%The following lines below should remain unchanged%%
sched = findResource('scheduler', 'type', 'generic');
set(sched, 'ClusterMatlabRoot', '/util/matlab/R2010b');
set(sched, 'HasSharedFilesystem', false); 
set(sched, 'ClusterOsType', 'unix');
set(sched, 'GetJobStateFcn', @getJobStateFcn);
set(sched, 'DestroyJobFcn', @destroyJobFcn);
clusterHost = 'edge.ccr.buffalo.edu';

%%modify below%%
% set(sched, 'DataLocation', '/home/nbarlow2/local'); %linux example
% set(sched, 'DataLocation', '/Users/natebarlow/matlabfiles'); %mac example
set(sched, 'DataLocation', 'C:\Users\nbarlow2\Matlab'); %windows example
remoteDataLocation = '/panasas/scratch/matlab';
ppn=2; % processors per node; ppn=2 (old Dell nodes, Myrinet) or ppn=8 (IBM or new Dell nodes, Infiniband)
time='01:00:00';
email='username@buffalo.edu'; %an e-mail will be sent here when the job completes

%%comment below and set nodeFlag manually if you wish
if ppn<=2
    nodeFlag='GM'; %old Dell nodes, Myrinet
else
    nodeFlag='IB'; %IBM or new Dell nodes, Infiniband
end


%%Don't modify this line%%
set(sched, 'ParallelSubmitFcn', ...
    {@parallelSubmitFcn, clusterHost, remoteDataLocation,ppn,time,nodeFlag,email});

%%Set Maximum Workers (labs)
pjob = createParallelJob(sched,'MinimumNumberOfWorkers',1,'MaximumNumberOfWorkers',4)

%submit a calculation using a built-in matlab function, like rand
createTask(pjob,@rand, 1,{3}); 

%%...or provide your own function 
%set(pjob, 'FileDependencies',{'dosomething.m'}) ;
%createTask(pjob, @dosomething,2); % dosomething.m has 2 outputs

submit(pjob)