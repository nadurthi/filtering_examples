%% decoding the axis
% 1-p
% 2M-cM
% 3- scaled -con
%% creating the input file for the bertini program
% dimension
N=5;
% no of mce used from start in order 
k=13;
% the axis used in order of variables used
cut_axis=[1,24,22,24];
nopt_axis=[];
% calculate the number of points in each axis
for i=1:1:length(cut_axis)
   if  cut_axis(i)==1
       nopt_axis(i)=2*N;
   end
   if  floor(cut_axis(i)/10)==2
       M=rem(cut_axis(i),10);
       nopt_axis(i)=2^M*nchoosek(N,M);
   end
   if  cut_axis(i)==3
       nopt_axis(i)=N*2^N;
   end
end

% symbolize the variables
syms r1 r2 r3  r4 w1 w2 w3 w4

%% generating the moment constraint equations

data = importdata('real_finite_solutions',' ',1);
data=data.data;
%%specify the number of variables
nvars=7;
nroots=size(data,1)/nvars;
%% number of roots
sols = zeros(nvars,nroots);

%% separate the roots
j=1;
for i=1:nvars:size(data,1)
    sols(:,j)=data(i:i+nvars-1,1);
    j=j+1;
end
%% removing the unnecessary solutions
rmm=[];
%for cut8-G dim 5 upto 13 mce only
if N==5 && k==13
    for i=1:1:nroots
        if sum(sign(sols(1:6,i)))<6
         rmm=horzcat(rmm,i);
        end
    end
   sols(:,rmm)=[]; 
   sols=sols';
    save('ssols','sols')
end

            