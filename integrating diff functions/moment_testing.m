%checking the rules 
%dim
clear
clc
N=6;
%Moment
M=8;
% [Xcut,wcut]=uniform_sigma_pts(-1*ones(1,N),ones(1,N),M);
[Xcut,wcut]=smolyak_sparse_grid(N,(M+2)/2,'GLgn');
[Xgl,wgl] = GLeg_pts((M+2)/2*ones(1,N), -1*ones(1,N), ones(1,N));
%calculating all the moments
A=[];
if N==2
A=vertcat(A,[sum(wcut.*Xcut(:,1).^2),sum(wgl.*Xgl(:,1).^2)]);
A=vertcat(A,[sum(wcut.*Xcut(:,1).^4),sum(wgl.*Xgl(:,1).^4)]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^2.*Xcut(:,2).^2)),sum(wgl.*(Xgl(:,1).^2.*Xgl(:,2).^2))]);
A=vertcat(A,[sum(wcut.*Xcut(:,1).^6),sum(wgl.*Xgl(:,1).^6)]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^4.*Xcut(:,2).^2)),sum(wgl.*(Xgl(:,1).^4.*Xgl(:,2).^2))]);
A=vertcat(A,[sum(wcut.*Xcut(:,1).^8),sum(wgl.*Xgl(:,1).^8)]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^6.*Xcut(:,2).^2)),sum(wgl.*(Xgl(:,1).^6.*Xgl(:,2).^2))]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^4.*Xcut(:,2).^4)),sum(wgl.*(Xgl(:,1).^4.*Xgl(:,2).^4))]);
end
if N==3
A=vertcat(A,[sum(wcut.*Xcut(:,1).^2),sum(wgl.*Xgl(:,1).^2)]);
A=vertcat(A,[sum(wcut.*Xcut(:,1).^4),sum(wgl.*Xgl(:,1).^4)]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^2.*Xcut(:,2).^2)),sum(wgl.*(Xgl(:,1).^2.*Xgl(:,2).^2))]);
A=vertcat(A,[sum(wcut.*Xcut(:,1).^6),sum(wgl.*Xgl(:,1).^6)]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^4.*Xcut(:,2).^2)),sum(wgl.*(Xgl(:,1).^4.*Xgl(:,2).^2))]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^2.*Xcut(:,2).^2.*Xcut(:,3).^2)),sum(wgl.*(Xgl(:,1).^2.*Xgl(:,2).^2.*Xgl(:,3).^2))]);
A=vertcat(A,[sum(wcut.*Xcut(:,1).^8),sum(wgl.*Xgl(:,1).^8)]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^6.*Xcut(:,2).^2)),sum(wgl.*(Xgl(:,1).^6.*Xgl(:,2).^2))]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^4.*Xcut(:,2).^4)),sum(wgl.*(Xgl(:,1).^4.*Xgl(:,2).^4))]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^4.*Xcut(:,2).^2.*Xcut(:,3).^2)),sum(wgl.*(Xgl(:,1).^4.*Xgl(:,2).^2.*Xgl(:,3).^2))]);
end
if N>3
A=vertcat(A,[sum(wcut.*Xcut(:,1).^2),sum(wgl.*Xgl(:,1).^2)]);
A=vertcat(A,[sum(wcut.*Xcut(:,1).^4),sum(wgl.*Xgl(:,1).^4)]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^2.*Xcut(:,2).^2)),sum(wgl.*(Xgl(:,1).^2.*Xgl(:,2).^2))]);
A=vertcat(A,[sum(wcut.*Xcut(:,1).^6),sum(wgl.*Xgl(:,1).^6)]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^4.*Xcut(:,2).^2)),sum(wgl.*(Xgl(:,1).^4.*Xgl(:,2).^2))]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^2.*Xcut(:,2).^2.*Xcut(:,3).^2)),sum(wgl.*(Xgl(:,1).^2.*Xgl(:,2).^2.*Xgl(:,3).^2))]);
A=vertcat(A,[sum(wcut.*Xcut(:,1).^8),sum(wgl.*Xgl(:,1).^8)]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^6.*Xcut(:,2).^2)),sum(wgl.*(Xgl(:,1).^6.*Xgl(:,2).^2))]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^4.*Xcut(:,2).^4)),sum(wgl.*(Xgl(:,1).^4.*Xgl(:,2).^4))]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^4.*Xcut(:,2).^2.*Xcut(:,3).^2)),sum(wgl.*(Xgl(:,1).^4.*Xgl(:,2).^2.*Xgl(:,3).^2))]);
A=vertcat(A,[sum(wcut.*(Xcut(:,1).^2.*Xcut(:,2).^2.*Xcut(:,3).^2.*Xcut(:,4).^2)),sum(wgl.*(Xgl(:,1).^2.*Xgl(:,2).^2.*Xgl(:,3).^2.*Xgl(:,4).^2))]);
end
disp('dim=')
N
disp('MOM=')
M
A