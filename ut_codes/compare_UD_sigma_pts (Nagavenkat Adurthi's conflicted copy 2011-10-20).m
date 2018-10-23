% uniform distribution poly integration
clear
clc
f=@(x)sqrt(1+sum(x.^2,2)*1e-3).^8;


%% Using Gauss Legendre points
[xgl,wgl] = get_colocation(11*[1 1 1 1], -1*ones(1,4), ones(1,4));
xgl=transform_points_uniform(xgl);
sum_gl=sum(wgl.*f(xgl));
Ngl=length(wgl);
%% Using CUT8- Uniform Distribution
[xcut8_ud,wcut8_ud] = uniform_sigma_pts(-1*ones(1,4),ones(1,4),8);
xcut8_ud=transform_points_uniform(xcut8_ud);
sum_cut8_ud=sum(wcut8_ud.*f(xcut8_ud));
Ncut8_ud=length(wcut8_ud);
%% Using clenshaw Curtis quadrature
% ND = 4; % Total number of parameters
% m = [-1 -1 -1 -1]; % Kind of distribution for input parameters. uniform for vent radius and grain size while gaussian for vent vel and sigma
% N = 4;
% Np = 2*N+1; % Number of quadrature points along each direction
% [xclen,wclen,pw] = GenerateQuadPoints_clcur(m,ND,Np*ones(1,ND));
load xclen
load wclen
sum_clen=sum(wclen.*f(xclen));
Nclen=length(wclen);


[sum_gl,sum_clen,sum_cut8_ud]
[Ngl,Nclen,Ncut8_ud]