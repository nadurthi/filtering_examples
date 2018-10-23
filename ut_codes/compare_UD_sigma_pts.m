% uniform distribution poly integration
clear
clc
f=@(x)sqrt(1+sum(x.^2,2)).^8;
f2=@(x)sqrt(1+sum(x.^2,2)).^9;
f3=@(x)sqrt(1+sum(x.^2,2)).^10;

%% Using Gauss Legendre points
[xgl,wgl] = get_colocation(21*[1 1 1 1], -1*ones(1,4), ones(1,4));
xgl=transform_points_uniform(xgl);
sum_gl=sum(wgl.*f(xgl));
sum_gl2=sum(wgl.*f2(xgl));
sum_gl3=sum(wgl.*f3(xgl));
Ngl=length(wgl);
%% Using CUT8- Uniform Distribution
[xcut8_ud,wcut8_ud] = uniform_sigma_pts(-1*ones(1,4),ones(1,4),8);
% xwcut8_ud=[xcut8_ud,wcut8_ud];
% save(['BentCUTCCUniform' num2str(161) 'pts'],'xwcut8_ud','-ascii','-double');
xcut8_ud=transform_points_uniform(xcut8_ud);
sum_cut8_ud=sum(wcut8_ud.*f(xcut8_ud));
sum_cut8_ud2=sum(wcut8_ud.*f2(xcut8_ud));
sum_cut8_ud3=sum(wcut8_ud.*f3(xcut8_ud));
Ncut8_ud=length(wcut8_ud);
%% Using clenshaw Curtis quadrature
% ND = 4; % Total number of parameters
% m = [-1 -1 -1 -1]; % Kind of distribution for input parameters. uniform for vent radius and grain size while gaussian for vent vel and sigma
% N = 4;
% Np = 2*N+1; % Number of quadrature points along each direction
% [xclen,wclen,pw] = GenerateQuadPoints_clcur(m,ND,Np*ones(1,ND));
load xxclen
load wwclen
sum_clen=sum(wint.*f(xint));
sum_clen2=sum(wint.*f2(xint));
sum_clen3=sum(wint.*f3(xint));
Nclen=length(wint);


A=[sum_gl,100*abs(sum_clen-sum_gl)/sum_gl,100*abs(sum_cut8_ud-sum_gl)/sum_gl;...
    sum_gl2,100*abs(sum_clen2-sum_gl2)/sum_gl2,100*abs(sum_cut8_ud2-sum_gl2)/sum_gl2;...
    sum_gl3,100*abs(sum_clen3-sum_gl3)/sum_gl3,100*abs(sum_cut8_ud3-sum_gl3)/sum_gl3]
[Ngl,Nclen,Ncut8_ud]