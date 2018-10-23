function y = reEntryDynamics_mc(~,x)

% process noise
Q=blkdiag(0,0,2.4064e-5,2.4064e-5,0);
beta0 = -0.59783;
R0 = 6374;
Gm0 = 3.986e5;
H0 = 13.406;
%% Compute gravity and drag coefficient
G = computeG(x(1:2),Gm0);
D = computeD(x,beta0,R0,H0);

%% State dynamics
y(1,1) = x(3);
y(2,1) = x(4);
y(3,1) = D*x(3) + G*x(1)+sqrt(2.4064e-5)*randn;
y(4,1) = D*x(4) + G*x(2)+sqrt(2.4064e-5)*randn;
y(5,1) = 0;

function G = computeG(r,Gm0)
R = sqrt(r(1)*r(1) + r(2)*r(2));
G = -Gm0/(R^3);

function D = computeD(x,beta0,R0,H0)
V = sqrt(x(3)*x(3) + x(4)*x(4));
beta = beta0 * exp(x(5));
R = sqrt(x(1)*x(1) + x(2)*x(2));
D = beta*exp((R0 - R)/H0)*V;
