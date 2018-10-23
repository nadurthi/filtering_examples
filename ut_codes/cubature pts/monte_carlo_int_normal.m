function [x,w]=monte_carlo_int_normal(mu,P,N)
% N is the number of samples
n=length(P);
mu = mu';
R = sqrtm(P);
x = repmat(mu,N,1) + randn(N,n)*R;
w=(1/N)*ones(N,1);
end