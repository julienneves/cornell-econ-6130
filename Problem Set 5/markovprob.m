function [Gl, TM, CDF, sz] = markovprob(mue, p, s, N, m)
% markovprob - function
%% Arguments:
% mue = intercept of AR(1) process;
% p = slope coeff. of AR(1) process;
% s = std. dev. of residuals in AR(1) process;
% N = # of grid points for the ‘z’ variable;
% m = Density of the grid for ‘z’ variable;
%% CODE:
sz = s / ((1-p^2)^(1/2)); % Std. Dev. of z.
zmin = -m * sz + mue/(1-p);
zmax = m * sz + mue/(1-p);
z = linspace(zmin,zmax,N); % Grid Points
%% Transition Matrix:
TM = zeros(N,N); % Transition Matrix
w = z(N) - z(N-1);
for j = 1:N;
TM(j,1) = cdf('norm',(z(1)+w/2-mue-p*z(j))/s,0,1);
TM(j,N) = 1 - cdf('norm',(z(N)-w/2-mue-p*z(j))/s,0,1);
for k = 2:N-1;
TM(j,k) = cdf('norm',(z(k)+w/2-mue-p*z(j))/s,0,1)-cdf('norm',(z(k)-w/2-mue-p*z(j))/s,0,1);
end
end
%% Cumulative Distribution Function:
CDF = cumsum(TM,2);
%% Invariant Distribution:
%% Grids:
Gl = exp(z');
fprintf('If we have a lognormal var. (log(z)) in AR(1) process, \n');
fprintf('To make the interval finer at the lower end and coarser at the upper end.\n');