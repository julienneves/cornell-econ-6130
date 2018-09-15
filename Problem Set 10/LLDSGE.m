function [ logL ] = LLDSGE( THETA )
%loglikelihood_DSGE Likelihood function for State space model
%   State space model:
%   X(t) = A*X(t-1) + C*eps(t)
%   Z(t) = D*X(t)
%   By Julien Neves
global Z

% Solve New Keynesian Model
[A,C,D,eu] = NKBC_model( THETA, 'data');

% Set starting values
X0 = zeros(size(A,2),1);    % set starting X
P0 = dlyap(A,C*C'); % set starting for the variance

% Compute the Kalman Filter
[~, ~, ~, Z_tilde, Omega] = kfilter(Z, A, C, D, 0, X0, P0 );

% Initialize loglikelihood
T = length(Z);
logL = -T/2*log(2*pi)*size(Z,1);

for t = 1:T
    % Update loglikelihood
    logL = logL - 1/2*log(det(Omega{t})) - 1/2* Z_tilde{t}'/Omega{t}* Z_tilde{t};
end
% if imaginary parts or not identified model set likelihood to small value
if (imag(logL)~=0)||(min(eu)==0)
    logL=-9e+200;
end

end

