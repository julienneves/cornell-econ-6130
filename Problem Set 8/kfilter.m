function [ X_post, P_post, X_prior, Z_tilde, Omega] = kfilter(Z, A, C, D, S_vv, X0, P0 )
%KFILTER Compute the Kalman Filter for a VAR(1) process
% 	X[t+1] = AX[t] + Cu[t+1]
%   Z[t+1] = DX[t+1] + v[t+1]

% Get size of observations
T = length(Z);

% Allocate space for estimates
X_prior = cell(1,T);
X_post = cell(1,T+1);
K = cell(1,T);
Z_tilde = cell(1,T);
Omega = cell(1,T);
P_post = cell(1,T+1);

% Set starting values
X_post{1} = X0;
P_post{1} = P0;

for t = 1:T
    % Compute X_{t+1|t}
    X_prior{t+1}  = A*X_post{t} ;
    
    % Compute P_{t+1|t}
    P_prior = A * P_post{t} * A' + C * C';
    
    % Compute K_{t+1}
    K{t+1} = P_prior*D'*(D*P_prior*D'+S_vv)^-1;
    
    % Innovation
    Z_tilde{t} = Z(:,t) - D*X_prior{t+1};
    Omega{t} = D*P_prior*D'+S_vv ;
    
    % Compute X_{t+1|t+1}
    X_post{t+1}  = A*X_post{t}  + K{t+1}  * Z_tilde{t};
       
    % Compute P_{t+1|t+1}
    P_post{t+1} = P_prior- K{t+1}*D*P_prior;
end

% Convert post and prior estimates to matrices
X_prior = cell2mat(X_prior);
X_post = cell2mat(X_post);
end