
function [val_fun, pol_fun] = bellman( k_grid, state_grid, alpha, beta, delta, prob)
% bellman Value Function Iteration with log utility
%   [val_fun, pol_fun] = bellman( k_grid, state_grid, alpha, beta, delta, prob) returns 
%   val_fun, value function and pol_fun, policy function

crit = 1;	% Initialize convergence criterion
tol = 1e-6;	% Convergence tolerance

%% Grid
nk = size(k_grid,1);
ns = size(state_grid,1);
% Create grid

% Empty
val_temp = zeros(nk,ns);  % Initialize temporary value function vector
val_fun = zeros(nk,ns);	% Initialize value function vector
pol_fun = zeros(nk,ns);	% Initialize policy function vector


%% Value function iteration
while crit>tol ;
    % Iterate on k
    for j = 1:ns
    for i=1:nk   
        c = exp(state_grid(j))* k_grid(i)^alpha + (1-delta)*k_grid(i) - k_grid; % Compute consumption for kt
        utility_c = log(c); % Compute utility for every ct
        utility_c(c<=0) = -Inf;	% Set utility to -Inf for c<=0
        [val_fun(i,j),pol_fun(i,j)] = max(utility_c + beta* val_temp*prob(j,:)');   % Solve Bellman equation
    end
    end
    crit = max(abs(val_fun-val_temp));  % Compute convergence criterion
    val_temp = val_fun;	% Update value function
end

end