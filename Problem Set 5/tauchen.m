function [prob, state_grid] = tauchen(ns, phi, sig_y ,dev)
%markovchain Discretize AR(1) process
%   [prob, cum_prob, state_grid] = tauchen(ns, phi, sig_y ,dev) returns
%   prob, transition matrix,
%   cum_prob, cummulative distribution
%   state_grid, value of states

state_grid = linspace(-dev*sig_y,dev*sig_y,ns)';    % Set state grid
d = mean(diff(state_grid));
% Compute y^j_t+1 - rho y^i_t
state_change = repmat(state_grid,1,ns)'  - rho*repmat(state_grid,1,ns);


prob = zeros(ns,ns);    % Initialize transition matrix
% Compute transition matrix according to Tauchen (1986)
prob(:,1) = cdf('norm',state_change(:,1)+d/2,0,sig_eps);
prob(:,end) = 1- cdf('norm',state_change(:,end) - d/2,0,sig_eps);
prob(:,2:end-1) = cdf('norm',state_change(:,2:end-1)+d/2,0,sig_eps)-cdf('norm',state_change(:,2:end-1)-d/2,0,sig_eps);

end
