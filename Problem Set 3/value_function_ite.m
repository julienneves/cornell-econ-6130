%% Value Function Iteration
clear, clc;

%% Set parameters
n = 1000;	% Size of grid

delta = 0.8;	% Depreciation rate
alpha = 0.3;    % Capital share of income
beta = 0.6;     % Discount factor

crit = 1;	% Initialize convergence criterion
tol = 1e-10;	% Convergence tolerance

%% Grid
k_star = ((1/beta+delta-1)/alpha)^(1/(alpha-1)); % Steady state
k_max = k_star + 0.8*k_star;	% Upper bound
k_min = k_star - 0.8*k_star;	% Lower bound
k_grid = linspace(k_min,k_max,n)';	% Create grid

% Empty
val_temp = zeros(n,1);  % Initialize temporary value function vector
val_fun = zeros(n,1);	% Initialize value function vector
pol_fun = zeros(n,1);	% Initialize policy function vector

ite = 0;    % Initialize iteration counter

%% Value function iteration
while crit>tol ;
    % Iterate on k
    for i=1:n   
        c = k_grid(i)^alpha + (1-delta)*k_grid(i) - k_grid; % Compute consumption for kt
        utility_c = log(c); % Compute utility for every ct
        utility_c(c<=0) = -Inf;	% Set utility to -Inf for c<=0
        [val_fun(i),pol_fun(i)] = max(utility_c + beta*val_temp);   % Solve bellman equation
    end
    crit = max(abs(val_fun-val_temp));  % Compute convergence criterion
    val_temp = val_fun;	% Update value function
    ite = ite + 1 % Print iteration number
end

%% Plot figures
% Value function
figure
plot(k_grid,val_fun) 

xlabel('k')
ylabel('v(k)')

print('plot_value','-dpng')

% Policy function
figure
plot(k_grid,k_grid(pol_fun))
axis equal
hold on
plot(k_grid,k_grid, '--b')

xlabel('k')
ylabel('g(k)')
legend('Policy Function','45 degree', 'Location', 'best')

print('plot_policy','-dpng')
