%% Value Function Iteration
clear, clc;

%% Set parameters
nk = 500;	% Size of grid
ns = 7;
dev = 3;

delta = 0.1;	% Depreciation rate
alpha = 0.3;    % Capital share of income
beta = 0.96;     % Discount factor

phi = 0.98; % AR coefficient
sig_y = sqrt(.1); % Standard deviation of y_t

%% Markov (Tauchen)
N = 10000;  % Number of simulation

sample = ones(N,3); % Initialize mean, std, autocorr. placeholder

[prob, state_grid] = tauchen(ns, phi, sig_y ,dev); % Transition matrix
for i = 1:N 
chain = markovchain(prob, 1000, 4); % Generate Markov chain
chain = state_grid(chain); % Map markov chain to state values

acf = autocorr(chain,1); % Compute autocorrelation

sample(i,:) = [mean(chain),std(chain), acf(2)]; % Store values
end

mean(sample) % Compute average of mean, std, autocorr.

%% Value function
k_max = 5;	% Upper bound
k_min = 1;	% Lower bound
k_grid = linspace(k_min,k_max,nk)';	% Create grid

% Compute value/policy function
[val_fun, pol_fun] = bellman( k_grid, state_grid, alpha, beta, delta, prob);

%% Data
url = 'https://fred.stlouisfed.org/';
c = fred(url);

GDP = fetch(c,'GDPC1'); % Fetch GDP from FRED
CON = fetch(c,'PCECC96');   % Fetch consumption from FRED
INV = fetch(c,'GPDIC1');    % Fetch investment from FRED

gdp_data = log(GDP.Data(:,2));  % Log of GDP
c_data = log(CON.Data(:,2));    % Log of consumption
i_data = log(INV.Data(:,2));    % Log of investment

[~,gdp_data] = hpfilter(gdp_data,1600); % Extract cyclical component of GDP
[~,i_data] = hpfilter(i_data,1600); % Extract cyclical component of C
[~,c_data] = hpfilter(c_data,1600); % Extract cyclical component of I

%% Simulation
T = size(c_data,1)+1;   % Set size of Markov Chain to match data

y_sim = markovchain(prob, T, 4);    % Generate Markov chain

k_sim = ones(T,1);  % Initialize capital vector
k_sim(1) = round(nk/2); % Set starting value to the middle of grid

for t = 2:T
    k_sim(t)= pol_fun(k_sim(t-1),y_sim(t-1));   % Compute k' from k
end

y_sim = state_grid(y_sim); % Map markov chain to state values
k_sim = k_grid(k_sim);  % Map markov chain to capital values

gdp_sim = exp(y_sim).* k_sim.^alpha;    % Compute GDP
gdp_sim = gdp_sim(1:end-1); % Drop last value

i_sim= k_sim(2:end) -(1-delta)*k_sim(1:end-1);  % Compute investment

c_sim = gdp_sim - i_sim;    % Compute consumption

c_sim = log(c_sim); % Log of c
i_sim = log(i_sim); % Log of i
gdp_sim = log(gdp_sim); % Log of GDP

c_sim = c_sim-mean(c_sim);  % Detrend c
i_sim = i_sim- mean(i_sim); % Detrend i
gdp_sim = gdp_sim - mean(gdp_sim);  % Detrend GDP

std_sim = std([c_sim,gdp_sim,i_sim]) % Compute standard value of simulation
std_data = std([c_data,gdp_data,i_data]) % Compute standard value of data
corr(c_sim,c_data)  % Compute the correlation for consumption
corr(gdp_sim,gdp_data)  % Compute the correlation for GDP
corr(i_sim,i_data)  % Compute the correlation for investment

%% Plot figures
% Value function
figure
plot(k_grid,val_fun) 

xlabel('k')
ylabel('v(k)')
legend('-3\sigma','-2\sigma','-1\sigma','Policy function','+1\sigma','+2\sigma','+3\sigma', 'Location', 'best')

print('plot_value','-dpng')

% Policy function
figure
plot(k_grid,k_grid(pol_fun))
axis equal
hold on
plot(k_grid,k_grid, '--b')

xlabel('k')
ylabel('g(k)')
legend('-3\sigma','-2\sigma','-1\sigma','Policy function','+1\sigma','+2\sigma','+3\sigma','45 degree', 'Location', 'best')

print('plot_policy','-dpng')

% Markov
figure
subplot(1,3,1); histogram(sample(:,1)); title('Mean');
subplot(1,3,2); histogram(sample(:,2)); title('Standard Deviation');
subplot(1,3,3); histogram(sample(:,3)); title('Autocorrelation');

print('plot_markov','-dpng')

% Simulation
figure
plot(GDP.Data(:,1),[c_sim,gdp_sim,i_sim])
hold on
plot(GDP.Data(:,1),[c_data,gdp_data,i_data], '--')

datetick('x')
legend('C - Simulation', 'GDP - Simulation','I - Simulation','C - Data','GDP - Data','I - Data', 'Location', 'best')
print('plot_sim','-dpng')
