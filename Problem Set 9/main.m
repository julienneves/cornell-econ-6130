%% New Keynesian Model - Simulated Annealing
% Based off Tilahun Emiru's and Kris Nimark's codes.
% By Julien Neves

%% Housekeeping
close all;
clear all;
warning off all;

%% Part (1)
global Z
% Set handle for FRED data
url = 'https://fred.stlouisfed.org/';
c = fred(url);

% Set dates for sample period
startdate = '01/01/1983';
enddate = '12/01/2007';

CPI = fetch(c,'USACPIALLQINMEI','09/01/1982',enddate); % fetch CPI from FRED
GDP = fetch(c,'GDPC1',startdate,enddate); % fetch GDP from FRED
INT = fetch(c,'IRSTFR01USQ156N',startdate,enddate); % fetch rate from FRED
UNR = fetch(c,'LRUN64TTUSQ156S',startdate,enddate); % fetch rate from FRED

pi_data = diff(log(CPI.Data(:,2))); % log[price(t)] - log[price(t-1)]
gdp_data = log(GDP.Data(:,2)); % log of GDP
i_data = log(1+INT.Data(:,2)/100); % log of nominal interest rate
n_data = log(1-UNR.Data(:,2)/100); % log of employment rate

[~, pi_data] = hpfilter(pi_data,1600); % extract cyclical component of pi
[~,gdp_data] = hpfilter(gdp_data,1600); % extract cyclical component of y
[~,i_data] = hpfilter(i_data,1600); % extract cyclical component of i
[~,n_data] = hpfilter(n_data,1600); % extract cyclical component of n

% Combine data
Z = [pi_data';gdp_data'; i_data'; n_data'];

% Plot FRED detrended data
figure(5)
subplot(2,2,1); plot(INT.Data(:,1),Z(1,:)); datetick('x','yyyy'); title('Inflation');
subplot(2,2,2); plot(INT.Data(:,1),Z(2,:)); datetick('x','yyyy'); title('Output');
subplot(2,2,3); plot(INT.Data(:,1),Z(3,:)); datetick('x','yyyy'); title('Nominal Interest Rates');
subplot(2,2,4); plot(INT.Data(:,1),Z(4,:)); datetick('x','yyyy'); title('Labor');
saveas(gcf, 'data.png')

%% Part (2)
% Calibration
sigma     = 2; % CRRA parameter.
beta      = 0.99; % discount factor
phi       = 3; % inverse of elastiicity of labor supply
eps       = 5; % elasticity of substitution between goods i and j
phi_pi    = 1.5; % taylor rule parameter
phi_y     = 0.5; % taylor rule parameter
theta     = 0.75; % degree of price stickiness
alpha     = 0.3; % production function parameter
rho_v     = 0.5; % persistence parameter
rho_a     = 0.8; % persistence parameter
rho_z     = 0.7; % persistence parameter
rho_u     = 0.5; % persistence parameter
sigma_v   = 0.01; % standard deviation
sigma_a   = 0.008; % standard deviation
sigma_z   = 0.03; % standard deviation
sigma_u   = 0.01; % standard deviation

%% Part (5)
% Set starting value
THETA = [sigma; beta; phi; eps; phi_pi; phi_y; theta; alpha;
    rho_v; rho_a; rho_z; rho_u; sigma_v; sigma_a; sigma_z; sigma_u];

LB= [0  0 1  1  1 0 0 0 0 0 0 0 0 0 0 0]'; % lower bound
UB= [10 1 10 25 5 5 1 1 1 1 1 1 10 10 10 10]'; % upper bound

% Simulated Annealing - Matlab
options = optimoptions(@simulannealbnd, 'Display', 'diagnose','PlotFcn', @saplotx, 'MaxTime', 600, 'FunctionTolerance', 1e-2, 'MaxStallIterations',1000);
logL = @loglikelihood_DSGE;
xhat = simulannealbnd(logL, THETA, LB, UB, options);

% % Simulated Annealing - Class
% sa_t= 5;
% sa_rt=.3;
% sa_nt=5;
% sa_ns=5;
% xhat = simannb( 'loglikelihood_DSGE', THETA, LB, UB, sa_t, sa_rt, sa_nt, sa_ns, 1);

thetalabel = ['sigma   '; 'beta    '; 'phi     '; 'eps     '; 'phi_pi  '; 'phi_y   '; 'theta   '; 'alpha   ';
    'rho_v   '; 'rho_a   '; 'rho_z   '; 'rho_u   '; 'sigma_v '; 'sigma_a '; 'sigma_z '; 'sigma_u '];
disp('ML estimate of THETA')
disp([thetalabel, num2str(xhat)])

%% Part (6)

% Compute Neq Keynesian Model with estimated theta
[A, C, D ,~, T] = nkbc_model(xhat, 'impulse');

% Compute impulse responses
time = 0:10;   % set time horizon
col = T;   % start impulse matrix
for j=1:length(time)
    resp(:,:,j)=D*col; % compute observations
    col=A*col;  % compute next period states
end

for i = 1:4
    % Extract impulse responses for observations
    resp_pi(:,i)=squeeze(resp(1,i,:));
    resp_y(:,i)=squeeze(resp(2,i,:));
    resp_yg(:,i)=squeeze(resp(3,i,:));
    resp_i(:,i)=squeeze(resp(4,i,:));
    resp_n(:,i)=squeeze(resp(5,i,:));
    
    % Plot Impulse Responses
    figure(i)
    subplot(2,3,1); plot(time,resp_pi(:,i),'-O'); title('Inflation'); grid on;
    subplot(2,3,2); plot(time,resp_y(:,i),'-O'); title('Output'); grid on;
    subplot(2,3,3); plot(time,resp_yg(:,i),'-O'); title('Output Gap'); grid on;
    subplot(2,3,4); plot(time,resp_i(:,i),'-O'); title('Nominal Interest Rate'); grid on;
    subplot(2,3,5); plot(time,resp_n(:,i),'-O'); title('Labor'); grid on;
end

% Print Impulse Responses - Monetary Shock
figure(1)
saveas(gcf, 'impulse_monetary.png')
% Print Impulse Responses - Productivity Shock
figure(2)
saveas(gcf, 'impulse_prod.png')
% Print Impulse Responses - Demand Shock
figure(3)
saveas(gcf, 'impulse_demand.png')
% Print Impulse Responses - Cost-push Shock
figure(4)
saveas(gcf, 'impulse_cost.png')

%% Part (7)
% Solve New Keynesian Model
[A,C,D] = nkbc_model( xhat, 'data');

% Set starting values
X0 = zeros(size(A,2),1);    % set starting X
P0 = dlyap(A,C*C'); % set starting for the variance

% Compute the Kalman Filter
[ X_post, P_post, X_prior, Z_tilde, Omega] = kfilter(Z, A, C, D, 0, X0, P0 );

figure(6)
plot(Z(2,:));
hold
plot(X_post(1,:),'--'); 
legend('Output','Output Gap', 'Location', 'best'); grid on;
saveas(gcf, 'kalman.png')
