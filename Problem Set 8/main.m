%% Housekeeping
close all;
clear all;
% Code based on Tilahun's code for NKBC model.

%% Part (a)
% Calibration
THETA = [2 0.99 3 5 1.5 0.5 0.75 0.3 0.8 0.5 0.7 0.008 0.01 0.03];

% Solve New Keynesian Model
[A,~,C]=gensys(Gamma0,Gamma1,Cons,Psi,Pi);

%% Part (b)
% Set up measurement matrix Z(t)=D*S(t)
D = [0 1 0 0 0 0 0 0 0; % inflation
    1 0 0 0 0 psi_ya 0 0 0; % output
    1 0 0 0 0 0 0 0 0;  % output gap
    0 0 0 1 0 0 0 0 0;  % nominal interest rate
    1/(1-alpha) 0 0 0 0 (psi_ya-1)/(1-alpha) 0 0 0]; % labor

% Compute impulse responses
time = 0:10;   % set time horizon
col = C;    % start impulse matrix
for j=1:size(time,2)
    resp(:,:,j)=D*col; % compute observations
    col=A*col;  % compute next period states
end

for i = 1:3
    % Extract impulse responses for observations
    resp_pi(:,i)=squeeze(resp(1,i,:));
    resp_y(:,i)=squeeze(resp(2,i,:));
    resp_yg(:,i)=squeeze(resp(3,i,:));
    resp_i(:,i)=squeeze(resp(4,i,:));
    resp_n(:,i)=squeeze(resp(5,i,:));
    
    % Plot Impulse Responses
    figure(i)
    subplot(2,3,1)
    plot(time,resp_pi(:,i),'-O')
    title('Inflation')
    subplot(2,3,2)
    plot(time,resp_y(:,i),'-O')
    title('Output')
    subplot(2,3,3)
    plot(time,resp_yg(:,i),'-O')
    title('Output Gap')
    subplot(2,3,4)
    plot(time,resp_i(:,i),'-O')
    title('Nominal Interest Rate')
    subplot(2,3,5)
    plot(time,resp_n(:,i),'-O')
    title('Labor')
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

%% Part (c)
% State space model:
% X(t) = A*X(t-1) + C*eps(t)
% Z(t) = D*X(t)

% Set up measurement matrix Z(t) = D*X(t)
D = [0 1 0 0 0 0 0 0 0;
    1 0 0 0 0 psi_ya 0 0 0;
    0 0 0 1 0 0 0 0 0];

time = 0:200;   % set time horizon
X(:,1)= zeros(9,1); % starting value for states

error = randn(3,size(time,2)); % generate epsilon

% Compute state space model recursively
for t=1:size(time,2)-1
    X(:,t+1)=A*X(:,t)+C*error(:,t+1); % X(t+1) = A*X(t) + C*eps(t+1)
    Z(:,t+1)= D*X(:,t+1);   % Z(t+1) = D*X(t+1)
end

% Compute model standard deviations
sig_model = sqrt(var(Z'))

% Plot simulated state space model
figure(4)
subplot(3,1,1)
plot(time,Z(1,:))
title('Inflation')
subplot(3,1,2)
plot(time,Z(2,:))
title('Output')
subplot(3,1,3)
plot(time,Z(3,:))
title('Nominal Interest Rate')
saveas(gcf, 'simulation.png')

%% Part (d)
% Set handle for FRED data
url = 'https://fred.stlouisfed.org/';
c = fred(url);

% Set dates for sample period
startdate = '01/01/1970';
enddate = '12/01/2016';

CPI = fetch(c,'USACPIALLQINMEI','09/01/1969',enddate); % fetch CPI from FRED
GDP = fetch(c,'GDPC1',startdate,enddate); % fetch GDP from FRED
INT = fetch(c,'IRSTFR01USQ156N',startdate,enddate); % fetch rate from FRED

pi_data = diff(log(CPI.Data(:,2))); % log[price(t)] - log[price(t-1)]
gdp_data = log(GDP.Data(:,2)); % log of GDP
i_data = log(1+INT.Data(:,2)/100); % log of nominal interest rate

[~, pi_data] = hpfilter(pi_data,1600); % extract cyclical component of pi
[~,gdp_data] = hpfilter(gdp_data,1600); % extract cyclical component of y
[~,i_data] = hpfilter(i_data,1600); % extract cyclical component of i

% Combine data
Z=[pi_data';gdp_data'; i_data'];

% Compute model standard deviations
sig_data = sqrt(var(Z'))

% Plot FRED detrended data
figure(5)
subplot(3,1,1)
plot(INT.Data(:,1),Z(1,:))
datetick('x','yyyy')
title('Inflation')
subplot(3,1,2)
plot(INT.Data(:,1),Z(2,:))
datetick('x','yyyy')
title('Output')
subplot(3,1,3)
plot(INT.Data(:,1),Z(3,:))
datetick('x','yyyy')
title('Nominal Interest Rates')
saveas(gcf, 'data.png')

%% Part (e)
% State space model:
% X(t) = A*X(t-1) + C*eps(t)
% Z(t) = D*X(t)

% Set up measurement matrix Z(t) = D*X(t)
D = [0 1 0 0 0 0 0 0 0;
    1 0 0 0 0 psi_ya 0 0 0;
    0 0 0 1 0 0 0 0 0];

% Set starting values
X0 = zeros(9,1);    % set starting X
P0 = dlyap(A,C*C'); % set starting for the variance

% Compute the Kalman Filter
[ X_post, P_post, X_prior, P_prior, K] = kfilter(Z, A, C, D, 0, X0, P0 );

% Plot output and output gap
figure(6)
plot(X_post(1,:))
hold
plot(Z(2,:),'--')
legend('Output Gap', 'Output', 'Location', 'best')
saveas(gcf, 'kalman.png')DPC1',startdate,enddate); % fetch GDP from FRED
INT = fetch(c,'IRSTFR01USQ156N',startdate,enddate); % fetch rate from FRED

pi_data = diff(log(CPI.Data(:,2))); % log[price(t)] - log[price(t-1)]
gdp_data = log(GDP.Data(:,2)); % log of GDP
i_data = log(1+INT.Data(:,2)/100); % log of nominal interest rate

[~, pi_data] = hpfilter(pi_data,1600); % extract cyclical component of pi
[~,gdp_data] = hpfilter(gdp_data,1600); % extract cyclical component of y
[~,i_data] = hpfilter(i_data,1600); % extract cyclical component of i

% Combine data
Z=[pi_data';gdp_data'; i_data'];

% Compute model standard deviations
sig_data = sqrt(var(Z'))

% Plot FRED detrended data
figure(5)
subplot(3,1,1)
plot(INT.Data(:,1),Z(1,:))
datetick('x','yyyy')
title('Inflation')
subplot(3,1,2)
plot(INT.Data(:,1),Z(2,:))
datetick('x','yyyy')
title('Output')
subplot(3,1,3)
plot(INT.Data(:,1),Z(3,:))
datetick('x','yyyy')
title('Nominal Interest Rates')
saveas(gcf, 'data.png')

%% Part (e)
% State space model:
% X(t) = A*X(t-1) + C*eps(t)
% Z(t) = D*X(t)

% Set up measurement matrix Z(t) = D*X(t)
D = [0 1 0 0 0 0 0 0 0;
    1 0 0 0 0 psi_ya 0 0 0;
    0 0 0 1 0 0 0 0 0];

% Set starting values
X0 = zeros(9,1);    % set starting X
P0 = dlyap(A,C*C'); % set starting for the variance

% Compute the Kalman Filter
[ X_post, P_post, X_prior, P_prior, K] = kfilter(Z, A, C, D, 0, X0, P0 );

% Plot output and output gap
figure(6)
plot(X_post(1,:))
hold
plot(Z(2,:),'--')
legend('Output Gap', 'Output', 'Location', 'best')
saveas(gcf, 'kalman.png')