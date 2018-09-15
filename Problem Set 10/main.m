%% New Keynesian Model - Simulated Annealing
% Based of Kris Nimark's code.
% Modified by Julien Neves

%% Housekeeping
close all;
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

% Calibration
sigma     = 2; % CRRA parameter.
beta      = 0.99; % discount factor
phi       = 4; % inverse of elastiicity of labor supply
eps       = 5; % elasticity of substitution between goods i and j
phi_pi    = 1.5; % taylor rule parameter
phi_y     = 0.125; % taylor rule parameter
theta     = 0.75; % degree of price stickiness
alpha     = 0.33; % production function parameter
rho_v     = 0.7; % persistence parameter
rho_a     = 0.95; % persistence parameter
rho_z     = 0.2; % persistence parameter
rho_u     = 0.5; % persistence parameter
sigma_v   = 0.01; % standard deviation
sigma_a   = 0.008; % standard deviation
sigma_z   = 0.008; % standard deviation
sigma_u   = 0.01; % standard deviation

% Set starting value
THETA = [sigma; beta; phi; eps; phi_pi; phi_y; theta; alpha;
    rho_v; rho_a; rho_z; rho_u; sigma_v; sigma_a; sigma_z; sigma_u];

LB= [0  0 1  1  1 0 0 0 0 0 0 0 0 0 0 0]'; % lower bound
UB= [10 1 10 25 5 5 1 1 1 1 1 1 10 10 10 10]'; % upper bound

x=THETA;

%Number of draws
S=5e5;%Number of draws in MCMC
burnin=0.2*S;% Fraction of MCMC disregarded as "burnin sample"
J=1000;%Number of draws from MCMC used to simulate posterior functions of theta
epseye=1e-6; %scale up or down to tune acceptance ratio
adaptive=1; %set to 1 to use adaptive proposal density

%--------------------------------------------------------------------------
% Initializes the MH algorithm
%--------------------------------------------------------------------------

%Initializes the proposal variance.
%Set up so that the first candidate draw is always accepted
lpostdraw = -9e+200;
bdraw=x;%Initial draw
vscale=diag(abs(theta))*1d-4+1e-5*eye(length(x));%Initial covariacne of proposal density (not really important)

thetaMCMC=zeros(length(x),S);%Store all draws in thetaMCMC

%Matrices that keep track of switches and drwas outside LB and UB
OutsideProp=zeros(S,1);%Keep track of proprtion of draws outside parmater bounds
SwitchesProp=zeros(S,1);%Keep track of proprtion of switches (accepted draws)
%Number of draws outside parameter boundaries
q=0;
%Number of switches (acceptances)
pswitch=0;
%Iteration counter
iter=0;

%--------------------------------------------------------------------------
% MH algorithm starts here
%--------------------------------------------------------------------------

tic
for iter=1:S
    % Draw from proposal density Theta*_{t+1} ~ N(Theta_{t},vscale)
    bcan = bdraw + norm_rnd(vscale);
    
    if min(bcan > LB)==1 && min(bcan < UB)==1
        lpostcan = log_prior_DSGE(bcan)+LLDSGE(bcan);
        laccprob = lpostcan-lpostdraw;
    else
        laccprob=-9e+200;%assign very low value to reject
        q=q+1;
    end
    
    %Accept candidate draw with log prob = laccprob, else keep old draw
    if log(rand)<laccprob
        lpostdraw=lpostcan;
        bdraw=bcan;
        pswitch=pswitch+1;
    end
    
    thetaMCMC(:,iter)=bdraw; %add accepted new draw (or the previous value if candidate was rejected) to chain
    
    OutsideProp(iter)=q/iter;
    SwitchesProp(iter)=pswitch/iter;
    
    if mod(iter,10000)==0 && iter > 10000 %do this every 10000 draws if iter > 10000
        disp(['iter: ',num2str(iter)]);
        disp(['acceptance rate: ',num2str(SwitchesProp(iter))]);
        if adaptive==1
            vscale=3d-3*cov(thetaMCMC(:,1000:iter)')+1e-10*eye(16);%update the covariance of proposal density (the second component is just to avoid possible singualrites o
        end
    end
end
toc

thetaMCMC_raw=thetaMCMC;
thetaMCMC=thetaMCMC(:,burnin:end); % disregard burnin sample and keep only evey 100th draw

%% Part (2)
% Plot raw data
plotraw(thetaMCMC_raw,1);
plotraw(thetaMCMC,0);

% Plot cumulative mean
convcheck(thetaMCMC);

%% Part (3)
plotpost(thetaMCMC,0);
close all

%% Impulse resonses / Variance decomposition
time = 0:25;
ra = max(size(thetaMCMC));
FF = ceil(ra.*rand(J,1));

upper=ceil(J*.95);
med=ceil(J*.5);
lower=ceil(J*.05);

% Compute impulse responses
for j = 1:J
theta=thetaMCMC(:,FF(j,1));

% Compute Neq Keynesian Model with estimated theta
[A, C, D ,~, T] = NKBC_model(theta, 'impulse');
col = T;   % start impulse matrix
sigma = D*dlyap(A,C*C')*D';

for s=1:length(time)
    resp(:,:,s)=D*col; % compute observations
    col=A*col;  % compute next period states
end

for i = 1:4
    % Extract impulse responses for observations
    resp_pi(:,i,j)=squeeze(resp(1,i,:));
    resp_y(:,i,j)=squeeze(resp(2,i,:));
    resp_yg(:,i,j)=squeeze(resp(3,i,:));
    resp_i(:,i,j)=squeeze(resp(4,i,:));
    resp_n(:,i,j)=squeeze(resp(5,i,:));
    
    sigma_i = D*dlyap(A,C(:,i)*C(:,i)')*D';
    var_decomp(:,i,j)=diag(sigma_i)./diag(sigma);
end
end

for i = 1:4
    % Extract impulse responses for observations
    resp_pi(:,i,:)=sort(resp_pi(:,i,:),3);
    resp_y(:,i,:)=sort(resp_y(:,i,:),3);
    resp_yg(:,i,:)=sort(resp_yg(:,i,:),3);
    resp_i(:,i,:)=sort(resp_i(:,i,:),3);
    resp_n(:,i,:)=sort(resp_n(:,i,:),3);
    
    % Plot Impulse Responses
    figure(i)
    subplot(2,3,1); 
    plot(time,resp_pi(:,i,upper),':','color','black','LineWidth',2); hold on; 
    plot(time,resp_pi(:,i,med),'color','black','LineWidth',2); hold on; 
    plot(time,resp_pi(:,i,lower),':','color','black','LineWidth',2); hold on; 
    title('Inflation'); grid on;
    subplot(2,3,2); 
    plot(time,resp_y(:,i,upper),':','color','black','LineWidth',2); hold on; 
    plot(time,resp_y(:,i,med),'color','black','LineWidth',2); hold on; 
    plot(time,resp_y(:,i,lower),':','color','black','LineWidth',2); hold on; 
    title('Output'); grid on;
    subplot(2,3,3); 
    plot(time,resp_yg(:,i,upper),':','color','black','LineWidth',2); hold on; 
    plot(time,resp_yg(:,i,med),'color','black','LineWidth',2); hold on; 
    plot(time,resp_yg(:,i,lower),':','color','black','LineWidth',2); hold on; 
    title('Output Gap'); grid on;
    subplot(2,3,4);
    plot(time,resp_i(:,i,upper),':','color','black','LineWidth',2); hold on; 
    plot(time,resp_i(:,i,med),'color','black','LineWidth',2); hold on; 
    plot(time,resp_i(:,i,lower),':','color','black','LineWidth',2); hold on; 
    title('Nominal Interest Rate'); grid on;
    subplot(2,3,5);
    plot(time,resp_n(:,i,upper),':','color','black','LineWidth',2); hold on; 
    plot(time,resp_n(:,i,med),'color','black','LineWidth',2); hold on; 
    plot(time,resp_n(:,i,lower),':','color','black','LineWidth',2); hold on; 
    title('Labor'); grid on;
    hold off;
    
    figure(i+4)
    subplot(2,3,1); [y, x] = hist(squeeze(var_decomp(1,i,:)),20); plot(x,y); axis([0 1 0 500 ]); axis 'auto y'; title('Inflation');
    subplot(2,3,2); [y, x] = hist(squeeze(var_decomp(2,i,:)),20); plot(x,y); axis([0 1 0 500 ]); axis 'auto y'; title('Output');
    subplot(2,3,3); [y, x] = hist(squeeze(var_decomp(3,i,:)),20); plot(x,y); axis([0 1 0 500 ]); axis 'auto y'; title('Output Gap');
    subplot(2,3,4); [y, x] = hist(squeeze(var_decomp(4,i,:)),20); plot(x,y); axis([0 1 0 500 ]); axis 'auto y'; title('Nominal Interest Rate');
    subplot(2,3,5); [y, x] = hist(squeeze(var_decomp(5,i,:)),20); plot(x,y); axis([0 1 0 500 ]); axis 'auto y'; title('Labor'); 
end


%% Part (4)
% Print Impulse Responses - Monetary Shock
figure(1)
print('impulse_monetary','-dpng')
% Print Impulse Responses - Productivity Shock
figure(2)
print('impulse_prod','-dpng')
% Print Impulse Responses - Demand Shock
figure(3)
print('impulse_demand','-dpng')
% Print Impulse Responses - Cost-push Shock
figure(4)
print('impulse_cost','-dpng')


%% Part (5)
% Print Variance Decomposition - Monetary Shock
figure(5)
print('var_monetary','-dpng')
% Print Variance Decomposition - Productivity Shock
figure(6)
print('var_prod','-dpng')
% Print Variance Decomposition - Demand Shock
figure(7)
print('var_demand','-dpng')
% Print Variance Decomposition - Cost-push Shock
figure(8)
print('var_cost','-dpng')


%% Part (6)
disp(['Probability(Labor<0 in period 8): ', num2str(mean(resp_n(8,2,:)<0))]);
