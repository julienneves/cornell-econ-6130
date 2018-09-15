% Set up and estimate miniture DSGE model
clc
clear all
close all
global Z
load('Z'); %load data. Order of variables: Inflation, output, interest rate, labor supply

figure(1)
subplot(2,2,1);
plot(Z(1,:),'linewidth',2);title('Output','fontsize',16);
subplot(2,2,2);
plot(Z(2,:),'linewidth',2);title('Inflation','fontsize',16)
subplot(2,2,3);
plot(Z(3,:),'linewidth',2);title('Nominal Interest Rate','fontsize',16)
subplot(2,2,4);
plot(Z(4,:),'linewidth',2);title('Labor ','fontsize',16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial values of structural parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau       = 3;  % Since I am using sigma for the standard deviaion of the shocks, I am using tau to denote the CRRA parameter.
beta      = 0.99; % discount factor
theta     = 3/4; % degree of price stickiness
phi_pi    = 1.5; % taylor rule parameter
phi_y     = 0.125; % taylor rule parameter
varphi    = 3; %inverse of elastiicity of labor supply
alpha     = 1/3; %production function parameter
eps       = 5; % elasticity of substitution between goods i and j in the consumption basket
rho_v     = 0.5; %persistence parameter
rho_a     = 0.75; %persistence parameter
sigma_v   = 0.02; %standard deviation
sigma_a   = (0.012^2*(1-rho_a^2))^.5; %standard deviation of innovation to a_t


THETA=[tau,beta,theta,phi_pi,phi_y,varphi,alpha,eps,rho_v,rho_a,sigma_v,sigma_a;]';%Starting value for parameter vector
LB=[0,0,0,1,0,1,0,1,zeros(1,4)]';%Lower bound for parameter vector
UB=[10,1,1,5,5,10,1,25,1,1,10,10]';%Upper bound for parameter vector
x=THETA;

% x=THETA;
sa_t= 5;
sa_rt=.3;
sa_nt=5;
sa_ns=5;
% warning off all;

[xhat]=simannb( 'LLDSGE', x, LB, UB, sa_t, sa_rt, sa_nt, sa_ns, 1);

%--------------------------------------------------------------------------
thetalabel=['tau    ';'beta   ';'theta  ';'phi_pi ';'phi_y  ';'varphi ';'alpha  ';'eps    ';'rho_v  ';'rho_a  ';'sigmav ';'sigmaa ';];
disp('ML estimate of THETA')
disp([thetalabel, num2str(xhat)])
%--------------------------------------------------------------------------
a_hat = filter_gap_DSGE(xhat,Z);
figure(2)
plot(Z(1,:),'linewidth',2);hold on;plot(a_hat,'linewidth',2)
legend('Output','Output gap')