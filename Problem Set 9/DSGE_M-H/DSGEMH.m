% Set up and estimate miniture DSGE model
clc
clear all
close all
global Z
load('Z'); %load data. Order of variables: Output, inflation, interest rate, labor supply

%--------------------------------------------------------------------------
% Set control parameters
%--------------------------------------------------------------------------

%Number of draws
J=5e5;
epseye=2e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial values of structural parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau       = 4;  % Since I am using sigma for the standard deviaion of the shocks, I am using tau to denote the CRRA parameter.
beta      = 0.99; % discount factor
theta     = 2/4; % degree of price stickiness
phi_pi    = 3; % taylor rule parameter
phi_y     = 1; % taylor rule parameter
varphi    = 6; %inverse of elasticity of labor supply
alpha     = 1/3; %production function parameter
eps       = 10; % elasticity of substitution between goods i and j in the consumption basket
rho_v     = 0.7; %persistence parameter
rho_a     = 0.95; %persistence parameter
sigma_v   = 0.02; %standard deviation
sigma_a   = 0.007; %standard deviation of innovation to a_t


THETA=[tau,beta,theta,phi_pi,phi_y,varphi,alpha,eps,rho_v,rho_a,sigma_v,sigma_a;]';%Starting value for parameter vector
LB=[0,0,0,1,0,1,0,1,zeros(1,4)]';%Lower bound for parameter vector
UB=[10,1,1,5,5,10,1,25,1,1,10,10]';%Upper bound for parameter vector

x=THETA;

%--------------------------------------------------------------------------
% Initializes the MH algorithm
%--------------------------------------------------------------------------

%Initializes the proposal variance.
%Set up so that the first candidate draw is always accepted
lpostdraw = -9e+200;
bdraw=x;
vscale=diag(abs(THETA))*epseye+1e-12*eye(length(x));
% vscale=1e-20*eye(length(THETA));

bb_=zeros(length(x),J);%Store all draws in bb_
%Matrices that keep track of switches and drwas outside LB and UB
OutsideProp=zeros(J,1);
SwitchesProp=zeros(J,1);
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
for iter=1:J
    iter=iter+1;
    % Draw from proposal density Theta*_{t+1} ~ N(Theta_{t},vscale)
    bcan = bdraw + norm_rnd(vscale);
    
    if min(bcan > LB)==1
        if min(bcan < UB)==1
            %             lpostcan = LLDSGE(bcan);%Uncomment for improper uniform priors
            lpostcan = log_prior_DSGE(bcan)+LLDSGE(bcan);%switch on for use of priors
            %                         lpostcan = log_prior_DSGE(bcan);%Uncomment for prior predictive analysis
            laccprob = lpostcan-lpostdraw;
        else
            laccprob=-9e+200;
            q=q+1;
        end
    else
        laccprob=-9e+200;
        q=q+1;
    end
    
    %Accept candidate draw with log prob = laccprob, else keep old draw
    if log(rand)<laccprob
        lpostdraw=lpostcan;
        bdraw=bcan;
        pswitch=pswitch+1;
    end
    
    bb_(:,iter)=bdraw;
    
    OutsideProp(iter)=q/iter;
    SwitchesProp(iter)=pswitch/iter;
    if iter >= 1000 && mod(iter,1000)==0
        
        disp(['iter: ',num2str(iter)]);
        disp(['acceptance rate: ',num2str(SwitchesProp(iter))]);
        
        iter
    end
    
    if iter >= 10000 && mod(iter,1000)==0
        vscale=1e-2*cov(bb_(:,8000:iter)');
    end
    %
end
toc

figure(1)
subplot(2,2,1);
plot(Z(1,:),'linewidth',2);title('Output','fontsize',16);
subplot(2,2,2);
plot(Z(2,:),'linewidth',2);title('Inflation','fontsize',16)
subplot(2,2,3);
plot(Z(3,:),'linewidth',2);title('Nominal Interest Rate','fontsize',16)
subplot(2,2,4);
plot(Z(4,:),'linewidth',2);title('Labor ','fontsize',16)

figure
bb_=bb_(:,50:end);
for j=1:12;
    subplot(4,3,j);
    hist(bb_(j,:),50);
end

convcheck(bb_(:,200000:end));

figure
plotpost(bb_(:,200000:end),0)
%%