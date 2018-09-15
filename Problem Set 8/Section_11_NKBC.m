% Tilahun Emiru
% November 10, 2017


close all;
clear all;

% Calibration
tau       = 2;  % Since I am using sigma for the standard deviaion of the shocks, I am using tau to denote the CRRA parameter.  
beta      = 0.95; % discount factor
theta     = 2/3; % degree of price stickiness
phi_pi    = 1.5; % taylor rule parameter
phi_y     = 0.5/4; % taylor rule parameter
varphi    = 0.25; %inverse of elastiicity of labor supply
alpha     = 1/3; %production function parameter
eps       = 6; % elasticity of substitution between goods i and j in the consumption basket
rho_v     = 0.5; %persistence parameter
rho_a     = 0.9; %persistence parameter
sigma_v   = 0.25; %standard deviation
sigma_a   = 1.16; %standard deviation

% -------------------------------------------------------------------------
% Solving a linear rational expectations model through the method of
% undetermined coefficients
% -------------------------------------------------------------------------

v(1) = sigma_v;
a(1) = sigma_a;
for i = 2:12
    v(i) = rho_v*v(i-1);
end

% coefficients obtained from the solution
psi_yna = (1+varphi)/(tau*(1-alpha)+varphi+alpha);
kappa   = (1/theta-1)*(1-beta*theta)*((1-alpha)/(1-alpha+alpha*eps))*(tau+((varphi+alpha)/(1-alpha)));


phi_na = (1+varphi)/(tau*(1-alpha)+varphi+alpha);
mu     = log (eps/(eps-1));
theta_yn = -((1-alpha)*(mu - log(1-alpha)))/(tau*(1-alpha)+varphi+alpha);
lambda_v = ((1-beta*rho_v)*((tau*(1-rho_v)+phi_y))+kappa*(phi_pi-rho_v))^(-1);

% use the coeffients to map the monetary policy shocks to the outcome
% variabels 
y_gap = -(1-beta*rho_v)*lambda_v*v;
pi_t = -kappa*lambda_v*v;
r_t = tau*(1-rho_v)*(1-beta*rho_v)*lambda_v*v;
i_t = (tau*(1-rho_v)*(1-beta*rho_v)-rho_v*kappa)*lambda_v*v;

figure(1)
subplot(3,2,1); plot(y_gap, '-O','Linewidth',1.4); xlabel('Output Gap'); grid on;
subplot(3,2,2); plot(pi_t, '-O','Linewidth',1.4); xlabel('Inflation'); grid on;
subplot(3,2,3); plot(i_t, '-O','Linewidth',1.4); xlabel('Nominal Rate'); grid on;
subplot(3,2,4); plot(r_t, '-O','Linewidth',1.4); xlabel('Real Rate'); grid on;
subplot(3,2,5); plot(v, '-O','Linewidth',1.4); xlabel('V'); grid on;

%%
% -------------------------------------------------------------------------
% Solving a linear rational expectations model using gensys
% -------------------------------------------------------------------------

% Creating coefficient matrices for gensys
Gamma0 = [1 0 (tau^(-1)) (-tau^(-1)) 0 0 -1 (-tau^(-1));
          -kappa 1 zeros(1,5) -beta;
          -phi_y -phi_pi 1 0 -1 zeros(1,3);
          zeros(1,3) 1 0 tau*psi_yna*(1-rho_a) 0 0;
          zeros(1,4) 1 zeros(1,3);
          zeros(1,5) 1 0 0;
          1 zeros(1,7);
          0 1 zeros(1,6)];
      
Gamma1 = [zeros(4,8);
          zeros(1,4) rho_v zeros(1,3);
          zeros(1,5) rho_a 0 0;
          zeros(1,6) 1 0;
          zeros(1,7) 1];
Psi    = [zeros(4,2);
          eye(2);
          zeros(2,2)];
Pi     = [zeros(6,2);
          eye(2)];      
C      = zeros(8,1);

% Solving the model using gensys - before you use the gensys function, make
% sure you have the functions qzdiv and qzswitch functions in the current folder you are
% working in. Both of these functions are available from Chris Sims webpage

[T,c,R,fmat,fwt,ywt,gev,eu,loose]=gensys(Gamma0,Gamma1,C,Psi,Pi);


% Plotting the Impulse Responses 
irhor = 12; % impulse response horizon
time = [1:irhor];
col = R;    % R is a (2x8) matrix, the first column represents the effect of 
            % monetary policy shocks on the eight variables and the second 
            % will be the effect of technology shocks.  
for j=1:irhor
    resp(:,:,j)=col; % this defines an array -- the first matrix is the 
                     % initial impact matrix; the second matrix will be the 
                     % effect after one 
    col=T*col; % this equation updates the effect of the shocks for each 
               % subsequent period. 
end

% the matlab function squeeze converts arrays into vectors or matrices. 
resp1y(:,1)=squeeze(resp(1,1,:)); % this picks element (1,1) of each of the matrices in the 12 dimensional array and convert them into a vector
resp2y(:,1)=squeeze(resp(1,2,:)); % this picks element (1,2) of each of the matrices in the 12 dimensional array and convert them into a vector
resp1pi(:,1)=squeeze(resp(2,1,:));
resp2pi(:,1)=squeeze(resp(2,2,:));
resp1i(:,1)=squeeze(resp(3,1,:));
resp2i(:,1)=squeeze(resp(3,2,:));
resp1r(:,1)=squeeze(resp(4,1,:));
resp2r(:,1)=squeeze(resp(4,2,:));
resp1v(:,1)=squeeze(resp(5,1,:));
resp2v(:,1)=squeeze(resp(5,2,:));
resp1a(:,1)=squeeze(resp(6,1,:));
resp2a(:,1)=squeeze(resp(6,2,:));
resp1ye(:,1)=squeeze(resp(7,1,:));
resp2ye(:,1)=squeeze(resp(7,2,:));
resp1pie(:,1)=squeeze(resp(8,1,:));
resp2pie(:,1)=squeeze(resp(8,2,:));


% Monetary Policy Shocks
figure(2)
subplot(4,2,1)
plot(time,resp1y(:,1),'-O', 'linewidth', 1.4)
title('Output Gap')
subplot(4,2,3)
plot(time,resp1i(:,1),'-O', 'linewidth', 1.4)
title('Nominal Interest Rates')
subplot(4,2,7)
plot(time,resp1a(:,1),'-O', 'linewidth', 1.4); ylim([-0.2 0.2]);
title('a_t')
subplot(4,2,5)
plot(time,resp1ye(:,1),'-O', 'linewidth', 1.4)
title('E(y)')
subplot(4,2,2)
plot(time,resp1pi(:,1),'-O', 'linewidth', 1.4)
title('Inflation')
subplot(4,2,4)
plot(time,resp1r(:,1),'-O', 'linewidth', 1.4); ylim([-0.2 0.2])
title('Real Interest Rate')
subplot(4,2,8)
plot(time,resp1v(:,1),'-O', 'linewidth', 1.4)
title('\nu_t')
subplot(4,2,6)
plot(time,resp1ye(:,1),'-O', 'linewidth', 1.4)
title('E(\pi)')
saveas(gcf, 'section11.png')
