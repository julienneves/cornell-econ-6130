function [A, C, D, eu, R] = NKBC_model( THETA, type)
%NKBC_model New Keynesian Model with cost-push shocks
%   State space model:
%   X(t) = A*X(t-1) + C*eps(t)
%   Z(t) = D*X(t)
%   By Julien Neves

%% Matrix A and C
% Calibration
sigma     = THETA(1); % CRRA parameter.
beta      = THETA(2); % discount factor
phi       = THETA(3); % inverse of elastiicity of labor supply
eps       = THETA(4); % elasticity of substitution between goods i and j
phi_pi    = THETA(5); % taylor rule parameter
phi_y     = THETA(6); % taylor rule parameter
theta     = THETA(7); % degree of price stickiness
alpha     = THETA(8); % production function parameter
rho_v     = THETA(9); % persistence parameter
rho_a     = THETA(10); % persistence parameter
rho_z     = THETA(11); % persistence parameter
rho_u     = THETA(12); % persistence parameter
sigma_v   = THETA(13); % standard deviation
sigma_a   = THETA(14); % standard deviation
sigma_z   = THETA(15); % standard deviation
sigma_u   = THETA(16); % standard deviation

% Compute the coefficients
rho = -log(beta);
lambda = (1-theta)*(1-beta*theta)*(1-alpha)/(theta*(1-alpha+alpha*eps));
kappa = lambda*(sigma + (phi+alpha)/(1-alpha));
psi_ya = (1+phi)/(sigma*(1-alpha)+phi+alpha);

% State: 'y','x'; 'pi'; 'r^e'; 'i'; 'v'; 'a'; 'z'; 'u'; 'E(x)'; 'E(pi)'
Gamma0 = [kappa -kappa 0 0 0 0 0 0 -1 0 0;
    0 -kappa 1 0 0 0 0 0 -1 0 -beta;
    0 1 0 -1/sigma 1/sigma 0 0 0 0 -1 -1/sigma;
    -phi_y 0 -phi_pi 0 1 -1 -phi_y*psi_ya 0 0 0 0;
    0 0 0 1 0 0 sigma*(1-rho_a)*psi_ya -(1-rho_z) 0 0 0;
    0 0 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0;
    0 1 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0];

Gamma1 = [0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 rho_v 0 0 0 0 0;
          0 0 0 0 0 0 rho_a 0 0 0 0;
          0 0 0 0 0 0 0 rho_z 0 0 0;
          0 0 0 0 0 0 0 0 rho_u 0 0;
          0 0 0 0 0 0 0 0 0 1 0;
          0 0 0 0 0 0 0 0 0 0 1];

Psi = [ 0 0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0 1 0 0;]';

Pi = [0 0 0 0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 0 0 0 1]';

Cons = [0 0 0 0 0 0 0 0 0 0 0]';

% Solve New Keynesian Model
[A,~,R,~,~,~,~,eu,~]=gensys(Gamma0,Gamma1,Cons,Psi,Pi);

C = R*[sigma_v 0 0 0;
       0 sigma_a 0 0;
       0 0 sigma_z 0;
       0 0 0 sigma_u];

%% Matrix D
% Set up measurement matrix Z(t)=D*S(t)
if strcmp(type, 'data')
    % State: 'y','x'; 'pi'; 'r^e'; 'i'; 'v'; 'a'; 'z'; 'u'; 'E(x)'; 'E(pi)'
    D = [0 0 1 0 0 0 0 0 0 0 0; % inflation
        1 0 0 0 0 0 psi_ya 0 0 0 0; % output
        0 0 0 0 1 0 0 0 0 0 0;  % nominal interest rate
        1/(1-alpha) 0 0 0 0 0 -(1-psi_ya)/(1-alpha) 0 0 0 0]; % labor
    
elseif strcmp(type, 'impulse')
    % State: 'y','x'; 'pi'; 'r^e'; 'i'; 'v'; 'a'; 'z'; 'u'; 'E(x)'; 'E(pi)'
    D = [0 0 1 0 0 0 0 0 0 0 0; % inflation
        1 0 0 0 0 0 psi_ya 0 0 0 0; % output
        1 0 0 0 0 0 0 0 0 0 0;  % output gap
        0 0 0 0 1 0 0 0 0 0 0;  % nominal interest rate
        1/(1-alpha) 0 0 0 0 0 -(1-psi_ya)/(1-alpha) 0 0 0 0]; % labor
else
    warning('Measurement matrix missing')
end

end

