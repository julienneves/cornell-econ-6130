%% Setup
% Set number of periods
T = 100;
time = 0:T;

% Set state AR(1) matrices
A = [0.95 0; 0 0.5]; 
C = [1 0; -0.5 1];

%% Part a)-c)
disp('Part a-c')

% Set measurement matrix
D = [1 0; 0 1];

% Set measurement error matrix
S_vv = eye(2);

% Generate random errors
v = randn(2,T+1);
u = randn(2,T+1);

% Set states and observations matrices
X = zeros(2,T+1);
Z = zeros(2,T+1);

% Generate states and observations
for t=1:T
X(:,t+1) = A * X(:,t) + C * u(:,t);
Z(:,t+1) = D * X(:,t+1) + S_vv * v(:,t+1);
end

% Part b)
% Compute state variances
var_x = reshape((eye(4)-kron(A,A))^-1*reshape(C*C',4,1),2,2)
var_x_sample = cov(X')

% Compute observations variances
var_z = D*var_x*D' + S_vv
var_z_sample = cov(Z')

% Part c)
% Set starting values
X0 = zeros(2,1);
P0 = var_x;

% Compute Kalman Filter
[X_post, P_post] = kfilter(Z, A, C, D, S_vv, X0, P0 );

% Plot Kalman Filter estimates
figure
subplot(2,1,1)
plot(time, X(1,:),time, X_post(1,:), '--')
xlabel('t')
ylabel('X_1')
legend('X','X_{t|t}')
subplot(2,1,2)
plot(time, X(2,:),time, X_post(2,:), '--')
xlabel('t')
ylabel('X_2')
legend('X','X_{t|t}')

print('plot_c','-dpng')

% part d)

var_err = P_post
var_err_sample = cov((X-X_post)')

%% Part e)
disp('Part e')
% Part a)
% Set measurement matrix
D = [1 1];

% Generate observations
Z = D * X;

% Part b)
% Compute observations variances
var_z = D*var_x*D'
var_z_sample = cov(Z')

% Part c)
% Set starting values
X0 = zeros(2,1);
P0 = var_x;

% Compute Kalman Filter
X_post = kfilter(Z, A, C, D, 0, X0, P0 );

% Plot Kalman Filter estimates
figure
subplot(2,1,1)
plot(time, X(1,:),time, X_post(1,:), '--')
xlabel('t')
ylabel('X_1')
legend('X','X_{t|t}')
subplot(2,1,2)
plot(time, X(2,:),time, X_post(2,:), '--')
xlabel('t')
ylabel('X_2')
legend('X','X_{t|t}')

print('plot_e','-dpng')

%% Part f)
disp('Part f')
% Part a)
% Set measurement matrix
D = [1 0];

% Generate observations
Z = D * X

% Part b)
% Compute observations variances
var_z = D*var_x*D'
var_z_sample = cov(Z')

% Part c)
% Set starting values
X0 = zeros(2,1);
P0 = var_x;

% Compute Kalman Filter
[X_post, P_post] = kfilter(Z, A, C, D, 0, X0, P0 );

% Plot Kalman Filter estimates
figure
subplot(2,1,1)
plot(time, X(1,:),time, X_post(1,:), '--')
xlabel('t')
ylabel('X_1')
legend('X','X_{t|t}')
subplot(2,1,2)
plot(time, X(2,:),time, X_post(2,:), '--')
xlabel('t')
ylabel('X_2')
legend('X','X_{t|t}')

print('plot_f','-dpng')