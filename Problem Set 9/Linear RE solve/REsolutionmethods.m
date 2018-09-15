% Solving simple New Keyensian model
%Kristoffer Nimark, Cornell % university

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear var
clc
keep_plot=1;
if keep_plot ==0
    close all
end
shock=2;%1=prod, 2=demand, 3=monetary policy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign values to parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_a=.9; %   Persistence of potential output
r_z=.9; %   Persistence of demand shocks
r_v=.9; %   Persistence of monetary policy shocks

b=.95;  %  Discount rate (beta)
sig = 1; % Curvature in consumption
varphi= 5 ; %curvature in labor supply
epsilon= 9 ; % CES aggregator elasticity

fi_pi=1.5 ; % Taylor rule parameter on inflation
fi_y=0.125 ; % Taylor rule parameter on output
theta=0.75 ; %Calvo parameter

alpha = .5; % Labor share in production function
lambda = ((1-theta)*(1-theta*b)/theta)* ( (1-alpha)/ (1-alpha+epsilon*alpha));
k = lambda * (sig + (varphi + (alpha + varphi)/(1-alpha)));
psi_ya = (1 + varphi)/(sig*(1-alpha) + varphi + alpha);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soderlind of stable/unstable decoupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cutoff=.999999; %Define the cutoff for stable vs unstable eignevalues (should be just below unity).

% Define model matrices
A0=[1,0,0,0,0;
    0,1,0,0,0;
    0,0,1,0,0;
    0,0,0,b,0;
    0,0,0,-sig,1;];

A1=[r_a,0,0,0,0;
    0,r_z,0,0,0;
    0,0,r_v,0,0;
    0,0,0,1,-k;
    sig*(fi_y*psi_ya +sig*(1-r_a)*psi_ya),-sig*(1-r_z),sig,sig*fi_pi,1+sig*fi_y;];

C1= [eye(3);zeros(2,3);];

A=A0\A1;
C=A0\C1;
C=C(1:3,1:3);

egen = abs(eig(A)) < cutoff;

n1=3; %Number of predetermined variables
n2=2; % Number of jump variables
n = n1 + n2; %Total number of variables

%MatLab, complex generalized Schur decomposition
[S,T,Qa,Z] = qz(eye(size(A)),A); %MatLab: I=Q'SZ' and A=Q'TZ'; Paul S:  I=QSZ' and A=QTZ',%but Q isn't used
[S,T,Qa,Z] = reorder(S,T,Qa,Z);   % reordering of generalized eigenvalues, T(i,i)/S(i,i), in ascending order
logcon = abs(diag(T)) <= (abs(diag(S))*cutoff);  %1 for stable eigenvalue

if sum(logcon) < n1
    warning('Too few stable roots: no stable solution');
    M = NaN; G = NaN;J0 = NaN;
    return;
elseif sum(logcon) > n1
    warning('Too many stable roots: inifite number of stable solutions');
    M = NaN; G = NaN;J0 = NaN;
    return;
end


Stt = S(1:n1,1:n1);
Zkt = Z(1:n1,1:n1);
Zlt = Z(n1+1:n,1:n1);
Ttt = T(1:n1,1:n1);

if cond(Zkt) > 1e+14
    warning('Zkt is singular: rank condition for solution not satisfied');
    M = NaN; G = NaN;J0 = NaN;
    return;
end

Zkt_1 = inv(Zkt);         %inverting
Stt_1 = inv(Stt);


M = real(Zkt*Stt_1*Ttt*Zkt_1);       %x1(t+1) = M*x1(t) + e(t+1)
G = real(Zlt*Zkt_1) ;              %x2(t) = G*x1(t)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display solution/output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc %display time passed since "tic"

display('Stable/unstable eigenvalue decoupling');

G

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct equilibrium functions for output, labor, interest rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_pi=G(1,:);
G_yg = G(2,:);
G_y = G(2,:) + [psi_ya ,0,0;];
G_n = (1/(1-alpha))*(G_y - [1,0,0;]);
G_i= fi_y*G_y + fi_pi*G_pi + [0,0,1;];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Impulse responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

periods = 50;
IRF = zeros(6,periods);
for t=1:periods
    IRF(1:5,t)=[G_pi;G_yg;G_y;G_n;G_i;]*M^(t-1)*C(:,shock);
    if t==1
        IRF(6,t)=G_pi*M^(t-1)*C(:,shock);
    else
        IRF(6,t)=IRF(6,t-1)+G_pi*M^(t-1)*C(:,shock);
    end
end

if keep_plot==0
    figure(1)
    subplot(2,3,1);
    plot(IRF(1,:),'linewidth',2);title('Inflation','fontsize',16);
    subplot(2,3,2);
    plot(IRF(2,:),'linewidth',2);title('Output Gap','fontsize',16)
    subplot(2,3,3);
    plot(IRF(3,:),'linewidth',2);title('Output','fontsize',16)
    subplot(2,3,4);
    plot(IRF(4,:),'linewidth',2);title('Labor inputs','fontsize',16)
    subplot(2,3,5);
    plot(IRF(5,:),'linewidth',2);title('Nominal Interest Rate','fontsize',16)
    subplot(2,3,6);
    plot(IRF(6,:),'linewidth',2);title('Price Level','fontsize',16)
    
else
    figure(1)
    hold on;subplot(2,3,1);
    hold on;plot(IRF(1,:),'linewidth',2,'linestyle',':');title('Inflation','fontsize',16);
    hold on;subplot(2,3,2);
    hold on;plot(IRF(2,:),'linewidth',2,'linestyle',':');title('Output Gap','fontsize',16)
    hold on;subplot(2,3,3);
    hold on;plot(IRF(3,:),'linewidth',2,'linestyle',':');title('Output','fontsize',16)
    hold on;subplot(2,3,4);
    hold on;plot(IRF(4,:),'linewidth',2,'linestyle',':');title('Labor inputs','fontsize',16)
    hold on;subplot(2,3,5);
    hold on;plot(IRF(5,:),'linewidth',2,'linestyle',':');title('Nominal Interest Rate','fontsize',16)
    hold on;subplot(2,3,6);
    hold on;plot(IRF(6,:),'linewidth',2,'linestyle',':');title('Price Level','fontsize',16)
end










