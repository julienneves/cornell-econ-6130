% Set up and estimate SVAR(1)
clc
clear all
close all
load Z_raw;% Order Interest Rate, Inflation, Output
%--------------------------------------------------------------------------
% Set control parameters
%--------------------------------------------------------------------------

%Number of draws
S=3e6;%Number of draws in MCMC
burnin=0.2*S;% Fraction of MCMC disregarded as "burnin sample"
J=500;%Number of draws from MCMC used to simulate posterior functions of theta
epseye=1e-6; %scale up or down to tune acceptance ratio
adaptive=1; %set to 1 to use adaptive proposal density
dates=1958:0.25:2017.5;

%post Volcker - pre-crisis sample
% Z_raw=Z_raw(:,100:end-40);
% dates=dates(100:end-40);

% %Pre-crisis sample
Z_raw=Z_raw(:,1:end-40);
dates=dates(1:end-40);

% Detrend data
[r_cycle, r_trend] = hp_detrend(0.01*(Z_raw(1,:)'), 1600);
[pi_cycle, pi_trend] = hp_detrend(0.01*(Z_raw(2,:)'), 1600);
[y_cycle, y_trend] = hp_detrend(log(Z_raw(3,:)'), 1600);

Z=[r_cycle';pi_cycle';y_cycle';];
%plot the data
figure(1)
subplot(3,1,1);
plot(dates,Z_raw(1,:),'linewidth',2);title('Federal Funds Rate','fontsize',16);axis([dates(1) dates(end) 0 500 ]); axis 'auto y';
subplot(3,1,2);
plot(dates,Z_raw(2,:),'linewidth',2);title('Inflation','fontsize',16);axis([dates(1) dates(end) 0 500 ]); axis 'auto y';
subplot(3,1,3);
plot(dates,Z_raw(3,:),'linewidth',2);title('Real GDP','fontsize',16);axis([dates(1) dates(end) 0 500 ]); axis 'auto y';

figure(2)
subplot(3,1,1);
plot(dates,Z(1,:),'linewidth',2);title('Federal Funds Rate','fontsize',16);axis([dates(1) dates(end) 0 500 ]); axis 'auto y';
subplot(3,1,2);
plot(dates,Z(2,:),'linewidth',2);title('Inflation','fontsize',16);axis([dates(1) dates(end) 0 500 ]); axis 'auto y';
subplot(3,1,3);
plot(dates,Z(3,:),'linewidth',2);title('Real GDP','fontsize',16);axis([dates(1) dates(end) 0 500 ]); axis 'auto y';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting values from OLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=length(Z);
Y=Z(:,2:end);
X=Z(:,1:end-1);
beta=Y*X'/(X*X');
Omega_ols= (1/(T-1))*(Y-beta*X)*(Y-beta*X)';
C_ols=chol(Omega_ols)';
Phi_ols=beta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate posterior using Metroplis Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta=[C_ols(1,1),C_ols(2,1),C_ols(2,2),C_ols(3,1),C_ols(3,2),C_ols(3,3),reshape(Phi_ols,1,9);]';%Starting value for parameter vector
LB=-100*ones(15,1);LB(1)=0;LB(3)=0;LB(6)=0;
UB=100*ones(15,1);
x=theta;


%--------------------------------------------------------------------------
% Initializes the MH algorithm
%--------------------------------------------------------------------------

%Initializes the proposal variance.
%Set up so that the first candidate draw is always accepted
lpostdraw = -9e+200;
bdraw=x;%Initial draw
vscale=diag(abs(theta))*1d-3+1e-5*eye(length(x));%Initial covariacne of proposal density (not really important)

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
    iter;
    % Draw from proposal density Theta*_{t+1} ~ N(Theta_{t},vscale)
    bcan = bdraw + norm_rnd(vscale);
    
    if min(bcan > LB)==1 %check bounds
        if min(bcan < UB)==1
            lpostcan = LLSVAR_P1(bcan,Z);
            laccprob = lpostcan-lpostdraw;
        else
            laccprob=-9e+200; %assign very low value to reject
            q=q+1;
        end
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %use for adaptive M-H
    if mod(iter,10000)==0 && iter > 10000 %do this every 10000 draws if iter > 10000
        disp(['iter: ',num2str(iter)]);
        disp(['acceptance rate: ',num2str(SwitchesProp(iter))]);
        if adaptive==1
            vscale=3d-1*cov(thetaMCMC(:,1000:iter)')+1e-10*eye(15);%update the covariance of proposal density (the second component is just to avoid possible singualrites o
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
toc

disp(['iter: ',num2str(iter)]);
disp(['acceptance rate: ',num2str(SwitchesProp(iter))]);

thetaMCMC=thetaMCMC(:,burnin:100:end); % disregard burnin sample and keep only evey 100th draw
convcheck(thetaMCMC);
plotpost(thetaMCMC,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find and plot probability intervals for IRFS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
periods=25;
ra = max(size(thetaMCMC));
FF = ceil(ra.*rand(J,1));

IMPR1=zeros(3,periods,J);
IMPR2=zeros(3,periods,J);
IMPR3=zeros(3,periods,J);

VARDECOMP1=zeros(3,J);
VARDECOMP2=zeros(3,J);
VARDECOMP3=zeros(3,J);
count=0;
D=[eye(3),zeros(3,9);];
for j=1:J
    theta=thetaMCMC(:,FF(j,1));
    
    C=[theta(1),0,0;
        theta(2),theta(3),0;
        theta(4),theta(5),theta(6);];
    
    Phi=reshape(theta(7:15),3,3);
    
    IR1=zeros(3,periods);
    IR2=zeros(3,periods);
    IR3=zeros(3,periods);
    for jj=0:periods-1
        IR1(:,jj+1)= (Phi^jj)*C(:,1);%MP shock
        IR2(:,jj+1)= (Phi^jj)*C(:,2);%Inflation shock
        IR3(:,jj+1)= (Phi^jj)*C(:,3);%Output Shock
    end
    
    count = count + (IR1(3,3)< 0); %add count if FFR on GDP (3rd variable) after 2 quarters (i.e. in period 3) is < 0
    
    IMPR1(:,:,j)=IR1;
    IMPR2(:,:,j)=IR2;
    IMPR3(:,:,j)=IR3;
    
    %     Variance decomp
    sigx=dlyap(Phi,C*C');
    sigx1=dlyap(Phi,C(:,1)*C(:,1)');
    sigx2=dlyap(Phi,C(:,2)*C(:,2)');
    sigx3=dlyap(Phi,C(:,3)*C(:,3)');
    VARDECOMP1(:,j)=diag(sigx1)./diag(sigx);
    VARDECOMP2(:,j)=diag(sigx2)./diag(sigx);
    VARDECOMP3(:,j)=diag(sigx3)./diag(sigx);
    %
end


display('Prob(FFR on GDP t+2 < 0)')
count/J

ImpSort1=sort(IMPR1,3);
ImpSort2=sort(IMPR2,3);
ImpSort3=sort(IMPR3,3);

upper=ceil(J*.95);
med=ceil(J*.5);
lower=ceil(J*.05);

figure(5)

subplot(3,3,1);
xlabel('FFR on FFR')
hold on;
plot(ImpSort1(1,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort1(1,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort1(1,:,lower),':','color','black','LineWidth',3);
hold on;

subplot(3,3,4);
xlabel('FFR on Inf')
hold on;
plot(ImpSort1(2,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort1(2,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort1(2,:,lower),':','color','black','LineWidth',3);
hold on;

subplot(3,3,7);
xlabel('FFR on GDP')
hold on;
plot(ImpSort1(3,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort1(3,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort1(3,:,lower),':','color','black','LineWidth',3);
hold on;

%Plot IRF to inflation shock
subplot(3,3,2);
xlabel('Inf on FFR')
hold on;
plot(ImpSort2(1,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort2(1,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort2(1,:,lower),':','color','black','LineWidth',3);
hold on;

subplot(3,3,5);
xlabel('Inf on Inf')
hold on;
plot(ImpSort2(2,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort2(2,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort2(2,:,lower),':','color','black','LineWidth',3);
hold on;

subplot(3,3,8);
xlabel('Inf on GDP')
hold on;
plot(ImpSort2(3,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort2(3,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort2(3,:,lower),':','color','black','LineWidth',3);
hold on;


%Plot IRF to inflation shock
subplot(3,3,2);
xlabel('Inf on FFR')
hold on;
plot(ImpSort2(1,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort2(1,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort2(1,:,lower),':','color','black','LineWidth',3);
hold on;

subplot(3,3,5);
xlabel('Inf on Inf')
hold on;
plot(ImpSort2(2,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort2(2,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort2(2,:,lower),':','color','black','LineWidth',3);
hold on;

subplot(3,3,8);
xlabel('Inf on GDP')
hold on;
plot(ImpSort2(3,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort2(3,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort2(3,:,lower),':','color','black','LineWidth',3);
hold on;



%Plot IRF to GDP shock
subplot(3,3,3);
xlabel('GDP on FFR')
hold on;
plot(ImpSort3(1,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort3(1,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort3(1,:,lower),':','color','black','LineWidth',3);
hold on;

subplot(3,3,6);
xlabel('GDP on Inf')
hold on;
plot(ImpSort3(2,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort3(2,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort3(2,:,lower),':','color','black','LineWidth',3);
hold on;

subplot(3,3,9);
xlabel('GDP on GDP')
hold on;
plot(ImpSort3(3,:,upper),':','color','black','LineWidth',3);
hold on;
plot(ImpSort3(3,:,med),'color','black','LineWidth',3);
hold on;
plot(ImpSort3(3,:,lower),':','color','black','LineWidth',3);
hold on;


%plot unsorted IRFs
figure(6)
subplot(3,1,1);
plot(reshape(IMPR1(1,:,:),periods,J));
xlabel('FFR on FFR')
subplot(3,1,2);
plot(reshape(IMPR1(2,:,:),periods,J));
xlabel('FFR on Inf')
subplot(3,1,3);
plot(reshape(IMPR1(3,:,:),periods,J));
xlabel('FFR on GDP')

%
figure(7)
subplot(3,3,1)
hist(VARDECOMP1(1,:),50);axis([0 1 0 500 ]); axis 'auto y';
xlabel('FFR shocks of FFR var')
subplot(3,3,4);
hist(VARDECOMP1(2,:),50);axis([0 1 0 500 ]); axis 'auto y';
xlabel('FFR shocks of Inf var')
subplot(3,3,7)
hist(VARDECOMP1(3,:),50);axis([0 1 0 500 ]); axis 'auto y';
xlabel('FFR shocks of GDP var')


subplot(3,3,2)
hist(VARDECOMP2(1,:),50);axis([0 1 0 500 ]); axis 'auto y';
xlabel('Inf shocks of FFR var')
subplot(3,3,5);
hist(VARDECOMP2(2,:),50);axis([0 1 0 500 ]); axis 'auto y';
xlabel('Inf shocks of Inf var')
subplot(3,3,8)
hist(VARDECOMP2(3,:),50);axis([0 1 0 500 ]); axis 'auto y';
xlabel('Inf shocks of GDP var')


subplot(3,3,3)
hist(VARDECOMP3(1,:),50);axis([0 1 0 500 ]); axis 'auto y';
xlabel('GDP shocks of FFR var')
subplot(3,3,6);
hist(VARDECOMP3(2,:),50);axis([0 1 0 500 ]); axis 'auto y';
xlabel('GDP shocks of Inf var')
subplot(3,3,9)
hist(VARDECOMP3(3,:),50);axis([0 1 0 500 ]); axis 'auto y';
xlabel('GDP shocks of GDP var')


