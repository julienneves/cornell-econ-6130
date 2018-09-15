% Kalman filter applications
% Code to play around with Kalman filter for state space system defined as
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	X[t] = AX[t-1] + Cu[t]
%
%   Z[t] = DX[t] + Ru[t]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all

periods=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bivariate state, 2 dim signal vector scalar signal
A=[0.95, 0; 0,0.5;];
C=eye(2);C(2,1)=-0.5;
D=eye(2);
R=eye(2);

%Bivariate state, signal sum of x1 and x2
A=[0.95, 0; 0,0.5;];
C=eye(2);C(2,1)=-0.5;
D=[1 1];
R=eye(1);

% %Bivariate state, signal is x1 + noise
A=[0.95, 0; 0,0.5;];
C=eye(2);C(2,1)=-0.5;
D=[1 0];
R=eye(1);

[K,P,Kss,Pss,Z]=kalman_filter_sim(A,C,D,R,periods);
