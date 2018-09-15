function [LL]=LLSVAR_P1(theta,Z)

Y=Z(:,2:end);
X=Z(:,1:end-1);

C=[theta(1),0,0;
    theta(2),theta(3),0;
    theta(4),theta(5),theta(6);];


Phi=reshape(theta(7:15),3,3);

CC=C*C';

RES=Y-Phi*X;
T=length(Y);
LL=-(T*2/2)*log(2*pi)-(T/2)*log(det(CC))-0.5*trace((RES'/CC)*RES);
if max(abs(eig(Phi))) >= 1  %disregard draws that imply non-stationary VAR
    ll=-1d20;
end