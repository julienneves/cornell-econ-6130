function x=convcheck(MCMC)
n=size(MCMC,1);
m=size(MCMC,2);
x=[];
for j=1:10000:m
    X=diag(cov(MCMC(:,100:j)'));    
    x=[x X];
end
%
sqrn=n^.5;
figure
for j=1:n;
    subplot(ceil(sqrn),ceil(sqrn),j);
    plot(x(j,:));
end

figure
for j=1:n;
%     subplot(ceil(sqrn),ceil(sqrn),j);
subplot(4,3,j);
   plot(MCMC(j,10:end));
     if j==1;
        xlabel({'\sigma '});
    end;
    if j==2;
        xlabel('\beta ');
    end;
    if j==3;
        xlabel('\theta ');
    end;
    if j==4;
        xlabel('\phi_{\pi} ');
    end;
    if j==5;
        xlabel('\phi_{y} ');
    end;
    if j==6;
        xlabel('\vartheta ');
    end;
   
    if j==7;
        xlabel('\alpha ');
    end;
    if j==8;
        xlabel('\epsilon ');
    end;
    if j==9;
        xlabel('\rho_{r} ');
    end;
    if j==10;
        xlabel('\rho_{a} ');
    end;
    if j==11;
        xlabel('\sigma_{v}');
    end;
    if j==12;
        xlabel('\sigma_{a}');
    end;
end