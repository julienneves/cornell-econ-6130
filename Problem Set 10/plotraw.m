function x=plotraw(MCMC,raw)
[n,m] = size(MCMC);
sqrn=n^.5;

figure('Name','MCMC')
for j=1:n;
%     subplot(ceil(sqrn),ceil(sqrn),j);
subplot(ceil(sqrn),ceil(sqrn),j);
   plot(MCMC(j,:));
    if j==1;
        xlabel({'\sigma '});
    end;
    if j==2;
        xlabel('\beta ');
    end;
    if j==3;
        xlabel('\phi ');
    end;
    if j==4;
        xlabel('\epsilon ');
    end;
    if j==5;
        xlabel('\phi_{\pi} ');
    end;
    if j==6;
        xlabel('\phi_{y} ');
    end;
    if j==7;
        xlabel('\theta ');
    end;
    if j==8;
        xlabel('\alpha ');
    end;
    if j==9;
        xlabel('\rho_{v} ');
    end;
    if j==10;
        xlabel('\rho_{a} ');
    end;
    if j==11;
        xlabel('\rho_{z}');
    end;
    if j==12;
        xlabel('\rho_{u}');
    end;
    if j==13;
        xlabel('\sigma_{v} ');
    end;
    if j==14;
        xlabel('\sigma_{a} ');
    end;
    if j==15;
        xlabel('\sigma_{z}');
    end;
    if j==16;
        xlabel('\sigma_{u}');
    end;
end
if raw == 1
    print('raw_MCMC','-dpng')
else
    print('burn_MCMC','-dpng')
end
end