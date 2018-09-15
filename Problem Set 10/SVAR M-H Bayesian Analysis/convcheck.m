function x=convcheck(MCMC)
n=size(MCMC,1);
m=size(MCMC,2);
x=[];
mu=[];
for j=10:10:m
    X=diag(cov(MCMC(:,1:j)'));
    x=[x X];
    MU=mean(MCMC(:,1:j),2);
    mu=[mu MU;];
end


sqrn=n^.5;


figure
for j=1:n;
    %     subplot(ceil(sqrn),ceil(sqrn),j);
    subplot(5,3,j);
    plot(MCMC(j,:),'linewidth',2);
    if j==1;
        xlabel({'c_{11} '});
    end;
    if j==2;
        xlabel('c_{21} ');
        title('Raw MCMC','fontsize',20)

    end;
    if j==3;
        xlabel('c_{22} ');
    end;
    if j==4;
        xlabel('c_{31} ');
    end;
    if j==5;
        xlabel('c_{32} ');
    end;
    if j==6;
        xlabel('c_{33} ');
    end;
    if j==7;
        xlabel('\phi_{11}');
    end;
    
    if j==8;
        xlabel({'\phi_{12} '});
    end;
    if j==9;
        xlabel('\phi_{13} ');
    end;
    if j==10;
        xlabel('\phi_{21} ');
    end;
    if j==11;
        xlabel('\phi_{22} ');
    end;
    if j==12;
        xlabel('\phi_{23} ');
    end;
    if j==13;
        xlabel('\phi_{31} ');
    end;
    if j==14;
        xlabel('\phi_{32}');
    end;
    if j==15;
        xlabel('\phi_{33}');
    end;
    
    
    
end
figure
for j=1:n;
     %     subplot(ceil(sqrn),ceil(sqrn),j);
    subplot(5,3,j);
    plot(mu(j,:),'linewidth',2);
    if j==1;
        xlabel({'c_{11} '});
    end;
    if j==2;
        xlabel('c_{21} ');
        title('Recursive mean of MCMC','fontsize',20)

    end;
    if j==3;
        xlabel('c_{22} ');
    end;
    if j==4;
        xlabel('c_{31} ');
    end;
    if j==5;
        xlabel('c_{32} ');
    end;
    if j==6;
        xlabel('c_{33} ');
    end;
    if j==7;
        xlabel('\phi_{11}');
    end;
    
    if j==8;
        xlabel({'\phi_{12} '});
    end;
    if j==9;
        xlabel('\phi_{13} ');
    end;
    if j==10;
        xlabel('\phi_{21} ');
    end;
    if j==11;
        xlabel('\phi_{22} ');
    end;
    if j==12;
        xlabel('\phi_{23} ');
    end;
    if j==13;
        xlabel('\phi_{31} ');
    end;
    if j==14;
        xlabel('\phi_{32}');
    end;
    if j==15;
        xlabel('\phi_{33}');
    end;
    
end

figure

for j=1:n;
    
     %     subplot(ceil(sqrn),ceil(sqrn),j);
    subplot(5,3,j);
    plot(x(j,:),'linewidth',2);
    if j==1;
        xlabel({'c_{11} '});
    end;
    if j==2;
        xlabel('c_{21} ');
        title('Recursive variance of MCMC','fontsize',20)
    end;
    if j==3;
        xlabel('c_{22} ');
    end;
    if j==4;
        xlabel('c_{31} ');
    end;
    if j==5;
        xlabel('c_{32} ');
    end;
    if j==6;
        xlabel('c_{33} ');
    end;
    if j==7;
        xlabel('\phi_{11}');
    end;
    
    if j==8;
        xlabel({'\phi_{12} '});
    end;
    if j==9;
        xlabel('\phi_{13} ');
    end;
    if j==10;
        xlabel('\phi_{21} ');
    end;
    if j==11;
        xlabel('\phi_{22} ');
    end;
    if j==12;
        xlabel('\phi_{23} ');
    end;
    if j==13;
        xlabel('\phi_{31} ');
    end;
    if j==14;
        xlabel('\phi_{32}');
    end;
    if j==15;
        xlabel('\phi_{33}');
    end;
    
    
end


