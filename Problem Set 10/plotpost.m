function plotpost(bb_,addplot)
%%
ysca=max(size(bb_));
Q=50; %HP smoothing parameter
% figure('yscaame','Calvo and Indexing','yscaumberTitle','off');
figure('Name','MCMC Distribution')
for j=1:length(bb_(:,1));
    hold on;
    [n,xout] = hist(bb_(j,:),100);
    n(1,1)=0;n(end,1)=0;
    xsca=(max(xout)-min(xout))/100;
    nn=hpfilter(n,Q);
    subplot(4,4,j);
    if addplot==1;
        AXX=axis;
        AX=[ min([min(xout),AXX(1,1)]) max([max(xout),AXX(1,2)]) 0 max([(max((nn*1.1))/(ysca*xsca)),AXX(1,4)])];
        plot(xout,nn/(ysca*xsca),'linewidth',2);axis(AX);
    else;
        plot(xout,nn/(ysca*xsca),'linewidth',2);axis([min(xout) max(xout) 0 max((nn*1.1))/(ysca*xsca)]);
    end
    
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
print('dist_MCMC','-dpng')
end
