function plotpost(bb_,addplot)
%%
ysca=max(size(bb_));
Q=50; %HP smoothing parameter
figure(8)
for j=1:length(bb_(:,1))
    hold on;
    [n,xout] = hist(bb_(j,:),100);
    n(1,1)=0;n(end,1)=0;
    xsca=(max(xout)-min(xout))/100;
    nn=hpfilter(n,Q);
       subplot(3,5,j);
    if addplot==1;
        AXX=axis;
        AX=[ min([min(xout),AXX(1,1)]) max([max(xout),AXX(1,2)]) 0 max([(max((nn*1.1))/(ysca*xsca)),AXX(1,4)])];
        plot(xout,nn/(ysca*xsca),'g','linewidth',2);axis(AX);
    else;
        plot(xout,nn/(ysca*xsca),'linewidth',2);axis([min(xout) max(xout) 0 max((nn*1.1))/(ysca*xsca)]);
    end
  if j==1;
        xlabel({'c_{11} '});
    end;
    if j==2;
        xlabel('c_{21} ');
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
