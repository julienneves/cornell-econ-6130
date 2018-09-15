url = 'https://fred.stlouisfed.org/';
c = fred(url);

format bank

GDP = fetch(c,'GDP');
CON = fetch(c,'PCEC');
INV = fetch(c,'GPDI');

gdp_data = log(GDP.Data(:,2));
c_data = log(CON.Data(:,2));
i_data = log(INV.Data(:,2));

[~,gdp_data] = hpfilter(gdp_data,1600);
[~,i_data] = hpfilter(i_data,1600);
[~,c_data] = hpfilter(c_data,1600);

% Simulation
figure
plot([c_sim,c_data,gdp_sim,gdp_data,i_sim,i_data])

xlabel('t')
legend('C - Simulation','C - Data', 'GDP - Simulation','GDP - Data','I - Simulation','I - Data', 'Location', 'best')
print('plot_sim','-dpng')


close(c)