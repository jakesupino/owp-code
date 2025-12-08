%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotFinalQC.m
% This script plots the final QC'd results for all three sites.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 5/15/2024
% Last updated: 1/19/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
cd([rootpath,'open-water-platform-data\north\cleaned\final-qc'])
load('finalQC.mat');
north = finalQC;

cd([rootpath,'open-water-platform-data\gull\cleaned\final-qc'])
load('finalQC.mat');
gull = finalQC;

cd([rootpath,'open-water-platform-data\south\cleaned\final-qc'])
load('finalQC.mat');
south = finalQC;

cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

% Find daily tidal range for each site
north_tidal = groupsummary(north,"datetime_utc","day","range","depth");
north_tidal.date = datetime(cellstr(north_tidal.day_datetime_utc),'TimeZone','UTC');
north_tidal = table2timetable(north_tidal);
north_tidal.GroupCount = [];
north_tidal.day_datetime_utc = [];
north_tidal.Properties.VariableNames = {'tidal'};

gull_tidal = groupsummary(gull,"datetime_utc","day","range","depth");
gull_tidal.date = datetime(cellstr(gull_tidal.day_datetime_utc),'TimeZone','UTC');
gull_tidal = table2timetable(gull_tidal);
gull_tidal.GroupCount = [];
gull_tidal.day_datetime_utc = [];
gull_tidal.Properties.VariableNames = {'tidal'};

south_tidal = groupsummary(south,"datetime_utc","day","range","depth");
south_tidal.date = datetime(cellstr(south_tidal.day_datetime_utc),'TimeZone','UTC');
south_tidal = table2timetable(south_tidal);
south_tidal.GroupCount = [];
south_tidal.day_datetime_utc = [];
south_tidal.Properties.VariableNames = {'tidal'};

% Find indices of deployment changes
ind_dep = find(diff(south.deployment) > 0);
dep = table([1;south.deployment(ind_dep+1)],[1;ind_dep+1]);
dep.Properties.VariableNames = {'depNum','ind'};

label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
    'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
    'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
    'Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};

%====Calculate the average deployment length===============================
dep_change = [south.datetime_utc(dep.ind);south.datetime_utc(end)];
for i = 2:length(dep_change)
    dt(i-1) = dep_change(i) - dep_change(i-1);
end

duration_avg = mean(dt);    % Mean deployment length in hh:mm:ss

duration_avg.Format = 'd';  % Mean deployment length in days --> Convert manually to weeks

%====Plot data=============================================================
% See Wong (2011) and https://www.rapidtables.com/convert/color/rgb-to-hex.html?r=86&g=180&b=233
north_clr = '#5D3A9B';
gull_clr = '#019E73';
south_clr = '#D55E00';

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(gull.datetime_utc,gull.depth,'.','Color',gull_clr,'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.depth,'.','Color',north_clr,'DisplayName','North')
plot(south.datetime_utc,south.depth,'.','Color',south_clr,'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Depth (m)')
legend('show','location','best')
title('Final QC''d data')

fig2 = figure(2);clf
fig2.WindowState = 'maximized';
plot(gull.datetime_utc,gull.salinity,'.','Color',gull_clr,'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.salinity,'.','Color',north_clr,'DisplayName','North')
plot(south.datetime_utc,south.salinity,'.','Color',south_clr,'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Salinity (psu)')
legend('show','location','best')
title('Final QC''d data')

fig3 = figure(3);clf
fig3.WindowState = 'maximized';
plot(gull.datetime_utc,gull.temperature,'.','Color',gull_clr,'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.temperature,'.','Color',north_clr,'DisplayName','North')
plot(south.datetime_utc,south.temperature,'.','Color',south_clr,'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Temperature (^oC)')
legend('show','location','best')
title('Final QC''d data')

fig4 = figure(4);clf
fig4.WindowState = 'maximized';
plot(gull.datetime_utc,gull.DOconc,'.','Color',gull_clr,'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.DOconc,'.','Color',north_clr,'DisplayName','North')
plot(south.datetime_utc,south.DOconc,'.','Color',south_clr,'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('DO conc (\mumol/L)')
legend('show','location','best')
title('Final QC''d data')

fig5 = figure(5);clf
fig5.WindowState = 'maximized';
plot(gull.datetime_utc,gull.pH,'.','Color',gull_clr,'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.pH,'.','Color',north_clr,'DisplayName','North')
plot(south.datetime_utc,south.pH,'.','Color',south_clr,'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('pH')
legend('show','location','best')
title('Final QC''d data')

fig6 = figure(6);clf
fig6.WindowState = 'maximized';
plot(gull.datetime_utc,gull.turbidity,'.','Color',gull_clr,'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.turbidity,'.','Color',north_clr,'DisplayName','North')
plot(south.datetime_utc,south.turbidity,'.','Color',south_clr,'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Turbidity (NTU)')
legend('show','location','best')
title('Final QC''d data')

fig7 = figure(7);clf
fig7.WindowState = 'maximized';
plot(gull.datetime_utc,gull.chla,'.','Color',gull_clr,'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.chla,'.','Color',north_clr,'DisplayName','North')
plot(south.datetime_utc,south.chla,'.','Color',south_clr,'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Chl a (RFU)')
legend('show','location','best')
title('Final QC''d data')

% %% Check influence of wind events
% 
% % Turbidity
% fig = figure;clf
% fig.WindowState = 'maximized';
% yyaxis left
% plot(era5Dat.datetime,era5Dat.wspd,'-')
% ylabel('Wind speed (m/s)')
% 
% yyaxis right
% plot(gull.datetime_utc,gull.turbidity,'.-','Color',gull_clr,'DisplayName','Gull')
% hold on
% plot(north.datetime_utc,north.turbidity,'.-','Color',north_clr,'DisplayName','North')
% plot(south.datetime_utc,south.turbidity,'.-','Color',south_clr,'DisplayName','South')
% xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('Turbidity (NTU)')
% legend('show','location','best')
% 
% title('Final QC''d data')
% 
% % Salinity
% fig = figure;clf
% fig.WindowState = 'maximized';
% yyaxis left
% plot(era5Dat.datetime,era5Dat.wspd,'-')
% ylabel('Wind speed (m/s)')
% 
% yyaxis right
% plot(gull.datetime_utc,gull.salinity,'.-','Color',gull_clr,'DisplayName','Gull')
% hold on
% plot(north.datetime_utc,north.salinity,'.-','Color',north_clr,'DisplayName','North')
% plot(south.datetime_utc,south.salinity,'.-','Color',south_clr,'DisplayName','South')
% xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% % ylabel('Chl a (RFU)')
% ylabel('S')
% legend('show','location','best')
% 
% title('Final QC''d data')


%====Save the plots========================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\open-water-platform\all-sites\Final QC']);

        saveas(fig1,'finalQC_depth.png')
        saveas(fig1,'finalQC_depth.fig')
        saveas(fig2,'finalQC_salinity.png')
        saveas(fig2,'finalQC_salinity.fig')
        saveas(fig3,'finalQC_temperature.png')
        saveas(fig3,'finalQC_temperature.fig')
        saveas(fig4,'finalQC_DOconc.png')
        saveas(fig4,'finalQC_DOconc.fig')
        saveas(fig5,'finalQC_pH.png')
        saveas(fig5,'finalQC_pH.fig')
        saveas(fig6,'finalQC_turbidity.png')
        saveas(fig6,'finalQC_turbidity.fig')
        saveas(fig7,'finalQC_chla.png')
        saveas(fig7,'finalQC_chla.fig')

        disp('Plots saved!')

    case 'No'
        disp('Plots not saved.')
end