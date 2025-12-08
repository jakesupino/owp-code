%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotDielSeasonalPatterns.m
% This script
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 5/23/2024
% Last updated: 7/11/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%==========================================================================
%   Import data
%==========================================================================
site = 'Gull'; % Input
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load('finalQC.mat')

%==========================================================================
%   Plot entire time series
%==========================================================================
% fig1 = figure(1);clf
% fig1.WindowState = 'maximized';
% % yyaxis left
% plot(finalQC.datetime_utc,finalQC.pH,'-')
% ylabel('pH (-)')
% yyaxis right
% plot(finalQC.datetime_utc,finalQC.DOsat,'-')
% yline(100,'k','linewidth',1.5)
% ylabel('DO sat (%)')
%%
%==========================================================================
%   Plot 1 week of data (in local time)
%==========================================================================
%====Determine day/night times=============================================
% Get day/night indices based on lat and lon and Hilary's indexDayNight function
lat = 39.08;
lon = -74.78;
time_in = finalQC.datetime_utc;
tol = 0;
UTCoffset = 0;  % Input datetime vector is in UTC
[dayind,nightind] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);

% Find the indices for when each day starts and stops
daystart = dayind(find(diff(dayind) > 1) + 1);
dayend = dayind(find(diff(dayind) > 1));

datetime_local = finalQC.datetime_utc;
datetime_local.TimeZone = 'America/New_York';

daystart_dt = datetime_local(daystart(1:end-1));
dayend_dt = datetime_local(dayend(2:end));

fig2 = figure(2);clf
t1 = tiledlayout(2,1,'TileSpacing','tight');
pos = [1 41 900 748];
fig2.Position = pos;

%====January 2023==========================================================
ind1 = 58506;
ind2 = ind1 + 1008;

% Find indices of closest matching datetimes to daystart and dayend times within window
dt1 = daystart_dt;
dt2 = datetime_local(ind1:ind2);
start_nearestind = interp1(dt2,1:length(dt2),dt1,'nearest');
start_nearestind = rmmissing(start_nearestind);

dt3 = dayend_dt;
dt4 = datetime_local(ind1:ind2);
end_nearestind = interp1(dt4,1:length(dt4),dt3,'nearest');
end_nearestind = rmmissing(end_nearestind);

nexttile
% Plot pH
yyaxis right
ymin = 8;
ymax = 8.6;
plot(datetime_local,finalQC.pH,'-')
ylim([ymin ymax])
ylabel('pH','FontSize',18)
% Plot DO
yyaxis left
ymin = 80;
ymax = 130;
x_points = [dt2(end_nearestind),dt2(end_nearestind),dt2(start_nearestind),dt2(start_nearestind)];
y_points = [repelem(ymin,7,1),repelem(ymax,7,2),repelem(ymin,7,1)];
for i = 1:length(start_nearestind)
    fill(x_points(i,:),y_points(i,:),rgb('lightgrey'),'LineStyle','none')
    hold on
end
plot(datetime_local(ind1:ind2),finalQC.DOsat(ind1:ind2),'-')
hold on
yline(100,'k','linewidth',1.5)
ylabel('DO sat (%)','FontSize',18)
xlim([datetime_local(ind1) datetime_local(ind2)])
ylim([ymin ymax])
box on

%====July 2023=============================================================
% Find indices of closest matching datetimes to daystart and dayend times within window
ind1 = 86405;
ind2 = ind1 + 1008;

dt1 = daystart_dt;
dt2 = datetime_local(ind1:ind2);
start_nearestind = interp1(dt2,1:length(dt2),dt1,'nearest');
start_nearestind = rmmissing(start_nearestind);

dt3 = dayend_dt;
dt4 = datetime_local(ind1:ind2);
end_nearestind = interp1(dt4,1:length(dt4),dt3,'nearest');
end_nearestind = rmmissing(end_nearestind);

nexttile
yyaxis right
ymin = 7;
ymax = 8;
plot(datetime_local,finalQC.pH,'-')
ylabel('pH','FontSize',18)
ylim([ymin ymax])
yyaxis left
ymin = 0;
ymax = 120;
x_points = [dt2(end_nearestind),dt2(end_nearestind),dt2(start_nearestind),dt2(start_nearestind)];
y_points = [repelem(ymin,7,1),repelem(ymax,7,2),repelem(ymin,7,1)];
for i = 1:length(start_nearestind)
    fill(x_points(i,:),y_points(i,:),rgb('lightgrey'),'LineStyle','none')
    hold on
end
plot(datetime_local(ind1:ind2),finalQC.DOsat(ind1:ind2),'-')
hold on
yline(100,'k','linewidth',1.5)
ylabel('DO sat (%)','FontSize',18)
xlim([datetime_local(ind1) datetime_local(ind2)])
ylim([ymin ymax])

% title(t1,num2str(site),'FontSize',28,'FontName','Helvetica','FontWeight','bold')

cd([rootpath,'figures\open-water-platform\all-sites'])
% Option to save figs
option = questdlg('Save figure?','Save Figure','Yes','No','Yes');
switch option
    case 'Yes'
        exportgraphics(t1,[site,'_dielPatterns.png'])
        saveas(t1,[site,'_dielPatterns.fig'])
        disp('Figure saved!')
    case 'No'
        disp('Figure not saved.')
end
%% Tiled plots of one site -- 3x2 format
site = 'North'; % Input
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat'])

%====Calculate [H+] and add to final QC table===========================
Hplus = 10.^(-finalQC.pH);
finalQC = addvars(finalQC,Hplus,'after','pH');

%====Do calculations====================================================
% (1) Compute monthly mean values, monthly means of daily maxima, and monthly means of daily minima
monAvg = groupsummary(finalQC,'datetime_utc','monthofyear','mean');
dielMax = groupsummary(finalQC,'datetime_utc','day','max');
dielMin = groupsummary(finalQC,'datetime_utc','day','min');
dielMax.day_datetime_utc = datetime(string(dielMax.day_datetime_utc));  % Convert day_datetime_utc variable to datetime
dielMin.day_datetime_utc = datetime(string(dielMin.day_datetime_utc));
dielMax_monAvg = groupsummary(dielMax,'day_datetime_utc','monthofyear','mean');
dielMin_monAvg = groupsummary(dielMin,'day_datetime_utc','monthofyear','mean');

% Convert monthly mean [H+] back into pH space
monAvg_pH = -log10(monAvg.mean_Hplus);
dielMax_monAvg_pH = -log10(dielMax_monAvg.mean_max_Hplus);
dielMin_monAvg_pH = -log10(dielMin_monAvg.mean_min_Hplus);

% (2a) Compute diel ranges
dielRange = groupsummary(finalQC,'datetime_utc','day','range');
dielRange.day_datetime_utc = datetime(string(dielRange.day_datetime_utc));

% Before computing monthly means of diel ranges, put pH in [H+] space
dielRange_Hplus = 10.^(-dielRange.range_pH);
dielRange.range_Hplus = dielRange_Hplus;

% (2b) Find monthly means and standard deviations of diel ranges
dielRange_monAvg = groupsummary(dielRange,'day_datetime_utc','monthofyear','mean');
% Convert monthly mean [H+] range back into pH space
% Simply take standard deviation from raw pH values -- see Jonathan Midwood's answer from https://www.researchgate.net/post/Is_it_possible_to_calculate_standard_deviation_for_pH
dielRange_monAvg_pH = -log10(dielRange_monAvg.mean_range_Hplus);
dielRange_monStd = groupsummary(dielRange,'day_datetime_utc','monthofyear','std');

%====Make the figure====================================================
fig3 = figure(3);clf
t3 = tiledlayout(3,2);
pos = [1 41 1536 748];
fig3.Position = pos;

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_temperature;
yu = dielMax_monAvg.mean_max_temperature;
yl = dielMin_monAvg.mean_min_temperature;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([0 30])
ylabel('Temperature (^oC)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_temperature;
y_sd = dielRange_monStd.std_range_temperature;
errorbar(x,y,y_sd,'.-','LineWidth',1);
% ylim([0 12])  % For site comparison
ylim([0 12])
ylabel('Diel range (^oC)')
box on

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_salinity;
yu = dielMax_monAvg.mean_max_salinity;
yl = dielMin_monAvg.mean_min_salinity;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([30 36])
ylabel('Salinity (psu)')
yyaxis right
ax = gca;
ax.YColor ="#0072BD";
yyaxis right
ax = gca;
ax.YColor = 'k';
ax.YTick = [];
box on

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg_pH;
yu = dielMax_monAvg_pH;
yl = dielMin_monAvg_pH;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([7.5 8.5])
ylabel('pH (-)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg_pH;
y_sd = dielRange_monStd.std_range_pH;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylabel('Diel range (-)')
box on

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_chla;
yu = dielMax_monAvg.mean_max_chla;
yl = dielMin_monAvg.mean_min_chla;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
% ylim([0 70]) % For site comparison
ylim([0 65])
ylabel('Chl a (RFU)')
yyaxis right
ax = gca;
ax.YColor ="#0072BD";
yyaxis right
ax = gca;
ax.YColor = 'k';
ax.YTick = [];
box on

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_DOsat;
yu = dielMax_monAvg.mean_max_DOsat;
yl = dielMin_monAvg.mean_min_DOsat;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([30 130])
ylabel('DO sat (%)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_DOsat;
y_sd = dielRange_monStd.std_range_DOsat;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylim([0 80])
ylabel('Diel range (%)')
box on

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_turbidity;
yu = dielMax_monAvg.mean_max_turbidity;
yl = dielMin_monAvg.mean_min_turbidity;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([0 800])
ylabel('Turbidity (NTU)')
ax = gca;
ax.YColor ="#0072BD";
yyaxis right
ax = gca;
ax.YColor = 'k';
ax.YTick = [];
box on

title(t3,num2str(site),'FontSize',32,'FontName','Helvetica','FontWeight','bold')
xlabel(t3,'Month','FontSize',20,'FontName','Helvetica')

cd([rootpath,'figures\open-water-platform\site-comparisons'])
% Option to save figs
option = questdlg('Save figure?','Save Figure','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'figures\open-water-platform\site-comparisons'])
        exportgraphics(t3,[site,'3x2.png'])
        saveas(t3,[site,'3x2.fig'])
        disp('Figure saved!')
    case 'No'
        disp('Figure not saved.')
end

%% Tiled plots of all sites
site = ["North","Gull","South"];

fig2 = figure(2);clf
t2 = tiledlayout(6,3,'TileSpacing','compact','TileIndexing','columnmajor');
fig2.WindowState = 'maximized';

for i = 1:length(site)
    cd([rootpath,'open-water-platform-data\',num2str(site(i)),'\cleaned\final-qc'])
    load([num2str(site(i)),'-cleaned.mat'])

    monAvg = groupsummary(finalQC,'datetime_utc','monthofyear','mean');
    dielMax = groupsummary(finalQC,'datetime_utc','day','max');
    dielMin = groupsummary(finalQC,'datetime_utc','day','min');
    dielMax.day_datetime_utc = datetime(string(dielMax.day_datetime_utc));  % Convert day_datetime_utc variable to datetime
    dielMin.day_datetime_utc = datetime(string(dielMin.day_datetime_utc));
    dielMax_monAvg = groupsummary(dielMax,'day_datetime_utc','monthofyear','mean');
    dielMin_monAvg = groupsummary(dielMin,'day_datetime_utc','monthofyear','mean');

    dielRange = groupsummary(finalQC,'datetime_utc','day','range');
    dielRange.day_datetime_utc = string(dielRange.day_datetime_utc);
    dielRange.day_datetime_utc = datetime(dielRange.day_datetime_utc);
    dielRange_monAvg = groupsummary(dielRange,'day_datetime_utc','monthofyear','mean');
    dielRange_monStd = groupsummary(dielRange,'day_datetime_utc','monthofyear','std');
    
    [ind val] = min(monAvg(:,5:end),[],1)
    [ind val] = max(monAvg(:,5:end),[],1)

    [ind val] = min(dielRange_monAvg(:,[7 9 10]),[],1)
    [ind val] = max(dielRange_monAvg(:,[7 9 10]),[],1)

    ax1 = nexttile;
    yyaxis left
    x = monAvg.monthofyear_datetime_utc;
    y = monAvg.mean_temperature;
    yu = dielMax_monAvg.mean_max_temperature;
    yl = dielMin_monAvg.mean_min_temperature;
    patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
    hold on
    plot(x,y,'.-');
    ylim([0 30])
    ylabel('Temperature (^oC)')
    yyaxis right
    x = dielRange_monAvg.monthofyear_day_datetime_utc;
    y = dielRange_monAvg.mean_range_temperature;
    y_sd = dielRange_monStd.std_range_temperature;
    errorbar(x,y,y_sd,'.-','LineWidth',1);
    ylim([0 12])
    ylabel('Diel range (^oC)')
    box on
    title(num2str(site(i)),'FontSize',20,'FontName','Helvetica','FontWeight','bold')

    ax2 = nexttile;
    yyaxis left
    x = monAvg.monthofyear_datetime_utc;
    y = monAvg.mean_pH;
    yu = dielMax_monAvg.mean_max_pH;
    yl = dielMin_monAvg.mean_min_pH;
    patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
    hold on
    plot(x,y,'.-');
    ylim([7.5 8.5])
    ylabel('pH (-)')
    yyaxis right
    x = dielRange_monAvg.monthofyear_day_datetime_utc;
    y = dielRange_monAvg.mean_range_pH;
    y_sd = dielRange_monStd.std_range_pH;
    errorbar(x,y,y_sd,'.-','LineWidth',1);
    ylim([0 0.8])
    ylabel('Diel range (-)')
    box on

    ax3 = nexttile;
    yyaxis left
    x = monAvg.monthofyear_datetime_utc;
    y = monAvg.mean_DOsat;
    yu = dielMax_monAvg.mean_max_DOsat;
    yl = dielMin_monAvg.mean_min_DOsat;
    patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
    hold on
    plot(x,y,'.-');
    % ylim([80 380])
    % ylabel('DO conc (\mumol/L)')
    ylim([30 130])
    ylabel('DO sat (%)')
    yyaxis right
    x = dielRange_monAvg.monthofyear_day_datetime_utc;
    y = dielRange_monAvg.mean_range_DOsat;
    y_sd = dielRange_monStd.std_range_DOsat;
    errorbar(x,y,y_sd,'.-','LineWidth',1);
    % ylim([0 200])
    % ylabel('Diel range (\mumol/L)')
    ylim([0 80])
    ylabel('Diel range (%)')
    box on

    ax4 = nexttile;
    yyaxis left
    x = monAvg.monthofyear_datetime_utc;
    y = monAvg.mean_salinity;
    yu = dielMax_monAvg.mean_max_salinity;
    yl = dielMin_monAvg.mean_min_salinity;
    patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
    hold on
    plot(x,y,'.-');
    ylim([30 36])
    ylabel('Salinity (psu)')
    yyaxis right
    ax = gca;
    ax.YColor ="#0072BD";
    yyaxis right
    ax = gca;
    ax.YColor = 'k';
    ax.YTick = [];
    box on

    ax5 = nexttile;
    yyaxis left
    x = monAvg.monthofyear_datetime_utc;
    y = monAvg.mean_chla;
    yu = dielMax_monAvg.mean_max_chla;
    yl = dielMin_monAvg.mean_min_chla;
    patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
    hold on
    plot(x,y,'.-');
    ylim([0 70])
    ylabel('Chl a (RFU)')
    % xlabel('Month')
    yyaxis right
    ax = gca;
    ax.YColor ="#0072BD";
    yyaxis right
    ax = gca;
    ax.YColor = 'k';
    ax.YTick = [];
    box on

    ax6 = nexttile;
    yyaxis left
    x = monAvg.monthofyear_datetime_utc;
    y = monAvg.mean_turbidity;
    yu = dielMax_monAvg.mean_max_turbidity;
    yl = dielMin_monAvg.mean_min_turbidity;
    patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
    hold on
    plot(x,y,'.-');
    ylim([0 800])
    ylabel('Turbidity (NTU)')
    xlabel('Month')
    ax = gca;
    ax.YColor ="#0072BD";
    yyaxis right
    ax = gca;
    ax.YColor = 'k';
    ax.YTick = [];
    box on

    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x')
end

% Option to save figs
option = questdlg('Save figure?','Save Figure','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'figures\open-water-platform\site-comparisons'])
        saveas(fig2,'AllSites_AllParams.png')
        saveas(fig2,'AllSites_AllParams.fig')
        disp('Figure saved!')
    case 'No'
        disp('Figure not saved.')
end

%% T, pH, DO -- all sites together
fig = figure(1);clf
tiledlayout(3,3)
fig.WindowState = 'maximized';

% North
site = 'North';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat'])

monAvg = groupsummary(finalQC,'datetime_utc','monthofyear','mean');
dielMax = groupsummary(finalQC,'datetime_utc','day','max');
dielMin = groupsummary(finalQC,'datetime_utc','day','min');
dielMax.day_datetime_utc = datetime(string(dielMax.day_datetime_utc));  % Convert day_datetime_utc variable to datetime
dielMin.day_datetime_utc = datetime(string(dielMin.day_datetime_utc));
dielMax_monAvg = groupsummary(dielMax,'day_datetime_utc','monthofyear','mean');
dielMin_monAvg = groupsummary(dielMin,'day_datetime_utc','monthofyear','mean');

dielRange = groupsummary(finalQC,'datetime_utc','day','range');
dielRange.day_datetime_utc = string(dielRange.day_datetime_utc);
dielRange.day_datetime_utc = datetime(dielRange.day_datetime_utc);
dielRange_monAvg = groupsummary(dielRange,'day_datetime_utc','monthofyear','mean');
dielRange_monStd = groupsummary(dielRange,'day_datetime_utc','monthofyear','std');

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_temperature;
yu = dielMax_monAvg.mean_max_temperature;
yl = dielMin_monAvg.mean_min_temperature;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([0 30])
ylabel('Monthly means (^oC)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_temperature;
y_sd = dielRange_monStd.std_range_temperature;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylim([0 12])
ylabel('Diel range (^oC)')
% xlabel('Month')
box on
title('Temperature','FontSize',20,'FontName','Helvetica','FontWeight','bold')

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_pH;
yu = dielMax_monAvg.mean_max_pH;
yl = dielMin_monAvg.mean_min_pH;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([7.5 8.5])
ylabel('Monthly means (-)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_pH;
y_sd = dielRange_monStd.std_range_pH;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylim([0 0.8])
ylabel('Diel range (-)')
% xlabel('Month')
box on
title('pH','FontSize',20,'FontName','Helvetica','FontWeight','bold')

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_DOsat;
yu = dielMax_monAvg.mean_max_DOsat;
yl = dielMin_monAvg.mean_min_DOsat;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([30 130])
ylabel('Monthly means (%)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_DOsat;
y_sd = dielRange_monStd.std_range_DOsat;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylim([0 80])
ylabel('Diel range (%)')
% xlabel('Month')
box on
title('DO % Saturation','FontSize',20,'FontName','Helvetica','FontWeight','bold')

% Gull
site = 'Gull';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat'])

monAvg = groupsummary(finalQC,'datetime_utc','monthofyear','mean');
dielMax = groupsummary(finalQC,'datetime_utc','day','max');
dielMin = groupsummary(finalQC,'datetime_utc','day','min');
dielMax.day_datetime_utc = datetime(string(dielMax.day_datetime_utc));  % Convert day_datetime_utc variable to datetime
dielMin.day_datetime_utc = datetime(string(dielMin.day_datetime_utc));
dielMax_monAvg = groupsummary(dielMax,'day_datetime_utc','monthofyear','mean');
dielMin_monAvg = groupsummary(dielMin,'day_datetime_utc','monthofyear','mean');

dielRange = groupsummary(finalQC,'datetime_utc','day','range');
dielRange.day_datetime_utc = string(dielRange.day_datetime_utc);
dielRange.day_datetime_utc = datetime(dielRange.day_datetime_utc);
dielRange_monAvg = groupsummary(dielRange,'day_datetime_utc','monthofyear','mean');
dielRange_monStd = groupsummary(dielRange,'day_datetime_utc','monthofyear','std');

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_temperature;
yu = dielMax_monAvg.mean_max_temperature;
yl = dielMin_monAvg.mean_min_temperature;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([0 30])
ylabel('Monthly means (^oC)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_temperature;
y_sd = dielRange_monStd.std_range_temperature;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylim([0 12])
ylabel('Diel range (^oC)')
% xlabel('Month')
box on

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_pH;
yu = dielMax_monAvg.mean_max_pH;
yl = dielMin_monAvg.mean_min_pH;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([7.5 8.5])
ylabel('Monthly means (-)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_pH;
y_sd = dielRange_monStd.std_range_pH;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylim([0 0.8])
ylabel('Diel range (-)')
% xlabel('Month')
box on

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_DOsat;
yu = dielMax_monAvg.mean_max_DOsat;
yl = dielMin_monAvg.mean_min_DOsat;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([30 130])
ylabel('Monthly means (%)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_DOsat;
y_sd = dielRange_monStd.std_range_DOsat;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylim([0 80])
ylabel('Diel  range (%)')
% xlabel('Month')
box on

% South
site = 'South';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat'])

monAvg = groupsummary(finalQC,'datetime_utc','monthofyear','mean');
dielMax = groupsummary(finalQC,'datetime_utc','day','max');
dielMin = groupsummary(finalQC,'datetime_utc','day','min');
dielMax.day_datetime_utc = datetime(string(dielMax.day_datetime_utc));  % Convert day_datetime_utc variable to datetime
dielMin.day_datetime_utc = datetime(string(dielMin.day_datetime_utc));
dielMax_monAvg = groupsummary(dielMax,'day_datetime_utc','monthofyear','mean');
dielMin_monAvg = groupsummary(dielMin,'day_datetime_utc','monthofyear','mean');

dielRange = groupsummary(finalQC,'datetime_utc','day','range');
dielRange.day_datetime_utc = string(dielRange.day_datetime_utc);
dielRange.day_datetime_utc = datetime(dielRange.day_datetime_utc);
dielRange_monAvg = groupsummary(dielRange,'day_datetime_utc','monthofyear','mean');
dielRange_monStd = groupsummary(dielRange,'day_datetime_utc','monthofyear','std');

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_temperature;
yu = dielMax_monAvg.mean_max_temperature;
yl = dielMin_monAvg.mean_min_temperature;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([0 30])
ylabel('Monthly means (^oC)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_temperature;
y_sd = dielRange_monStd.std_range_temperature;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylim([0 12])
ylabel('Diel range (^oC)')
xlabel('Month')
box on

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_pH;
yu = dielMax_monAvg.mean_max_pH;
yl = dielMin_monAvg.mean_min_pH;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([7.5 8.5])
ylabel('Monthly means (-)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_pH;
y_sd = dielRange_monStd.std_range_pH;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylim([0 0.8])
ylabel('Diel range (-)')
xlabel('Month')
box on

nexttile
yyaxis left
x = monAvg.monthofyear_datetime_utc;
y = monAvg.mean_DOsat;
yu = dielMax_monAvg.mean_max_DOsat;
yl = dielMin_monAvg.mean_min_DOsat;
patch([x; flip(x)] , [yl; flip(yu)], [.6 .7 .8])
hold on
plot(x,y,'.-');
ylim([30 130])
ylabel('Monthly means (%)')
yyaxis right
x = dielRange_monAvg.monthofyear_day_datetime_utc;
y = dielRange_monAvg.mean_range_DOsat;
y_sd = dielRange_monStd.std_range_DOsat;
errorbar(x,y,y_sd,'.-','LineWidth',1);
ylim([0 80])
ylabel('Diel range (%)')
xlabel('Month')
box on

cd([rootpath,'figures\open-water-platform\site-comparisons'])