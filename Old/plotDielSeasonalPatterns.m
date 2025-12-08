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
% yyaxis left
% plot(finalQC.datetime_utc,finalQC.pH,'-')
% ylabel('pH (-)')
% yyaxis right
% plot(finalQC.datetime_utc,finalQC.DOsat,'-')
% ylabel('DO sat (%)')

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

% daystart_dt = time_in(daystart(1:end-1));
% dayend_dt = time_in(dayend(2:end));
daystart_dt = datetime_local(daystart(1:end-1));
dayend_dt = datetime_local(dayend(2:end));

fig2 = figure(2);clf
t1 = tiledlayout(2,1);
pos = [1 41 900 748];
fig2.Position = pos;

%====January 2023==========================================================
ind1 = 58506;
ind2 = ind1 + 1008;

% Find indices of closest matching datetimes to daystart and dayend times within window
dt1 = daystart_dt;
% dt2 = finalQC.datetime_utc(ind1:ind2);
dt2 = datetime_local(ind1:ind2);
start_nearestind = interp1(dt2,1:length(dt2),dt1,'nearest');
start_nearestind = rmmissing(start_nearestind);

dt3 = dayend_dt;
% dt4 = finalQC.datetime_utc(ind1:ind2);
dt4 = datetime_local(ind1:ind2);
end_nearestind = interp1(dt4,1:length(dt4),dt3,'nearest');
end_nearestind = rmmissing(end_nearestind);

ymin = 8;
ymax = 8.7;

x_points = [dt2(end_nearestind),dt2(end_nearestind),dt2(start_nearestind),dt2(start_nearestind)];
y_points = [repelem(ymin,7,1),repelem(ymax,7,2),repelem(ymin,7,1)];

nexttile
yyaxis left
for i = 1:length(start_nearestind)
    fill(x_points(i,:),y_points(i,:),rgb('lightgrey'),'LineStyle','none')
    hold on
end
% plot(finalQC.datetime_utc,finalQC.pH,'-')
plot(datetime_local,finalQC.pH,'-')
ylim([ymin ymax])
ylabel('pH (-)','FontSize',18)
yyaxis right
% plot(finalQC.datetime_utc(ind1:ind2),finalQC.DOsat(ind1:ind2),'-')
plot(datetime_local(ind1:ind2),finalQC.DOsat(ind1:ind2),'-')
ylabel('DO sat (%)','FontSize',18)
% xlim([finalQC.datetime_utc(ind1) finalQC.datetime_utc(ind2)])
xlim([datetime_local(ind1) datetime_local(ind2)])
ylim([80 140])
box on

%====July 2023=============================================================
% Find indices of closest matching datetimes to daystart and dayend times within window
ind1 = 86405;
ind2 = ind1 + 1008;

dt1 = daystart_dt;
% dt2 = finalQC.datetime_utc(ind1:ind2);
dt2 = datetime_local(ind1:ind2);
start_nearestind = interp1(dt2,1:length(dt2),dt1,'nearest');
start_nearestind = rmmissing(start_nearestind);

dt3 = dayend_dt;
% dt4 = finalQC.datetime_utc(ind1:ind2);
dt4 = datetime_local(ind1:ind2);
end_nearestind = interp1(dt4,1:length(dt4),dt3,'nearest');
end_nearestind = rmmissing(end_nearestind);

ymin = 7.3;
ymax = 8;

x_points = [dt2(end_nearestind),dt2(end_nearestind),dt2(start_nearestind),dt2(start_nearestind)];
y_points = [repelem(ymin,7,1),repelem(ymax,7,2),repelem(ymin,7,1)];

nexttile
yyaxis left
for i = 1:length(start_nearestind)
    fill(x_points(i,:),y_points(i,:),rgb('lightgrey'),'LineStyle','none')
    hold on
end
% plot(finalQC.datetime_utc,finalQC.pH,'-')
plot(datetime_local,finalQC.pH,'-')
ylabel('pH (-)','FontSize',18)
ylim([ymin ymax])
yyaxis right
% plot(finalQC.datetime_utc(ind1:ind2),finalQC.DOsat(ind1:ind2),'-')
plot(datetime_local(ind1:ind2),finalQC.DOsat(ind1:ind2),'-')
ylabel('DO sat (%)','FontSize',18)
% xlim([finalQC.datetime_utc(ind1) finalQC.datetime_utc(ind2)])
xlim([datetime_local(ind1) datetime_local(ind2)])
ylim([0 120])

title(t1,num2str(site),'FontSize',28,'FontName','Helvetica','FontWeight','bold')

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