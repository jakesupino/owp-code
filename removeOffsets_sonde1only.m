%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% removeOffsets_sonde1only.m
% This script removes time offsets from UTC, and erroneous depth offsets.
% Refer to the Collaborative Lab Notebook, which documents which deployments
% have time and/or depth offsets.
%
% Run this script for deployments that only have data for Sonde 1 (BC)!!
%
% AUTHOR:
% Emily Chua
%
% DATE:
% 10/30/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

fig = uifigure;
site = uiconfirm(fig,"Select the platform","Site selection","Options",["gull","north","south"]);
close(fig)

% Load merged data for the platform
cd([rootpath,'open-water-platform-data\',site,'\original\merged'])
load(['alldeps-',site,'.mat'])

means1 = grpstats(sonde1_all,"deployment",{"mean"});

% For each sonde, find the mean depth across all deployments that aren't erroneously high
switch site
    case 'gull'
        alldepths1_avg = mean(means1.mean_depth([1 2 3:8 10 11 13])); % BC: Dep# 1, 2, 5-10, 12, 13, 15
        alldepths2_avg = mean(means2.mean_depth([1:6 8:10 13])); % ERDC: Dep# 1, 2, 5-8, 10-12, 15
    case 'north'
        alldepths1_avg = mean(means1.mean_depth([1:6 8 10 11])); % BC: Dep# 2, 6-10, 12, 14, 15
        alldepths2_avg = mean(means2.mean_depth([1:4 7 8 10 11])); % ERDC: Dep# 2, 6-8, 11, 12, 14, 15
    case 'south'
        alldepths1_avg = mean(means1.mean_depth([1:7 9:11])); % BC: Dep# 1, 2, 4-8, 10-12
end

%===Read in sonde data for a specific deployment===========================
cd([rootpath,'open-water-platform-data\',site,'\original\deployments'])

[fileName,dataPath] = uigetfile('*.mat');

load(fileName);

depNum = sonde1.deployment(1);

%===Read in USGS data======================================================
paramNames = ["agency","site_no","datetime_local","timezone","tidal_elev","qual-code"];
usgs = readtable('G:\My Drive\Postdoc\SMIIL\usgs-data\tidal-elev.txt','TextType','string');
usgs.Properties.VariableNames = paramNames;

usgs.datetime_local.TimeZone = 'America/New_York';
datetime_utc = table(datetime(usgs.datetime_local,'TimeZone','utc'),'VariableNames',"datetime_utc");
usgs = [datetime_utc,usgs];

% Convert [ft] to [m]
tidal_elev = usgs.tidal_elev/3.281;
usgs.tidal_elev = tidal_elev;

paramUnits = ["","","","","","m",""];
usgs.Properties.VariableUnits = paramUnits;

%===Plot raw sonde and USGS data===========================================
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
FontSize = 12;
LineWidth = 1;

fig1 = figure(1);
fig1.WindowState = 'maximized';
h3 = plot(usgs.datetime_utc,usgs.tidal_elev,'k');
hold on
h1 = plot(sonde1.datetime_utc,sonde1.depth,'Color',red);
    legend([h1 h3],'BC','USGS')
hold off
xlim([min(sonde1.datetime_utc) max(sonde1.datetime_utc)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde1.deployment(1)),' - Original'])
set(gca,'FontSize',FontSize,'LineWidth',LineWidth)
grid on

disp('Press enter to remove depth offset')
pause

%===Remove depth offset============================================
depth1_avg = mean(sonde1.depth,"omitnan");
% delta_depth1 = abs(depth1_avg - alldepths1_avg);
delta_depth1 = depth1_avg - alldepths1_avg;
depth1_adj = sonde1.depth - delta_depth1;

%===Adjust pressure data===========================================
p1_adj = sonde1.density.*1000*9.81.*depth1_adj/6894.76;

%===Fix time offset(s) from UTC============================================
% Find peaks in USGS data
dups = find(diff(usgs.datetime_utc)<minutes(6));    % Find duplicate USGS timepoints
usgs(dups,:) = [];                                  % Remove duplicate USGS timepoints
smoothed0 = smoothdata(usgs.tidal_elev);
[pks0,locs0] = findpeaks(smoothed0,'MinPeakDistance',115);

% Find peaks in sonde data; smooth first to remove local peaks
smoothed1 = smoothdata(sonde1.depth);
% Change this value depending on sampling interval (6, 10, or 12 min)
minPkDist = 55;   % 12 min
% minPkDist = 65;    % 10 min
% minPkDist = 115;   % 6 min
[pks1,locs1] = findpeaks(smoothed1,'MinPeakDistance',minPkDist);

% figure(3),clf
% plot(sonde1.datetime_utc,sonde1.depth)
% hold on
% plot(sonde1.datetime_utc,smoothed1)
% hold off
% legend('Original','Smoothed')

% Plot peaks in data to double-check
fig2=figure(2);
fig2.WindowState = 'maximized';
h0 = plot(usgs.datetime_utc,usgs.tidal_elev,'k');
hold on
plot(usgs.datetime_utc(locs0),usgs.tidal_elev(locs0),'ok')
h1 = plot(sonde1.datetime_utc,sonde1.depth,'Color',red);
plot(sonde1.datetime_utc(locs1),sonde1.depth(locs1),'o','Color',red)
hold off
legend([h0 h1],'USGS','BC')
xlim([min(sonde1.datetime_utc) max(sonde1.datetime_utc)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde1.deployment(1)),' - Identified Peaks'])

disp('Press enter to adjust times')
pause

% Find time offset between third sonde peak and corresponding USGS peak
diff1 = usgs.datetime_utc(locs0) - sonde1.datetime_utc(locs1(3));
ind1 = find(diff1==max(diff1(diff1<0)));
delta_t1 = diff1(ind1);
if abs(delta_t1) > hours(12)
    delta_t1 = diff1(ind1+1);
else
    % Do nothing
end

% Adjust UTC time
datetime_utc1_adj = sonde1.datetime_utc + delta_t1;

% Adjust local time
datetime_local1_adj = sonde1.datetime_local + delta_t1;

%===Plot adjusted data=============================================
fig3 = figure(3);
fig3.WindowState = 'maximized';
h3 = plot(usgs.datetime_utc,usgs.tidal_elev,'k');
hold on
h1 = plot(datetime_utc1_adj,depth1_adj,'Color',red);
hold off
legend([h1 h3],'BC','USGS')
xlim([min(datetime_utc1_adj) max(datetime_utc1_adj)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde1.deployment(1)),' - Adjusted Depth and Time'])
set(gca,'FontSize',FontSize,'LineWidth',LineWidth)

pause

%====Save adjusted data====================================================
% Save new tables with columns for UTC & local time, depth, and p replaced
% with adjusted data
sonde1.datetime_utc = datetime_utc1_adj;
sonde1.datetime_local = datetime_local1_adj;
sonde1.depth = depth1_adj;
sonde1.p = p1_adj;

%====Save created tables in .mat files=====================================
saveFilePath = ['open-water-platform-data\',site,'\adjusted\deployments'];

option = questdlg(['Save .mat file in ',saveFilePath,'?'],'Save File','Y','N','Y');
switch option
    case 'Y'
        cd([rootpath,saveFilePath])
        save(fileName,"sonde1","delta_depth1","delta_t1")
        disp('File saved!')
    case 'N'
        disp('File not saved.')
end