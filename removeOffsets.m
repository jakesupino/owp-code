%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% removeOffsets.m
% This script removes time offsets from UTC, and erroneous depth offsets.
% Refer to the Collaborative Lab Notebook, which documents which deployments
% have time and/or depth offsets.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% 10/12/2023
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
means2 = grpstats(sonde2_all,"deployment",{"mean"});

% For each sonde, find the mean depth across all deployments that aren't
% erroneously high, excluding zero values.
switch site
    case 'gull'
        alldepths1_avg = mean(nonzeros(means1.mean_depth([1 2 3:8 10 11 13]))); % BC: Dep# 1, 2, 5-10, 12, 13, 15
        alldepths2_avg = mean(nonzeros(means2.mean_depth([1:6 8:10 13]))); % ERDC: Dep# 1, 2, 5-8, 10-12, 15
    case 'north'
        alldepths1_avg = mean(nonzeros(means1.mean_depth([1:6 8 10 11]))); % BC: Dep# 2, 6-10, 12, 14, 15
        alldepths2_avg = mean(nonzeros(means2.mean_depth([1:4 7 8 10 11]))); % ERDC: Dep# 2, 6-8, 11, 12, 14, 15
    case 'south'
        alldepths1_avg = mean(nonzeros(means1.mean_depth([1:7 9:11]))); % BC: Dep# 1, 2, 4-8, 10-12
        alldepths2_avg = mean(nonzeros(means2.mean_depth([1:6 9 10 12]))); % ERDC: Dep# 1, 2, 5, 7-9, 12, 13, 15
end

%===Read in sonde data for a specific deployment===========================
cd([rootpath,'open-water-platform-data\',site,'\original\deployments'])

[fileName,dataPath] = uigetfile('*.mat');

load(fileName);

depNum = sonde1.deployment(1);

%===Read in USGS data======================================================
paramNames = ["agency","site_no","datetime_local","timezone","tidal_elev","qual-code"];
usgs = readtable([rootpath,'usgs-data\tidal-elev.txt'],'TextType','string');
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
if strcmp(site,'south') == 1 && (depNum == 4 || depNum == 6)
    % Don't plot sonde2 data for skipped deployments
    legend([h1 h3],'BC','USGS')
else
    h2 = plot(sonde2.datetime_utc,sonde2.depth,'Color',blue);
    legend([h1 h2 h3],'BC','ERDC','USGS')
end
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
depth1_avg = mean(nonzeros(sonde1.depth),"omitnan");
depth2_avg = mean(nonzeros(sonde2.depth),"omitnan");
delta_depth1 = depth1_avg - alldepths1_avg;
delta_depth2 = depth2_avg - alldepths2_avg;
depth1_adj = sonde1.depth - delta_depth1;
depth2_adj = sonde2.depth - delta_depth2;

%===Adjust pressure data===========================================
p1_adj = sonde1.density.*1000*9.81.*depth1_adj/6894.76;
p2_adj = sonde2.density.*1000*9.81.*depth2_adj/6894.76;

%===Fix time offset(s) from UTC============================================
% Some deployments have weird time offsets between the BC and ERDC sondes, assessed visually
% South Deployment 5
% Gull Deployment 8
if strcmp(site,'south') && depNum == 5
    t_offset = sonde2.datetime_utc(1) - sonde1.datetime_utc(1);
    sonde2.datetime_utc = sonde2.datetime_utc - t_offset;
elseif strcmp(site,'gull') && depNum == 8
    t_offset = abs(sonde2.datetime_utc(1) - sonde1.datetime_utc(1));
    sonde1.datetime_utc = sonde1.datetime_utc - t_offset;
% elseif strcmp
else
    % Do nothing
end

% Find peaks in USGS data
dups = find(diff(usgs.datetime_utc)<minutes(6));    % Find duplicate USGS timepoints
usgs(dups,:) = [];                                  % Remove duplicate USGS timepoints
smoothed0 = smoothdata(usgs.tidal_elev);
[pks0,locs0] = findpeaks(smoothed0,'MinPeakDistance',115);

% Find peaks in sonde data; smooth first to remove local peaks
smoothed1 = smoothdata(sonde1.depth);
smoothed2 = smoothdata(sonde2.depth);
% Change this value depending on sampling interval (6, 10, or 12 min)
minPkDist = 55;   % 12 min
% minPkDist = 65;    % 10 min
% minPkDist = 115;   % 6 min
[pks1,locs1] = findpeaks(smoothed1,'MinPeakDistance',minPkDist);
[pks2,locs2] = findpeaks(smoothed2,'MinPeakDistance',minPkDist);

% figure(3),clf
% plot(sonde1.datetime_utc,sonde1.depth)
% hold on
% plot(sonde1.datetime_utc,smoothed1)
% hold off
% legend('Original','Smoothed')
% 
% figure(4),clf
% plot(sonde2.datetime_utc,sonde2.depth)
% hold on
% plot(sonde2.datetime_utc,smoothed2)
% hold off
% legend('Original','Smoothed')

% Plot peaks in data to double-check
fig2= figure(2);
fig2.WindowState = 'maximized';
h0 = plot(usgs.datetime_utc,usgs.tidal_elev,'k');
hold on
plot(usgs.datetime_utc(locs0),usgs.tidal_elev(locs0),'ok')
h1 = plot(sonde1.datetime_utc,sonde1.depth,'Color',red);
plot(sonde1.datetime_utc(locs1),sonde1.depth(locs1),'o','Color',red)
h2 = plot(sonde2.datetime_utc,sonde2.depth,'Color',blue);
plot(sonde2.datetime_utc(locs2),sonde2.depth(locs2),'o','Color',blue)
hold off
legend([h0 h1 h2],'USGS','BC','ERDC')
xlim([min(sonde1.datetime_utc) max(sonde1.datetime_utc)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde1.deployment(1)),' - Identified Peaks'])
grid on

disp('Press enter to adjust times')
pause

% Find time offset between third sonde peak and corresponding USGS peak
diff1 = usgs.datetime_utc(locs0) - sonde1.datetime_utc(locs1(3));
ind1 = find(diff1==max(diff1(diff1<0)));
delta_t1 = diff1(ind1);
if abs(delta_t1) > hours(11.5)
    delta_t1 = diff1(ind1+1);
else
    % Do nothing
end

diff2 = usgs.datetime_utc(locs0) - sonde2.datetime_utc(locs2(3));
ind2 = find(diff2==max(diff2(diff2<0)));
delta_t2 = diff2(ind2);
if abs(delta_t2) > hours(11.5)
    delta_t2 = diff2(ind2+1);
else
    % Do nothing
end

% Adjust UTC time
datetime_utc1_adj = sonde1.datetime_utc + delta_t1;
datetime_utc2_adj = sonde2.datetime_utc + delta_t2;

% Adjust local time
datetime_local1_adj = sonde1.datetime_local + delta_t1;
datetime_local2_adj = sonde2.datetime_local + delta_t2;

%% Adjust for Daylight Saving Time
dt1 = diff(datetime_utc1_adj);
dst_end1 = find(dt1 > hours(1));
dst_start1 = find(dt1 < 0);
dt2 = diff(datetime_utc2_adj);
dst_end2 = find(dt2 > hours(1));
dst_start2 = find(dt2 < 0);

if ~isempty(dst_end1)
    % For times after DST ends, shift backwards by 1 hour
    datetime_utc1_adj(dst_end1+1:end) = datetime_utc1_adj(dst_end1+1:end) - hours(1);
elseif ~isempty(dst_start1)
    % For times after DST starts, shift forwards by 1 hour
    datetime_utc1_adj(dst_start1+1:end) = datetime_utc1_adj(dst_start1+1:end) + hours(1);
end

if ~isempty(dst_end2)
    % For times after DST ends, shift backwards by 1 hour
    datetime_utc2_adj(dst_end2+1:end) = datetime_utc2_adj(dst_end2+1:end) - hours(1);
elseif ~isempty(dst_start2)
    % For times after DST starts, shift forwards by 1 hour
    datetime_utc2_adj(dst_start2+1:end) = datetime_utc2_adj(dst_start2+1:end) + hours(1);
end


%===Plot adjusted data=============================================
fig3 = figure(4);
fig3.WindowState = 'maximized';
h3 = plot(usgs.datetime_utc,usgs.tidal_elev,'.k');
hold on
h1 = plot(datetime_utc1_adj,depth1_adj,'.','Color',red);
h2 = plot(datetime_utc2_adj,depth2_adj,'.','Color',blue);
hold off
legend([h1 h2 h3],'BC','ERDC','USGS')
xlim([min(datetime_utc1_adj) max(datetime_utc1_adj)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde1.deployment(1)),' - Adjusted Depth and Time'])
set(gca,'FontSize',FontSize,'LineWidth',LineWidth)
grid on

pause

%====Save adjusted data====================================================
% Save new tables with columns for UTC & local time, depth, and p replaced
% with adjusted data
sonde1.datetime_utc = datetime_utc1_adj;
sonde1.datetime_local = datetime_local1_adj;
sonde1.depth = depth1_adj;
sonde1.p = p1_adj;

sonde2.datetime_utc = datetime_utc2_adj;
sonde2.datetime_local = datetime_local2_adj;
sonde2.depth = depth2_adj;
sonde2.p = p2_adj;

%====Save created tables in .mat files=====================================
saveFilePath = ['open-water-platform-data\',site,'\adjusted\deployments'];

option = questdlg(['Save .mat file in ',saveFilePath,'?'],'Save File','Y','N','Y');
switch option
    case 'Y'
        cd([rootpath,saveFilePath])
        save(fileName,"sonde1","sonde2","delta_depth1","delta_depth2","delta_t1","delta_t2")
        disp('File saved!')
    case 'N'
        disp('File not saved.')
end