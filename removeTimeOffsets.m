%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% removeTimeOffsets.m
% This script removes time offsets from UTC. Refer to the Collaborative Lab Notebook,
% which documents which deployments have time offsets.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 10/12/2023
% Last amended: 6/26/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

fig = uifigure;
site = uiconfirm(fig,"Select the platform","Site selection","Options",["gull","north","south"]);
close(fig)

%===Read in sonde data for a specific deployment===========================
cd([rootpath,'open-water-platform-data\',site,'\original\deployments'])

[fileName,dataPath] = uigetfile('*.mat');

load(fileName);

depNum = sonde1.deployment(1);

%===Read in USGS data======================================================
paramNames = ["agency","site_no","datetime_local","timezone","tidal_elev","qual-code"];
usgs = readtable([rootpath,'physical-data\usgs\tidal-elev.txt'],'TextType','string');
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
h3 = plot(usgs.datetime_utc,usgs.tidal_elev,'.k');
hold on
h1 = plot(sonde1.datetime_utc,sonde1.depth,'.','Color',red);
if strcmp(site,'south') == 1 && (depNum == 4 || depNum == 6)
    % Don't plot sonde2 data for skipped deployments
    legend([h1 h3],'BC','USGS')
else
    h2 = plot(sonde2.datetime_utc,sonde2.depth,'.','Color',blue);
    legend([h1 h2 h3],'BC','ERDC','USGS')
end
hold off
xlim([min(sonde1.datetime_utc) max(sonde1.datetime_utc)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde1.deployment(1)),' - Original'])
set(gca,'FontSize',FontSize,'LineWidth',LineWidth)
grid on

%===Fix time offset(s) from UTC============================================
% Some deployments have weird time offsets between the BC and ERDC sondes, assessed visually
% Gull Deployment 8
% North Deployment 7
% South Deployment 5
if strcmp(site,'gull') && depNum == 8
    t_offset = abs(sonde2.datetime_utc(1) - sonde1.datetime_utc(1));
    sonde1.datetime_utc = sonde1.datetime_utc - t_offset;
elseif strcmp(site,'north') && depNum == 7
    t_offset = abs(sonde2.datetime_utc(1) - sonde1.datetime_utc(1));
    sonde2.datetime_utc = sonde2.datetime_utc - t_offset;
elseif strcmp(site,'south') && depNum == 5
    t_offset = sonde2.datetime_utc(1) - sonde1.datetime_utc(1);
    sonde2.datetime_utc = sonde2.datetime_utc - t_offset;
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
minPkDist = 55;   % 12 min (use for all other deployments)
% minPkDist = 65;    % 10 min
% minPkDist = 115;   % 6 min (use for Gull Deployment 1!!)
[pks1,locs1] = findpeaks(smoothed1,'MinPeakDistance',minPkDist);
[pks2,locs2] = findpeaks(smoothed2,'MinPeakDistance',minPkDist);
% Use next two lines for Gull Deployment 1 only!!
% [pks1,locs1] = findpeaks(sonde1.depth,'MinPeakDistance',minPkDist);
% [pks2,locs2] = findpeaks(sonde1.depth,'MinPeakDistance',minPkDist);

% Plot peaks in data to double-check
fig2 = figure(2);
fig2.WindowState = 'maximized';
h0 = plot(usgs.datetime_utc,usgs.tidal_elev,'.k');
hold on
plot(usgs.datetime_utc(locs0),usgs.tidal_elev(locs0),'ok')
h1 = plot(sonde1.datetime_utc,sonde1.depth,'.','Color',red);
plot(sonde1.datetime_utc(locs1),sonde1.depth(locs1),'o','Color',red)
h2 = plot(sonde2.datetime_utc,sonde2.depth,'.','Color',blue);
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

% Adjust for Daylight Saving Time
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
h3 = plot(usgs.datetime_utc,usgs.tidal_elev,'.k','MarkerSize',12);
hold on
h1 = plot(datetime_utc1_adj,sonde1.depth,'.','Color',red,'MarkerSize',12);
h2 = plot(datetime_utc2_adj,sonde2.depth,'.','Color',blue,'markersize',12);
hold off
legend([h1 h2 h3],'BC','ERDC','USGS')
xlim([min(datetime_utc1_adj) max(datetime_utc1_adj)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde1.deployment(1)),' - Adjusted Time'])
set(gca,'FontSize',FontSize,'LineWidth',LineWidth)
grid on

%====Save adjusted data====================================================
% Save new tables with columns for UTC & local time replaced with adjusted data
sonde1.datetime_utc = datetime_utc1_adj;
sonde1.datetime_local = datetime_local1_adj;

sonde2.datetime_utc = datetime_utc2_adj;
sonde2.datetime_local = datetime_local2_adj;

%====Save created tables in .mat files=====================================
saveFilePath = ['open-water-platform-data\',site,'\adjusted\deployments'];

option = questdlg(['Save .mat file in ',saveFilePath,'?'],'Save File','Y','N','Y');
switch option
    case 'Y'
        cd([rootpath,saveFilePath])
        save(fileName,"sonde1","sonde2","delta_t1","delta_t2")
        disp('File saved!')
    case 'N'
        disp('File not saved.')
end