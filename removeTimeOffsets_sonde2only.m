%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% removeTimeOffsets_sonde2only.m
% This script removes time offsets from UTC for deployments that only have sonde 2 data. 
% Refer to the Collaborative Lab Notebook, which documents which deployments have time offsets.
%
% Run this script for deployments that only have data for Sonde 2 (ERDC)!!
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 10/30/2023
% Last amended: 3/13/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

fig = uifigure;
site = uiconfirm(fig,"Select the platform","Site selection","Options",["gull","north","south"]);
close(fig)

%===Read in sonde data for a specific deployment===========================
cd([rootpath,'open-water-platform-data\',site,'\original\deployments'])

[fileName,dataPath] = uigetfile('*.mat');

load(fileName);

depNum = sonde2.deployment(1);

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
h2 = plot(sonde2.datetime_utc,sonde2.depth,'Color',blue);
legend([h2 h3],'ERDC','USGS')
hold off
xlim([min(sonde2.datetime_utc) max(sonde2.datetime_utc)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde2.deployment(1)),' - Original'])
set(gca,'FontSize',FontSize,'LineWidth',LineWidth)
grid on

%===Fix time offset(s) from UTC============================================
% Find peaks in USGS data
dups = find(diff(usgs.datetime_utc)<minutes(6));    % Find duplicate USGS timepoints
usgs(dups,:) = [];                                  % Remove duplicate USGS timepoints
smoothed0 = smoothdata(usgs.tidal_elev);
[pks0,locs0] = findpeaks(smoothed0,'MinPeakDistance',110);

% Find peaks in sonde data; smooth first to remove local peaks
smoothed2 = smoothdata(sonde2.depth);
% Change this value depending on sampling interval (6, 10, or 12 min)
minPkDist = 55;   % 12 min
% minPkDist = 65;    % 10 min
% minPkDist = 115;   % 6 min
[pks2,locs2] = findpeaks(smoothed2,'MinPeakDistance',minPkDist);

% Plot peaks in data to double-check
fig2= figure(2);
fig2.WindowState = 'maximized';
h0 = plot(usgs.datetime_utc,usgs.tidal_elev,'k');
hold on
plot(usgs.datetime_utc(locs0),usgs.tidal_elev(locs0),'ok')
h2 = plot(sonde2.datetime_utc,sonde2.depth,'Color',blue);
plot(sonde2.datetime_utc(locs2),sonde2.depth(locs2),'o','Color',blue)
hold off
legend([h0 h2],'USGS','ERDC')
xlim([min(sonde2.datetime_utc) max(sonde2.datetime_utc)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde2.deployment(1)),' - Identified Peaks'])

disp('Press enter to adjust times')
pause

% Find time offset between third sonde peak and corresponding USGS peak
diff2 = usgs.datetime_utc(locs0) - sonde2.datetime_utc(locs2(3));
ind2 = find(diff2==max(diff2(diff2<0)));
delta_t2 = diff2(ind2);
if abs(delta_t2) > hours(12)
    delta_t2 = diff2(ind2+1);
else
    % Do nothing
end

% Adjust UTC time
datetime_utc2_adj = sonde2.datetime_utc + delta_t2;

% Adjust local time
datetime_local2_adj = sonde2.datetime_local + delta_t2;

%===Plot adjusted data=============================================
fig3 = figure(3);
fig3.WindowState = 'maximized';
h3 = plot(usgs.datetime_utc,usgs.tidal_elev,'k');
hold on
h2 = plot(datetime_utc2_adj,sonde2.depth,'Color',blue);
hold off
legend([h2 h3],'ERDC','USGS')
xlim([min(datetime_utc2_adj) max(datetime_utc2_adj)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde2.deployment(1)),' - Adjusted Time'])
set(gca,'FontSize',FontSize,'LineWidth',LineWidth)

%====Save adjusted data====================================================
% Save new tables with columns for UTC & local time replaced with adjusted data
sonde2.datetime_utc = datetime_utc2_adj;
sonde2.datetime_local = datetime_local2_adj;

%====Save created tables in .mat files=====================================
saveFilePath = ['open-water-platform-data\',site,'\adjusted\deployments'];

option = questdlg(['Save .mat file in ',saveFilePath,'?'],'Save File','Y','N','Y');
switch option
    case 'Y'
        cd([rootpath,saveFilePath])
        save(fileName,"sonde2","delta_t2")
        disp('File saved!')
    case 'N'
        disp('File not saved.')
end