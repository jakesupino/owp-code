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

%===Read in sonde data=====================================================
clear all;close all;clc

site = 'gull'; % CHANGE THIS

cd(['G:\My Drive\Postdoc\SMIIL\raw-data\open-water-platform-data\',site]);

[fileName,dataPath] = uigetfile('*.mat');

load(fileName);

depNum = sonde1.deployment(1);

%===Read in USGS data======================================================
paramNames = ["agency","site_no","datetime_local","timezone","tidal_elev","qual-code"];
usgs = readtable('G:\My Drive\Postdoc\SMIIL\raw-data\usgs-data\tidal-elev.txt','TextType','string');
usgs.Properties.VariableNames = paramNames;

usgs.datetime_local.TimeZone = 'America/New_York';
datetime_utc = table(datetime(usgs.datetime_local,'TimeZone','utc'),'VariableNames',"datetime_utc");
usgs = [datetime_utc,usgs];

% Convert [ft] to [m]
tidal_elev = usgs.tidal_elev/3.281;
usgs.tidal_elev = tidal_elev;

paramUnits = ["","","","","","m",""];
usgs.Properties.VariableUnits = paramUnits;
%%
%===Plot raw sonde and USGS data===========================================
cd(['G:\My Drive\Postdoc\SMIIL\figures\open-water-platform-figures\',site])

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
title(['Deployment ',num2str(sonde1.deployment(1))])
set(gca,'FontSize',FontSize,'LineWidth',LineWidth)
grid on

%% Do this part manually, depending on the deployment

%===Remove depth offset====================================================
depth1_avg = mean(sonde1.depth);
depth2_avg = mean(sonde2.depth);
delta_depth = abs(depth2_avg - depth1_avg);

depth2_adj = sonde2.depth - delta_depth;

%===Plot adjusted data=====================================================
fig1 = figure(1);
fig1.WindowState = 'maximized';
h3 = plot(usgs.datetime_utc,usgs.tidal_elev,'k');
hold on
h1 = plot(sonde1.datetime_utc,sonde1.depth,'Color',red);
h2 = plot(sonde2.datetime_utc,depth2_adj,'Color',blue);
hold off
legend([h1 h2 h3],'BC','ERDC','USGS')
xlim([min(sonde1.datetime_utc) max(sonde1.datetime_utc)])
xlabel('UTC')
ylabel('Depth (m)')
title(['Deployment ',num2str(sonde1.deployment(1))])
set(gca,'FontSize',FontSize,'LineWidth',LineWidth)
grid on