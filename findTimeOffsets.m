%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findTimeOffset.m
% This script finds whether there is a time offset between the two sondes
% (BC and ERDC) by deployment.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% 10/11/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===Read in sonde data=====================================================
clear all;close all;clc

site = 'gull'; % CHANGE THIS

ds = fileDatastore(['G:\My Drive\Postdoc\SMIIL\raw-data\open-water-platform-data\',site],"ReadFcn",@load,"FileExtensions",'.mat');

dat = readall(ds);

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

%===Plot sonde and USGS data===============================================
cd(['G:\My Drive\Postdoc\SMIIL\figures\open-water-platform-figures\',site])

red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
FontSize = 12;
LineWidth = 1;

for i = 1:length(dat)
    fig1 = figure(1);
    fig1.WindowState = 'maximized';
    h3 = plot(usgs.datetime_utc,usgs.tidal_elev,'k'); 
    hold on
    h1 = plot(dat{i}.sonde1.datetime_utc,dat{i}.sonde1.depth,'Color',red);
    
    if strcmp(site,'south') == 1
        skipNum = [3,5]; % Deployments to skip
        if ~ismember(i,skipNum)
            h2 = plot(dat{i}.sonde2.datetime_utc,dat{i}.sonde2.depth,'Color',blue);
            legend([h1 h2 h3],'BC','ERDC','USGS')
        else
            % Don't plot sonde2 data for skipped deployments
            legend([h1 h3],'BC','USGS')
        end
    else
        h2 = plot(dat{i}.sonde2.datetime_utc,dat{i}.sonde2.depth,'Color',blue);
        legend([h1 h2 h3],'BC','ERDC','USGS')
    end
    hold off
    xlim([min(dat{i}.sonde1.datetime_utc) max(dat{i}.sonde1.datetime_utc)])
    xlabel('UTC')
    ylabel('Depth (m)')
    title(['Deployment ',num2str(dat{i}.sonde1.deployment(1))])
    set(gca,'FontSize',FontSize,'LineWidth',LineWidth)
    grid on
    saveas(fig1,['dep',num2str(dat{i}.sonde1.deployment(1)),'-depth.fig'])
    saveas(fig1,['dep',num2str(dat{i}.sonde1.deployment(1)),'-depth.png'])
    pause
end