%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge&PlotDeployments.m
% This script merges and plots the entire record of data (8 plots; one for 
% each parameter) from one open-water platform for both sondes (BC and ERDC).
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 9/14/2023
% Last amended: 1/24/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

fig = uifigure;
site = uiconfirm(fig,"Select the platform","Site selection","Options",["gull","north","south"]);
close(fig)

fig = uifigure;
dataset = uiconfirm(fig,"Plot original or adjusted data?","Dataset","Options",["original","adjusted"],"DefaultOption","adjusted");
close(fig)

switch dataset
    case "original"
        ds = fileDatastore([rootpath,'open-water-platform-data\',site,'\original\deployments'],"ReadFcn",@load,"FileExtensions",".mat");
    case "adjusted"
        ds = fileDatastore([rootpath,'open-water-platform-data\',site,'\adjusted\deployments'],"ReadFcn",@load,"FileExtensions",".mat");
end

dat = readall(ds);

sonde1_all = table();
sonde2_all = table();

for i = 1:length(dat)
    if strcmp(site,'south') == 1
        skipNum = 14;   % Skip sonde2 data for South Deployment 15
        if ~ismember(i,skipNum)
            sonde1_all = [sonde1_all;dat{i}.sonde1];
        end
    else
        sonde1_all = [sonde1_all;dat{i}.sonde1];
    end
end

for j = 1:length(dat)
    if strcmp(site,'south') == 1
        skipNum = [3,5];    % Skip sonde1 data for South Deployments 4 & 6
        if ~ismember(j,skipNum)
            sonde2_all = [sonde2_all;dat{j}.sonde2];
        end
    else
        sonde2_all = [sonde2_all;dat{j}.sonde2];
    end
end

% Format the site name for plotting
depSite = [upper(site(1)),site(2:end)];

% Global plotting settings
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = max(sonde1_all.datetime_utc(end), sonde2_all.datetime_utc(end));
Legend = {'BC','ERDC'};
XLabel = 'UTC';
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
FontSize = 14;
LineWidth = 1;

% Find indices of deployment changes
switch site
    case 'gull'
        ind_dep = find(diff(sonde1_all.deployment) > 0);
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
            'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
            'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
            'Deployment 15','Deployment 16'};
    case 'north'
        ind_dep = find(diff(sonde1_all.deployment) > 0);
        label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8',...
            'Deployment 9','Deployment 10','Deployment 11','Deployment 12',...
            'Deployment 13','Deployment 14','Deployment 15','Deployment 16'};
    case 'south'
        ind_dep15 = find(ismember(sonde2_all.deployment,15),1);   % BC sonde didn't have a Dep 15; find where it is in erdc data
        ind_dep15 = interp1(unique(sonde1_all.datetime_utc),1:length(unique(sonde1_all.datetime_utc)),sonde2_all.datetime_utc(ind_dep15),'nearest');
        ind_dep = find(diff(sonde1_all.deployment) > 0);
        ind_dep = sort([ind_dep;ind_dep15]);
        label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
            'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
            'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
            'Deployment 14','Deployment 15','Deployment 16'};
end

fig1 = figure(1);
fig1.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.depth,'.','Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.depth,'.','Color',blue)
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ax1 = gcf().CurrentAxes;
ylabel('Depth (m)')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylim([-0.5 4])                    % Use same y limits for comparing sites
legend(Legend)
xlabel(XLabel)
title(depSite)

fig2 = figure(2);
fig2.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.p,'.','Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.p,'.','Color',blue)
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ylabel('Pressure (psi)')
xlim([dt1 dt2])    
% ylim([0 30])
ylim([-1 8])
ax2 = gcf().CurrentAxes;
legend(Legend)
xlabel(XLabel)
title(depSite)

fig3 = figure(3);
fig3.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.salinity,'.','Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.salinity,'.','Color',blue)
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ylabel('Salinity (PSU)')
xlim([dt1 dt2])           
ylim([0 50])
ax3 = gcf().CurrentAxes;
legend(Legend)
xlabel(XLabel)
title(depSite)

fig4 = figure(4);
fig4.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.temperature,'.','Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.temperature,'.','Color',blue)
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ylabel('Temperature (^oC)')
xlim([dt1 dt2])           
ylim([-10 40])
ax4 = gcf().CurrentAxes;
legend(Legend)
xlabel(XLabel)
title(depSite)

fig5 = figure(5);
fig5.WindowState = 'maximized';
plot(sonde2_all.datetime_utc,sonde2_all.turbidity,'.','Color',blue)
hold on
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ylabel('Turbidity (NTU)')
xlim([dt1 dt2])           
ylim([0 8500])
ax5 = gcf().CurrentAxes;
legend(Legend{2})
xlabel(XLabel)
title(depSite)

fig6 = figure(6);
fig6.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.pH,'.','Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.pH,'.','Color',blue)
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ylabel('pH')
xlim([dt1 dt2])           
ylim([0 14])
ax6 = gcf().CurrentAxes;
legend(Legend)
xlabel(XLabel)
title(depSite)

fig7 = figure(7);
fig7.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.DO_conc,'.','Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.DO_conc,'.','Color',blue)
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ylabel('DO concentration (\mumol/L)')
xlim([dt1 dt2])           
ylim([0 700])
ax7 = gcf().CurrentAxes;
legend(Legend)
xlabel(XLabel)
title(depSite)

fig8 = figure(8);
fig8.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.ORP,'.','Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.ORP,'.','Color',blue)
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ylabel('ORP (mV)')
xlim([dt1 dt2])           
ylim([-500 600])
ax8 = gcf().CurrentAxes;
legend(Legend)
xlabel(XLabel)
title(depSite)

fig9 = figure(9);
fig9.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.chla,'.','Color',red)
hold on
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ylabel('Chl a (RFU)')
xlim([dt1 dt2])           
ylim([0 410])
ax9 = gcf().CurrentAxes;
legend(Legend{1})
xlabel(XLabel)
title(depSite)

% Set same properties for all figures
% set([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9],'FontSize',FontSize,'LineWidth',LineWidth,'XTick',XTick)
set([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9],'FontSize',FontSize,'LineWidth',LineWidth)


%====Save created tables in .mat files=====================================
switch dataset
    case "original"
        saveFilePath = ['open-water-platform-data\',site,'\original\merged'];
    case "adjusted"
        saveFilePath = ['open-water-platform-data\',site,'\adjusted\merged'];
end
option = questdlg(['Save .mat file in ',saveFilePath,'?'],'Save File','Y','N','Y');
cd([rootpath,saveFilePath])

switch option
    case 'Y'
        switch dataset
            case "original"
                save(['alldeps-',site,'-raw.mat'],"sonde1_all","sonde2_all")
            case "adjusted"
                save(['alldeps-',site,'-adj.mat'],"sonde1_all","sonde2_all")
        end
        disp('File saved!')
    case 'N'
        disp('File not saved!')
end

%===Option to save plots===================================================
switch dataset
    case "original"
        saveFilePath = ['figures\open-water-platform-figures\',site,'\original\merged'];
    case "adjusted"
        saveFilePath = ['figures\open-water-platform-figures\',site,'\adjusted\merged'];
end

option = questdlg(['Save plots as .png and .fig in ',saveFilePath,'?'],'Save plots','Y','N','Y');
cd([rootpath,saveFilePath])

switch option
    case 'Y'
        saveas(fig1,'alldeps-depth.png')
        saveas(fig1,'alldeps-depth.fig')

        saveas(fig2,'alldeps-pressure.png')
        saveas(fig2,'alldeps-pressure.fig')

        saveas(fig3,'alldeps-salinity.png')
        saveas(fig3,'alldeps-salinity.fig')

        saveas(fig4,'alldeps-temperature.png')
        saveas(fig4,'alldeps-temperature.fig')

        saveas(fig5,'alldeps-turbidity.png')
        saveas(fig5,'alldeps-turbidity.fig')

        saveas(fig6,'alldeps-pH.png')
        saveas(fig6,'alldeps-pH.fig')

        saveas(fig7,'alldeps-do.png')
        saveas(fig7,'alldeps-do.fig')

        saveas(fig8,'alldeps-orp.png')
        saveas(fig8,'alldeps-orp.fig')

        saveas(fig9,'alldeps-chla.png')
        saveas(fig9,'alldeps-chla.fig')

        disp('Plots saved!')
    case 'N'
        disp('Plots not saved!')
end