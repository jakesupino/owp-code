%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotAllDeployments.m
% This script plots the entire record of data (8 plots; one for each parameter)
% from one open-water platform for both sondes (BC and ERDC).
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% 9/27/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

site = 'gull'; % CHANGE THIS

cd(['G:\My Drive\Postdoc\SMIIL\raw-data\open-water-platform-data\',site])

ds = fileDatastore(['G:\My Drive\Postdoc\SMIIL\raw-data\open-water-platform-data\',site],"ReadFcn",@load,"FileExtensions",".mat");

dat = readall(ds);

sonde1_all = table();
sonde2_all = table();

for i = 1:length(dat)
    sonde1_all = [sonde1_all;dat{i}.sonde1];
end

for j = 1:length(dat)
    if strcmp(site,'south') == 1
        skipNum = [3,5];
        if ~ismember(j,skipNum)
            sonde2_all = [sonde2_all;dat{j}.sonde2];
        end
    else
        sonde2_all = [sonde2_all;dat{j}.sonde2];
    end
end

cd(['G:\My Drive\Postdoc\SMIIL\figures\open-water-platform-figures\',site])

% Format the site name for plotting
depSite = [upper(site(1)),site(2:end)];

% Global plotting settings
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = sonde1_all.datetime_utc(end);  
NumTicks = 13;
XTick = linspace(dt1,dt2,NumTicks);
XTickFormat = "M/yy";
Legend = {'BC','ERDC'};
XLabel = 'Month/Year';
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
FontSize = 12;
LineWidth = 1;

fig1 = figure(1);
fig1.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.depth,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.depth,'Color',blue)
hold off
ax1 = gcf().CurrentAxes;
ylabel('Depth (m)')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylim([0 14])                    % Use same y limits for comparing sites
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)
title(depSite)

fig2 = figure(2);
fig2.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.p,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.p,'Color',blue)
hold off
ylabel('Pressure (psi)')
xlim([dt1 dt2])              
ylim([0 25])
ax2 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)
title(depSite)

fig3 = figure(3);
fig3.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.salinity,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.salinity,'Color',blue)
hold off
ylabel('Salinity (PSU)')
xlim([dt1 dt2])           
ylim([0 50])
ax3 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)
title(depSite)

fig4 = figure(4);
fig4.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.temperature,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.temperature,'Color',blue)
hold off
ylabel('Temperature (^oC)')
xlim([dt1 dt2])           
ylim([0 40])
ax4 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)
title(depSite)

fig5 = figure(5);
fig5.WindowState = 'maximized';
plot(sonde2_all.datetime_utc,sonde2_all.turbidity,'Color',blue)
ylabel('Turbidity (NTU)')
xlim([dt1 dt2])           
ylim([0 8500])
ax5 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend{2})
xlabel(XLabel)
title(depSite)

fig6 = figure(6);
fig6.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.pH,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.pH,'Color',blue)
hold off
ylabel('pH')
xlim([dt1 dt2])           
ylim([0 14])
ax6 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)
title(depSite)

fig7 = figure(7);
fig7.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.DO_conc,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.DO_conc,'Color',blue)
hold off
ylabel('DO concentration (\mumol/L)')
xlim([dt1 dt2])           
ylim([0 700])
ax7 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)
title(depSite)

fig8 = figure(8);
fig8.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.ORP,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.ORP,'Color',blue)
hold off
ylabel('ORP (mV)')
xlim([dt1 dt2])           
ylim([-500 500])
ax8 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)
title(depSite)

fig9 = figure(9);
fig9.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.chla,'Color',red)
ylabel('Chl a (RFU)')
xlim([dt1 dt2])           
ylim([0 250])
ax9 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend{1})
xlabel(XLabel)
title(depSite)
set(ax9,'FontSize',FontSize,'LineWidth',LineWidth,'XTick',XTick)

% Set these properties for all figures
set([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9],'FontSize',FontSize,'LineWidth',LineWidth,'XTick',XTick)

cd(['G:\My Drive\Postdoc\SMIIL\figures\open-water-platform-figures\',site])

%===Option to save plots===================================================
option = questdlg('Save as .png and .fig?','Save plots','Y','N','Y');

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

        disp(['Plots saved!'])
    case 'N'
        disp(['Plots not saved!'])
end