%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge&PlotDeployments.m
% This script merges and the entire record of sonde data from one
% open-water platform. The option is given to merge and plot the raw "original" 
% or the time-offset "adjusted" data. If the "adjusted" data are chosen, depth 
% offsets for individual deployments are further adjusted by pinning the
% each deployment's mean depth to a global mean value.
% Merged results are then plotted for both sondes (BC and ERDC).
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 9/14/2023
% Last amended: 6/26/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

fig = uifigure;
site = uiconfirm(fig,"Select the platform","Site selection","Options",["gull","north","south"]);
close(fig)

fig = uifigure;
dataset = uiconfirm(fig,"Load original or time-adjusted data?","Dataset","Options",["original","adjusted"],"DefaultOption","adjusted");
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
location = 'northwest';
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
            'Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
    case 'north'
        ind_dep = find(diff(sonde1_all.deployment) > 0);
        label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8',...
            'Deployment 9','Deployment 10','Deployment 11','Deployment 12',...
            'Deployment 13','Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
    case 'south'
        ind_dep15 = find(ismember(sonde2_all.deployment,15),1);   % BC sonde didn't have a Dep 15; find where it is in ERDC data instead
        ind_dep15 = interp1(unique(sonde1_all.datetime_utc),1:length(unique(sonde1_all.datetime_utc)),sonde2_all.datetime_utc(ind_dep15),'nearest');
        ind_dep = find(diff(sonde1_all.deployment) > 0);
        ind_dep = sort([ind_dep;ind_dep15]);
        label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
            'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
            'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
            'Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
end

%===Remove depth offset from time-adjusted data============================
if strcmp(dataset,'adjusted')
    % Find mean depth for each deployment, excluding zero/negative values (data dropouts; mostly a problem in Deployment 1)
    depth1 = table(sonde1_all.deployment,sonde1_all.depth,'VariableNames',["deployment","depth"]);  
    ind_rm = find(depth1.depth <= 0);
    depth1.depth(ind_rm) = NaN;
    depth2 = table(sonde2_all.deployment,sonde2_all.depth,'VariableNames',["deployment","depth"]);
    ind_rm = find(depth2.depth <= 0);
    depth2.depth(ind_rm) = NaN;

    means1 = grpstats(depth1,"deployment","mean");
    means2 = grpstats(depth2,"deployment","mean");
    
    % For each sonde, find the mean depth across all deployments that aren't erroneously high, excluding zero values
    switch site
        case 'gull'
            alldepths1_avg = mean(nonzeros(means1.mean_depth([1 2 3:8 10 11 13 14 15 16]))); % BC: Dep# 1, 2, 5-10, 12, 13, 15, 16, 17, 18
            alldepths2_avg = mean(nonzeros(means2.mean_depth([1:6 8:10 13 14 15 16])));      % ERDC: Dep# 1, 2, 5-8, 10-12, 15, 16, 17, 18
        case 'north'
            alldepths1_avg = mean(nonzeros(means1.mean_depth([1:6 8 10 11 12 13 14])));   % BC: Dep# 2, 6-10, 12, 14, 15, 16, 17, 18
            alldepths2_avg = mean(nonzeros(means2.mean_depth([1:4 7 8 10 11 12 13 14]))); % ERDC: Dep# 2, 6-8, 11, 12, 14, 15, 16, 17, 18
        case 'south'
            alldepths1_avg = mean(nonzeros(means1.mean_depth([1:7 9:11 14 15 16])));    % BC: Dep# 1, 2, 4-8, 10-12, 16, 17, 18
            alldepths2_avg = mean(nonzeros(means2.mean_depth([1:6 9 10 12:14 15])));    % ERDC: Dep# 1, 2, 5, 7-9, 12, 13, 15, 16, 17, 18
    end

    % Sonde 1
    delta_depth1 = means1.mean_depth - alldepths1_avg;
    delta_depth1_rep = [];
    for i = 1:height(means1)
        % j = means1.deployment(i);
        dd = repelem(delta_depth1(i),means1.GroupCount(i))';  % Replicate depth differences for length of given deployment
        delta_depth1_rep = [delta_depth1_rep;dd];
    end
    depth1_adj = sonde1_all.depth - delta_depth1_rep;

    delta_depth2 = means2.mean_depth - alldepths2_avg;
    delta_depth2_rep = [];
    for i = 1:height(means2)
        dd = repelem(delta_depth2(i),means2.GroupCount(i))';  % Replicate depth differences for length of given deployment
        delta_depth2_rep = [delta_depth2_rep;dd];
    end
    depth2_adj = sonde2_all.depth - delta_depth2_rep;

    %===Adjust pressure data using adjusted depth data=========================
    p1_adj = sonde1_all.density.*1000*9.81.*depth1_adj/6894.76;
    p2_adj = sonde2_all.density.*1000*9.81.*depth2_adj/6894.76;

    %====Save adjusted data====================================================
    % Save new tables with columns for depth and p replaced with adjusted data
    sonde1_all.depth = depth1_adj;
    sonde1_all.p = p1_adj;
    sonde2_all.depth = depth2_adj;
    sonde2_all.p = p2_adj;
end

%===Plot the results=======================================================
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(sonde1_all.datetime_utc,sonde1_all.depth,'.','Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.depth,'.','Color',blue)
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ax1 = gcf().CurrentAxes;
ylabel('Depth (m)')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
if strcmp(dataset,'adjusted')
    ylim([-0.5 4])              % Use same y limits for comparing sites
end
legend(Legend,'Location',location)
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
if strcmp(dataset,'adjusted')
    ylim([0 30])
end
ylim([-1 8])
ax2 = gcf().CurrentAxes;
legend(Legend,'Location',location)
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
if strcmp(dataset,'adjusted')
    ylim([0 50])
end
ax3 = gcf().CurrentAxes;
legend(Legend,'Location',location)
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
if strcmp(dataset,'adjusted')
    ylim([-10 40])
end
ax4 = gcf().CurrentAxes;
legend(Legend,'Location',location)
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
if strcmp(dataset,'adjusted')
    ylim([0 8500])
end
ax5 = gcf().CurrentAxes;
legend(Legend{2},'Location',location)
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
if strcmp(dataset,'adjusted')
    ylim([0 14])
end
ax6 = gcf().CurrentAxes;
legend(Legend,'Location',location)
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
if strcmp(dataset,'adjusted')
    ylim([0 700])
end
ax7 = gcf().CurrentAxes;
legend(Legend,'Location',location)
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
if strcmp(dataset,'adjusted')
    ylim([-500 600])
end
ax8 = gcf().CurrentAxes;
legend(Legend,'Location',location)
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
if strcmp(dataset,'adjusted')
    ylim([0 410])
end
ax9 = gcf().CurrentAxes;
legend(Legend{1},'Location',location)
xlabel(XLabel)
title(depSite)

% Set same properties for all figures
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
        saveFilePath = ['figures\open-water-platform\',site,'\original\merged'];
    case "adjusted"
        saveFilePath = ['figures\open-water-platform\',site,'\adjusted\merged'];
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