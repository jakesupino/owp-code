%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotOneDeployment.m
% This script plots the data from one open-water platform for one deployment
% from both sondes (BC and ERDC).
%
% AUTHOR:
% Emily Chua
%
% DATE:
% 9/27/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

fig = uifigure;
site = uiconfirm(fig,"Select the platform","Site selection","Options",["gull","north","south"]);
close(fig)

cd(['G:\My Drive\Postdoc\SMIIL\open-water-platform-data\raw-data\',site])

[fileName,dataPath] = uigetfile('*.mat');

load(fileName);

% Extract the deployment number from the filename and make a column to add to data table
depNum = extractBetween(fileName,"dep","-");
depNum = str2double(depNum);

% Format the site name for plotting
depSite = [upper(site(1)),site(2:end)];

% Plot the data
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde

figure,clf
tl = tiledlayout(4,2,'TileSpacing','Tight');
title(tl,[depSite,' Deployment #',num2str(depNum),' - Both Sondes'])
xlabel(tl,'Time (UTC)')

switch depNum
    % Plots if both sondes telemetered (Deployments 1-2)
    case{1, 2}
        nexttile(1)
        plot(sonde1.datetime_utc,sonde1.depth,'color',red)
        hold on
        plot(sonde2.datetime_utc,sonde2.depth,'color',blue)
        hold off
        title('Depth')
        ylabel('Depth (m)')
        xlim('tight')

        % Add a global legend to label sondes
        lg = legend('BC','ERDC');
        lg.Layout.Tile = 'north';

        nexttile(2)
        plot(sonde2.datetime_utc,sonde2.pH,'color',blue)
        title('pH')
        ylabel('pH')
        xlim('tight')

        nexttile(3)
        plot(sonde1.datetime_utc,sonde1.temperature,'color',red)
        hold on
        plot(sonde2.datetime_utc,sonde2.temperature,'color',blue)
        title('Temperature')
        ylabel('Temperature (^oC)')
        xlim('tight')

        nexttile(4)
        plot(sonde1.datetime_utc,sonde1.DO_conc,'color',red)
        hold on
        plot(sonde2.datetime_utc,sonde2.DO_conc,'color',blue)
        hold off
        title('DO concentration')
        ylabel('DO (\mumol/L)')
        xlim('tight')

        nexttile(5)
        plot(sonde1.datetime_utc,sonde1.salinity,'color',red)
        hold on
        plot(sonde2.datetime_utc,sonde2.salinity,'color',blue)
        hold off
        title('Salinity')
        ylabel('Salinity (PSU)')
        xlim('tight')

        nexttile(7)
        plot(sonde2.datetime_utc,sonde2.turbidity,'Color',blue)
        title('Turbidity')
        ylabel('Turbidity (NTU)')
        xlim('tight')

        nexttile(8)
        plot(sonde1.datetime_utc,sonde1.chla,'color',red)
        title('Chl a')
        ylabel('Chl a (RFU)')
        xlim('tight')

    % Plots if one or more sondes were internally logged (Deployment 4 onwards)
    case{4,5,6,7,8,9,10,11,12,13,14,15}
        if strcmp(site,'south') == 1
            if depNum == 4 || depNum == 6
                nexttile(1)
                plot(sonde1.datetime_utc,sonde1.depth,'color',red)
                ax1 = gcf().CurrentAxes;
                title('Depth')
                ylabel('Depth (m)')
                xlim('tight');ylim('tight')

                % Add a global legend to label sondes
                lg = legend('BC');
                lg.Layout.Tile = 'north';

                if depNum == 6
                    nexttile(2)
                    plot(sonde1.datetime_utc,sonde1.pH,'color',red)
                    title('pH')
                    ylabel('pH')
                    xlim(ax1.XLim);ylim('tight')
                end

                nexttile(3)
                plot(sonde1.datetime_utc,sonde1.temperature,'color',red)
                title('Temperature')
                ylabel('Temperature (^oC)')
                xlim(ax1.XLim);ylim('tight')

                nexttile(4)
                plot(sonde1.datetime_utc,sonde1.DO_conc,'color',red)
                title('DO concentration')
                ylabel('DO (\mumol/L)')
                xlim(ax1.XLim);ylim('tight')

                nexttile(5)
                plot(sonde1.datetime_utc,sonde1.salinity,'color',red)
                title('Salinity')
                ylabel('Salinity (PSU)')
                xlim(ax1.XLim);ylim('tight')

                nexttile(8)
                plot(sonde1.datetime_utc,sonde1.chla,'color',red)
                title('Chl a')
                ylabel('Chl a (RFU)')
                xlim(ax1.XLim);ylim('tight')

            elseif depNum == 15
                nexttile(1)
                plot(sonde2.datetime_utc,sonde2.depth,'color',blue)
                ax1 = gcf().CurrentAxes;
                title('Depth')
                ylabel('Depth (m)')
                xlim('tight');ylim('tight')

                % Add a global legend to label sondes
                lg = legend('ERDC');
                lg.Layout.Tile = 'north';

                nexttile(2)
                plot(sonde2.datetime_utc,sonde2.pH,'color',blue)
                title('pH')
                ylabel('pH')
                xlim(ax1.XLim);ylim('tight')

                nexttile(3)
                plot(sonde2.datetime_utc,sonde2.temperature,'color',blue)
                title('Temperature')
                ylabel('Temperature (^oC)')
                xlim(ax1.XLim);ylim('tight')

                nexttile(4)
                plot(sonde2.datetime_utc,sonde2.DO_conc,'color',blue)
                title('DO concentration')
                ylabel('DO (\mumol/L)')
                xlim(ax1.XLim);ylim('tight')

                nexttile(5)
                plot(sonde2.datetime_utc,sonde2.salinity,'color',blue)
                title('Salinity')
                ylabel('Salinity (PSU)')
                xlim(ax1.XLim);ylim('tight')

                nexttile(6)
                plot(sonde2.datetime_utc,sonde2.ORP,'color',blue)
                title('ORP')
                ylabel('ORP (mV)')
                xlim(ax1.XLim);ylim('tight')

                nexttile(8)
                plot(sonde2.datetime_utc,sonde2.turbidity,'color',blue)
                title('Turbidity')
                ylabel('Turbidity (NTU)')
                xlim(ax1.XLim);ylim('tight')
            end

        else
            nexttile(1)
            plot(sonde1.datetime_utc,sonde1.depth,'color',red)
            hold on
            plot(sonde2.datetime_utc,sonde2.depth,'color',blue)
            hold off
            ax1 = gcf().CurrentAxes;
            title('Depth')
            ylabel('Depth (m)')
            xlim('tight');ylim('tight')

            % Add a global legend to label sondes
            lg = legend('BC','ERDC');
            lg.Layout.Tile = 'north';

            nexttile(2)
            plot(sonde1.datetime_utc,sonde1.pH,'color',red)
            hold on
            plot(sonde2.datetime_utc,sonde2.pH,'color',blue)
            title('pH')
            ylabel('pH')
            xlim(ax1.XLim);ylim('tight')

            nexttile(3)
            plot(sonde1.datetime_utc,sonde1.temperature,'color',red)
            hold on
            plot(sonde2.datetime_utc,sonde2.temperature,'color',blue)
            title('Temperature')
            ylabel('Temperature (^oC)')
            xlim(ax1.XLim);ylim('tight')

            nexttile(4)
            plot(sonde1.datetime_utc,sonde1.DO_conc,'color',red)
            hold on
            plot(sonde2.datetime_utc,sonde2.DO_conc,'color',blue)
            hold off
            title('DO concentration')
            ylabel('DO (\mumol/L)')
            xlim(ax1.XLim);ylim('tight')

            nexttile(5)
            plot(sonde1.datetime_utc,sonde1.salinity,'color',red)
            hold on
            plot(sonde2.datetime_utc,sonde2.salinity,'color',blue)
            hold off
            title('Salinity')
            ylabel('Salinity (PSU)')
            xlim(ax1.XLim);ylim('tight')

            nexttile(6)
            plot(sonde1.datetime_utc,sonde1.ORP,'color',red)
            hold on
            plot(sonde2.datetime_utc,sonde2.ORP,'color',blue)
            hold off
            title('ORP')
            ylabel('ORP (mV)')
            xlim(ax1.XLim);ylim('tight')

            nexttile(7)
            plot(sonde2.datetime_utc,sonde2.turbidity,'Color',blue)
            title('Turbidity')
            ylabel('Turbidity (NTU)')
            xlim(ax1.XLim);ylim('tight')

            nexttile(8)
            plot(sonde1.datetime_utc,sonde1.chla,'color',red)
            title('Chl a')
            ylabel('Chl a (RFU)')
            xlim(ax1.XLim);ylim('tight')
        end
end

cd(['G:\My Drive\Postdoc\SMIIL\figures\open-water-platform-figures\',site])