%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getWindspeeds.m
% This script compares the wind speed data measured at the Gull Island met 
% station and downloaded from ERA5 and gives the option to save the ERA5
% data as the wind speed dataset to be use in later analyses.
%
% AUTHOR:
% Jake Supino (based on work done by Emily Chua)
% 
% DATE:
% First created: 4/2/2024
% Last Updated: 10/24/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
% Gull met station data
cd([rootpath,'physical-data\gull-met'])
load('gullMetStn.mat')

% ERA5 wind speed data
cd([rootpath,'physical-data\ERA5'])
load era5.mat

%====Compare cleaned Gull met data and ERA5 data===========================
cd([rootpath,'figures\physical-data'])

fig4 = figure(4);clf
fig4.WindowState = 'maximized';
plot(metDat_cleaned.datetime_utc,metDat_cleaned.wspd,'.k','DisplayName','Gull Met Station')
hold on
plot(era5Dat.datetime,era5Dat.wspd,'.','DisplayName','ERA5')
hold off
ylabel('Wind speed (m/s)')
legend('show')

% Compare Gull and ERA5 hourly mean wind speeds
metDat_hourlyMean = retime(metDat_cleaned,'hourly','mean');

% Find data range in ERA5 data that spans Gull Met data
[~,start] = ismember(metDat_hourlyMean.datetime_utc(1),era5Dat.datetime);
[~,stop] = ismember(metDat_hourlyMean.datetime_utc(end),era5Dat.datetime);

fig5 = figure(5);clf
fig5.WindowState = 'maximized';
plot(metDat_hourlyMean.datetime_utc,metDat_hourlyMean.wspd,'.k','DisplayName','Gull Met Station - Hourly Means')
hold on
plot(era5Dat.datetime,era5Dat.wspd,'.','DisplayName','ERA5')
hold off
ylabel('Wind speed (m/s)')
legend('show')

% Plot linear regression between cleaned Gull Met and ERA5 data
tbl = table(metDat_hourlyMean.wspd,era5Dat.wspd(start:stop));
tbl.Properties.VariableNames = ["metDat","era5"];
mdl = fitlm(tbl.metDat,tbl.era5,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl.Rsquared.Ordinary,2);

fig6 = figure(6);clf
% fig6.WindowState = 'maximized';
h = plot(mdl,'marker','.');
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([0 ylim], [0 ylim],'--k')
legend('',[eqn,newline,'R^2 = ',R2],'1:1 line','Location','southeast')
xlabel('Gull hourly wind speed (m/s)')
ylabel('ERA5 reanalysis hourly wind speed (m/s)')
% title('Hourly Wind Speed (m/s)')
title('')
daspect([1 1 1])

% Find mean absolute difference between datasets
diff_wspd = diff([metDat_hourlyMean.wspd,era5Dat.wspd(start:stop)],1,2);
two_sigma.wspd = mean(abs(diff_wspd),'omitmissing');    % [m/s]

fig7 = figure(7);clf
fig7.WindowState = 'maximized';
plot(metDat_hourlyMean.datetime_utc,metDat_hourlyMean.wspd,'.k','DisplayName','Gull Met Station')
hold on
plot(era5Dat.datetime(start:stop),era5Dat.wspd(start:stop),'.','DisplayName','ERA5')
hold off
ylabel('Hourly wind speed (m/s)')
legend('show')

%====Option to save data===================================================
option = questdlg('Save data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'physical-data\final-dataset'])
        save('windspeed.mat','era5Dat','two_sigma')
        disp('Files saved!')
    case 'No'
        disp('Files not saved.')
end