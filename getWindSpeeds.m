%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getWindSpeeds.m
% This script creates output tables (.mat) of wind speed data (1) measured at
% the Gull Island met station and (2) downloaded from the NOAA CDO.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% January 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

cd('G:\My Drive\Postdoc\Work\SMIIL\physical-data\wind-speed')

%====Import wind speed data measured at Gull Met Station===================
dat1 = readtable('Gull_Met_210914_220625.csv');
varNames = ["datetime_utc","Tair","wspd","wdir","patm","rhumid","source"];
varUnits = ["","degC","m/s","","hPa","",""];
dat1.Properties.VariableNames = varNames;
dat1.Properties.VariableUnits = varUnits;

dat2 = readtable('Gull_Met_WL_220625_220824.xlsx');
varNames = ["datetime_utc","wspd","wdir","Tair","rhumid","patm"];
varUnits = ["","m/s","","degC","","hPa"];
dat2.Properties.VariableNames = varNames;
dat2.Properties.VariableUnits = varUnits;
dat2.datetime_utc = datetime(dat2.datetime_utc);

dat3 = readtable('Gull_Met_WL_220824_221008.xlsx');
dat3.Properties.VariableNames = varNames;
dat3.Properties.VariableUnits = varUnits;
dat3.datetime_utc = datetime(dat3.datetime_utc);

dat4 = readtable('Gull_Met_WL_221008_to_230207.xlsx');
dat4.Properties.VariableNames = varNames;
dat4.Properties.VariableUnits = varUnits;
dat4.datetime_utc = datetime(dat4.datetime_utc);

% Concatenate data tables
dat1 = dat1(:,{'datetime_utc','wspd','wdir','Tair','rhumid','patm'});   % Restructure table 1
metDat_orig = [dat1;dat2;dat3;dat4];
metDat_orig.datetime_utc.TimeZone = 'UTC';

metDat_orig = table2timetable(metDat_orig);

clear dat1 dat2 dat3 dat4 varNames varUnits

%====QC the met data=======================================================
% Visually assess the data
figure(2),clf;
t = tiledlayout(3,1);
ax1 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.patm,'.');ylabel('p_{atm}')
ax2 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.Tair,'.');ylabel('T_{air} (^oC)')
ax3 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.wspd,'.');ylabel('Wind speed (m/s)')
title(t,'Gull Island Met Station','fontsize',16)
t.TileSpacing = 'compact';

% Link the axes
linkaxes([ax1,ax2,ax3],'x')     

% Create a table of flags for each sonde
metDat_cleaned = metDat_orig;

% Start with everything as passing (flag = 1)
flags = ones(height(metDat_orig),width(metDat_orig));
flags = array2table(flags);
varNames = ["Tair","wspd","wdir","patm","rhumid"];
flags.Properties.VariableNames = varNames;

% Gross Range Test
% INPUTS
p_low = 900;      % Lower limit (hPa)
p_high = 1100;  % Upper limit (hPa)

ind_low = find(metDat_orig.patm < p_low);
ind_high = find(metDat_orig.patm > p_high);
ind_grossRange = [ind_low;ind_high];

% Clean data
% Discard points that failed Gross Range Test
metDat_cleaned(ind_grossRange,:) = {NaN};

% Flag all failed points
flags(ind_grossRange,:) = {4};

% Highlight the flagged points
figure(3),clf;
t = tiledlayout(3,1);
ax1 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.patm,'.','MarkerSize',6);ylabel('p_{atm}')
hold on
plot(metDat_orig.datetime_utc(ind_grossRange),metDat_orig.patm(ind_grossRange),'or','MarkerSize',8)
legend('Original Data','Flagged Points')
ax2 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.Tair,'.','MarkerSize',6);ylabel('T_{air} (^oC)')
hold on
plot(metDat_orig.datetime_utc(ind_grossRange),metDat_orig.Tair(ind_grossRange),'or','MarkerSize',8)
ax3 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.wspd,'.','MarkerSize',6);ylabel('Wind speed (m/s)')
hold on
plot(metDat_orig.datetime_utc(ind_grossRange),metDat_orig.wspd(ind_grossRange),'or','MarkerSize',8)
title(t,'Gull Island Met Station - Data QC','fontsize',16)
t.TileSpacing = 'compact';

% Link the axes
linkaxes([ax1,ax2,ax3],'x')     

%====Import climatological data downloaded from NOAA=======================
% https://www.ncdc.noaa.gov/cdo-web/datasets/GHCND/stations/GHCND:USW00013724/detail
noaaDat = readtable('noaa_atlanticCity.csv');
varNames = ["station","name","latitude","longitude","elevation","date","wspd_avg","precip","Tmax","Tmin"];
varUnits = ["","","deg","deg","m","","m/s","mm","degC","degC"];
noaaDat.Properties.VariableNames = varNames;
noaaDat.Properties.VariableUnits = varUnits;
noaaDat.date.TimeZone = 'UTC';

noaaDat = table2timetable(noaaDat);

clear varNames varUnits

%====Compare Gull and NOAA daily mean wind speeds==========================
% Find daily means of cleaned wind speeds
metDat_dailyMean = retime(metDat_cleaned,'daily','mean');

figure(4),clf
plot(metDat_dailyMean.datetime_utc,metDat_dailyMean.wspd,'.-')
hold on
plot(noaaDat.date,noaaDat.wspd_avg,'.-')
ylabel('Daily Mean Wind Speed (m/s)')
legend('Gull Met Station (cleaned)','NOAA - Atlantic City Marina')

% Find date range in NOAA data that spans Gull Met data
[~,start] = ismember(metDat_dailyMean.datetime_utc(1),noaaDat.date);
[~,stop] = ismember(metDat_dailyMean.datetime_utc(end),noaaDat.date);

% Plot linear regression between cleaned Gull Met and NOAA data
tbl = table(metDat_dailyMean.wspd,noaaDat.wspd_avg(start:stop));
tbl.Properties.VariableNames = ["metDat","noaa"];
mdl = fitlm(tbl.metDat,tbl.noaa,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl.Rsquared.Ordinary,2);

figure(5),clf
h = plot(mdl,'marker','.');
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([0 ylim], [0 ylim],'--k')
legend('Data',[eqn,newline,'R^2 = ',R2],'1:1 line')
xlabel('Gull Met Station (Cleaned)')
ylabel('NOAA - Atlantic City')
title('Daily Mean Wind Speed (m/s)')
daspect([1 1 1])

%====Save the cleaned data=================================================
cd('G:\My Drive\Postdoc\Work\SMIIL\physical-data\wind-speed')
save('windSpeed.mat','metDat_cleaned','metDat_orig','noaaDat')