%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getGullMetData.m
% This script creates an output table (.mat) of cleaned data from the Gull Met Station
% (which measures wind speed, atmospheric pressure, and air temperature).
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 5/30/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

cd([rootpath,'physical-data\gull-met'])

%====Import the Gull Met wind speed data===================================
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

dat5 = readtable('SMIIL_Gull_Met_WL_230207_to_230524.xlsx');
dat5.Properties.VariableNames = varNames;
dat5.Properties.VariableUnits = varUnits;
dat5.datetime_utc = datetime(dat5.datetime_utc);

% Concatenate data tables
dat1 = dat1(:,{'datetime_utc','wspd','wdir','Tair','rhumid','patm'});   % Restructure table 1
metDat_orig = [dat1;dat2;dat3;dat4;dat5];
metDat_orig.datetime_utc.TimeZone = 'UTC';

metDat_orig = table2timetable(metDat_orig);

clear dat1 dat2 dat3 dat4 dat5 varNames varUnits

%====QC the Gull Met data=======================================================
% Visually assess the data
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
t = tiledlayout(3,1);
ax1 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.patm,'.k');ylabel('p_{atm}')
ax2 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.Tair,'.k');ylabel('T_{air} (^oC)')
ax3 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.wspd,'.k');ylabel('Wind speed (m/s)')
title(t,'Gull Island Met Station - Original Data','fontsize',16)
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

% Gross Range Test - Atmospheric Pressure 
% INPUTS
p_low = 900;      % Lower limit (hPa)
p_high = 1100;    % Upper limit (hPa)

ind_low = find(metDat_orig.patm < p_low);
ind_high = find(metDat_orig.patm > p_high);
ind_grossRange = [ind_low;ind_high];

% Gross Range Test - Wind Speed
% INPUTS
u_low = 0;      % Lower limit (m/s)
u_high = 25;    % Upper limit (m/s) -- highest Tropical Storm Ian wind speed was 22.7 m/s

ind_low = find(metDat_orig.wspd < u_low);
ind_high = find(metDat_orig.wspd > u_high);
ind_grossRange = [ind_grossRange;ind_low;ind_high];

% Clean data
% Discard points that failed Gross Range Test
metDat_cleaned(ind_grossRange,:) = {NaN};

% Flag all failed points
flags(ind_grossRange,:) = {4};

cd([rootpath,'figures\physical-data'])

% Highlight the flagged points
fig2 = figure(2);clf
fig2.WindowState = 'maximized';
t = tiledlayout(3,1);
ax1 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.patm,'.k','MarkerSize',6,'DisplayName','Original Data')
ylabel('p_{atm}')
hold on
plot(metDat_orig.datetime_utc(ind_grossRange),metDat_orig.patm(ind_grossRange),'or','MarkerSize',8,'DisplayName','Flagged Points')
legend('show')
ax2 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.Tair,'.k','MarkerSize',6);ylabel('T_{air} (^oC)')
hold on
plot(metDat_orig.datetime_utc(ind_grossRange),metDat_orig.Tair(ind_grossRange),'or','MarkerSize',8)
ax3 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.wspd,'.k','MarkerSize',6);ylabel('Wind speed (m/s)')
hold on
plot(metDat_orig.datetime_utc(ind_grossRange),metDat_orig.wspd(ind_grossRange),'or','MarkerSize',8)
title(t,'Gull Island Met Station - Flagged Points','fontsize',16)
t.TileSpacing = 'compact';
linkaxes([ax1,ax2,ax3],'x')     

% Plot the cleaned data
fig3 = figure(3);clf
fig3.WindowState = 'maximized';
t = tiledlayout(3,1);
ax1 = nexttile;
plot(metDat_cleaned.datetime_utc,metDat_cleaned.patm,'.k','MarkerSize',6);ylabel('p_{atm}')
hold on
plot(metDat_cleaned.datetime_utc(ind_grossRange),metDat_cleaned.patm(ind_grossRange),'or','MarkerSize',8)
ax2 = nexttile;
plot(metDat_cleaned.datetime_utc,metDat_cleaned.Tair,'.k','MarkerSize',6);ylabel('T_{air} (^oC)')
hold on
plot(metDat_cleaned.datetime_utc(ind_grossRange),metDat_cleaned.Tair(ind_grossRange),'or','MarkerSize',8)
ax3 = nexttile;
plot(metDat_cleaned.datetime_utc,metDat_cleaned.wspd,'.k','MarkerSize',6);ylabel('Wind speed (m/s)')
hold on
plot(metDat_cleaned.datetime_utc(ind_grossRange),metDat_cleaned.wspd(ind_grossRange),'or','MarkerSize',8)
title(t,'Gull Island Met Station - Cleaned Data','fontsize',16)
t.TileSpacing = 'compact';
linkaxes([ax1,ax2,ax3],'x')     

%====Option to save data===================================================
option = questdlg('Save data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'physical-data\gull-met'])
        save('gullMetStn.mat','metDat_cleaned')
        disp('Files saved!')
    case 'No'
        disp('Files not saved.')
end