%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combineFinalData.m
% This script...
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 7/22/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import final QC'd data & diel analysis results========================
site = 'gull';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
gull_params = finalQC;
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
gull_metab = diel_dtd;

site = 'north';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
north_params = finalQC;
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
north_metab = diel_dtd;

site = 'south';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
south_params = finalQC;
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
south_metab = diel_dtd;

clearvars finalQC diel_dtd diel_obs

%====Import windspeed and PAR data=========================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

cd([rootpath,'physical-data\final-dataset'])
load('par.mat')

%====Calculate daily means for each site===================================
gull_dailyAvg = retime(gull_params,'daily','mean');
north_dailyAvg = retime(north_params,'daily','mean');
south_dailyAvg = retime(south_params,'daily','mean');

wspd_dailyAvg = retime(era5Dat,'daily','mean');
par_dailyAvg = retime(parDat,'daily','mean');

cd([rootpath,'figures\stats-analyses'])

%====Calculate monthly means for each site=================================
% Gull
dt2 = dateshift(gull_metab.daystart_dt,'start','day');
gull_metab.daystart_dt = dt2;
TT_gull = synchronize(gull_dailyAvg,wspd_dailyAvg,par_dailyAvg,gull_metab);
TT_gull = removevars(TT_gull,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
TT_gull(TT_gull.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
TT_gull(TT_gull.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Get month numbers
mo = month(TT_gull.datetime_utc);  

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);   
isSpring = ismember(mo,[3 4 5]);    
isSummer = ismember(mo,[6 7 8]);    
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

% Add a column for season number
TT_gull.season(indWinter) = 1;
TT_gull.season(indSpring) = 2;
TT_gull.season(indSummer) = 3;
TT_gull.season(indFall) = 4;
TT_gull.season = categorical(TT_gull.season);

% Add a column for site number
TT_gull.site = repelem(2,height(TT_gull))';
TT_gull.site = categorical(TT_gull.site);

% North
dt2 = dateshift(north_metab.daystart_dt,'start','day');
north_metab.daystart_dt = dt2;
TT_north = synchronize(north_dailyAvg,wspd_dailyAvg,par_dailyAvg,north_metab);
TT_north = removevars(TT_north,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
TT_north(TT_north.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
TT_north(TT_north.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Get month numbers
mo = month(TT_north.datetime_utc);  

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);   
isSpring = ismember(mo,[3 4 5]);    
isSummer = ismember(mo,[6 7 8]);    
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

% Add a column for season number
TT_north.season(indWinter) = 1;
TT_north.season(indSpring) = 2;
TT_north.season(indSummer) = 3;
TT_north.season(indFall) = 4;
TT_north.season = categorical(TT_north.season);

% Add a column for site number
TT_north.site = repelem(1,height(TT_north))';
TT_north.site = categorical(TT_north.site);

% South
dt2 = dateshift(south_metab.daystart_dt,'start','day');
south_metab.daystart_dt = dt2;
TT_south = synchronize(south_dailyAvg,wspd_dailyAvg,par_dailyAvg,south_metab);
TT_south = removevars(TT_south,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
TT_south(TT_south.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
TT_south(TT_south.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Get month numbers
mo = month(TT_south.datetime_utc);  

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);   
isSpring = ismember(mo,[3 4 5]);    
isSummer = ismember(mo,[6 7 8]);    
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

% Add a column for season number
TT_south.season(indWinter) = 1;
TT_south.season(indSpring) = 2;
TT_south.season(indSummer) = 3;
TT_south.season(indFall) = 4;
TT_south.season = categorical(TT_south.season);

% Add a column for site number
TT_south.site = repelem(3,height(TT_south))';
TT_south.site = categorical(TT_south.site);

%====Create a table to store predictor variables for all sites=============
tbl_gull = timetable2table(TT_gull);
tbl_gull.datetime_utc = [];

tbl_north = timetable2table(TT_north);
tbl_north.datetime_utc = [];

tbl_south = timetable2table(TT_south);
tbl_south.datetime_utc = [];

data = [tbl_north;tbl_gull;tbl_south];