%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepDielData.m
% This script preps the sonde data in the proper format to run in Beck
% et al. (2015)'s R code
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 2/8/2024
% Last updated: 9/24/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])

load('finalQC.mat')

dat = finalQC;

cd([rootpath,'physical-data\final-dataset'])
load windspeed.mat
load BP&AT.mat

% Retime physical data to same datetimes as sonde data
newTimes = dat.datetime_utc(1):minutes(10):dat.datetime_utc(end);
era5dat_rt = retime(era5Dat,newTimes,'previous');
BP_AT_dat_rt = retime(BP_AT_dat,newTimes,'previous');

% Horizontally concatenate sonde and physical data
dat = synchronize(dat,BP_AT_dat_rt,era5dat_rt);

% Check how many days of data remain in synchronized data
dataPoints = daysact(dat.datetime_utc(1),dat.datetime_utc);
wholeDays = length(unique(floor(dataPoints)));

%====Export data in proper format to run Beck's R code=====================
DO_mgL = dat.DOconc / 1000 * 31.999;    % Convert umol/L to mg/L
varNames = ["DateTimeStamp","Temp","Sal","DO_obs","ATemp","BP","WSpd","Tide"];
dat_R = table(dat.datetime_utc,dat.temperature,dat.salinity,DO_mgL,dat.Tair,dat.patm,dat.wspd,dat.depth,'VariableNames',varNames);
dat_R = rmmissing(dat_R);

% Mean water column depth = d + D [m]
% See Collab Lab Notebook - Table 2 for manual measurements for D for each site
switch site
    case 'Gull'
        H = mean(dat.depth,'omitnan') + 0.42;
    case 'North'
        H = mean(dat.depth,'omitnan') + 0.47;
    case 'South'
        H = mean(dat.depth,'omitnan') + 0.80;
end
disp(sprintf(['H = ',num2str(H),' m \n']))

%====Save the synchronized data============================================
option = questdlg('Save data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'\diel-method\owp-data\final-qc'])
        writetable(dat_R,[site,'_obs.csv'])
        save([site,'_obs.mat'],'dat')
        disp('Files saved!')
    case 'No'
        disp('Files not saved.')
end