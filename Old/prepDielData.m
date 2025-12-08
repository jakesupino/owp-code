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
% Last updated: 4/22/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])

load([site,'-cleaned.mat'])

dat = finalQC;

cd([rootpath,'physical-data\wind-speed'])
load windSpeed.mat

cd([rootpath,'physical-data\baro-pressure'])
load baroPress.mat

%====Replace Gull met station wind speed data with ERA5 wind speeds========
% Retime cleaned Gull met data to same datetimes as sonde data
newTimes = dat.datetime_utc(1):minutes(10):dat.datetime_utc(end);
metDat_rt = retime(metDat_cleaned,newTimes,'mean');

% Retime ERA5 wind speed data to same datetimes as sonde data
era5Dat_rt = retime(era5Dat,newTimes,'previous');

metDat_rt.wspd = [];
metDat_rt = [metDat_rt,era5Dat_rt];

%====Gap fill Gull met station air T data==================================/
% Retime Baro Pressure HOBO data to same datetimes as sonde data
bpDat_rt = retime(bpDat,newTimes,'previous');
% Cut off data once original dataset ends
endDate = bpDat.datetime_utc(end);
ind_end = find(ismember(bpDat_rt.datetime_utc,endDate));
bpDat_rt(ind_end:end,{'patm' 'Tair'}) = {NaN};

% Find where there are gaps in Gull met station air T data
for i = 1:height(metDat_rt.datetime_utc)
    nid = find(isnan(metDat_rt.Tair));
    id(1:length(nid),1) = nid;
end

% Replace missing Gull T_air data with HOBO Baro Pressure T_air data
ind_nan = find(isnan(metDat_rt.Tair));
metDat_rt.Tair(ind_nan) = bpDat_rt.Tair(ind_nan);

figure(1),clf
plot(metDat_rt.datetime_utc,metDat_rt.Tair,'.','MarkerSize',4,'DisplayName','Gull Met Station')
hold on
plot(metDat_rt.datetime_utc(ind_nan),metDat_rt.Tair(ind_nan),'og','MarkerSize',6,'LineWidth',1,'DisplayName','HOBO Water-Level Logger')
ylabel('T_{air} (^oC)')
legend('show','location','best')
title('Gap-Filled Air Temperature Data')

%====Gap fill Gull met station atmos p data================================
% Replace missing Gull p_atm data with HOBO Baro Pressure p_atm data
ind_nan = find(isnan(metDat_rt.patm));
metDat_rt.patm(ind_nan) = bpDat_rt.patm(ind_nan);

figure(2),clf
plot(metDat_rt.datetime_utc,metDat_rt.patm,'.','MarkerSize',4,'DisplayName','Gull Met Station')
hold on
plot(metDat_rt.datetime_utc(ind_nan),metDat_rt.patm(ind_nan),'og','MarkerSize',6,'LineWidth',1,'DisplayName','HOBO Water-Level Logger')
ylabel('p_{atm} (hPa)')
legend('show','location','best')
title('Gap-Filled Atmospheric Pressure Data')

% Horizontally concatenate sonde and met data (already have common time vector)
dat = synchronize(dat,metDat_rt);
% dat = rmmissing(dat,'DataVariables',{'DO_conc','wspd'});

% Check how many days of data remain in synchronized data
dataPoints = daysact(dat.datetime_utc(1),dat.datetime_utc);
wholeDays = length(unique(floor(dataPoints)));

%====Export data in proper format to run Beck's R code=====================
DO_mgL = dat.DOconc / 1000 * 31.999;    % Convert umol/L to mg/L
varNames = ["DateTimeStamp","Temp","Sal","DO_obs","ATemp","BP","WSpd","Tide"];
dat_R = table(dat.datetime_utc,dat.temperature,dat.salinity,DO_mgL,dat.Tair,dat.patm,dat.wspd,dat.depth,'VariableNames',varNames);
dat_R = rmmissing(dat_R);

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