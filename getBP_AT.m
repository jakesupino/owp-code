%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getBP_AT.m
% This script combines BP (barometric pressure) and AT (air temperature) data from 
% the Gull met station and the HOBO water-level logger at TWI to produce data (.mat)
% that span the sonde measurement record.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 1/2024
% Last updated: 10/24/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
% Gull met station data
cd([rootpath,'physical-data\gull-met'])
load('gullMetStn.mat')

% HOBO water-level logger
cd([rootpath,'physical-data\HOBO']);
load('HOBO.mat')

%====Sort timetable by date================================================
hoboDat = sortrows(hoboDat);

% Convert p_atm units from [kPa] to [hPa]
hoboDat.patm = hoboDat.patm .* 10;    % [hPa]
hoboDat.Properties.VariableUnits = ["hPa","degC"];

%====Retime both datasets==================================================
% Retime cleaned Gull met data to same datetimes as sonde data 
startTime = metDat_cleaned.datetime_utc(1);
endTime = hoboDat.datetime_utc(end);
newTimes = startTime:minutes(10):endTime;

metDat_rt = retime(metDat_cleaned,newTimes,'mean'); 

hoboDat_rt = retime(hoboDat,newTimes,'mean');

%====Determine uncertainties===============================================
% Find mean absolute difference between sensors
AT_all = [metDat_rt.Tair, hoboDat_rt.Tair];
diff_AT = diff(AT_all,1,2);
two_sigma.AT = mean(abs(diff_AT),'omitmissing');    % [degC]

BP_all = [metDat_rt.patm, hoboDat_rt.patm];
diff_BP = diff(BP_all,1,2);
two_sigma.BP = mean(abs(diff_BP),'omitmissing');    % [hPa]
%%
%====Gap fill Gull met station air T data==================================
% Replace missing Gull T_air data with HOBO Baro Pressure T_air data
ind_nan = find(isnan(metDat_rt.Tair));
metDat_rt.Tair(ind_nan) = hoboDat_rt.Tair(ind_nan);

figure(1),clf
plot(metDat_rt.datetime_utc,metDat_rt.Tair,'.','MarkerSize',4)
hold on
% plot(hoboDat_rt.datetime_utc,hoboDat_rt.Tair,'og','MarkerSize',6,'LineWidth',1)
plot(metDat_rt.datetime_utc(ind_nan),metDat_rt.Tair(ind_nan),'og','MarkerSize',6,'LineWidth',1)
ylabel('Air Temperature (^oC)')
legend('Gull Met Station','HOBO Water-Level Logger Dataset','location','southeast')
title('Gap-Filled Air Temperature Data')

%====Gap fill Gull met station atmos p data================================
% Replace missing Gull p_atm data with HOBO Baro Pressure p_atm data
ind_nan = find(isnan(metDat_rt.patm));
metDat_rt.patm(ind_nan) = hoboDat_rt.patm(ind_nan);

figure(2),clf
plot(metDat_rt.datetime_utc,metDat_rt.patm,'.','MarkerSize',4)
hold on
% plot(hoboDat_rt.datetime_utc,hoboDat_rt.patm,'og','MarkerSize',6,'LineWidth',1)
plot(metDat_rt.datetime_utc(ind_nan),metDat_rt.patm(ind_nan),'og','MarkerSize',6,'LineWidth',1)
ylabel('p_{atm} (hPa)')
legend('Gull Met Station','HOBO Water-Level Logger Dataset')
title('Gap-Filled Atmospheric Pressure Data')

%====Create data table with gap-filled BP and AT data======================
BP_AT_dat = metDat_rt;
% Remove wind data
BP_AT_dat.wspd = [];
BP_AT_dat.wdir = [];

%====Option to save data===================================================
option = questdlg('Save data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'physical-data\final-dataset'])
        save('BP&AT.mat','BP_AT_dat','two_sigma')
        disp('Files saved!')
    case 'No'
        disp('Files not saved.')
end
