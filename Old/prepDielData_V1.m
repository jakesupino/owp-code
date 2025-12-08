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
% Last updated: 4/12/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
prompt = {'Choose the platform'};
answer = questdlg(prompt,'Platform Selection','Gull','North','South','Cancel');
switch answer
    case 'Gull'
        cd([rootpath,'open-water-platform-data\gull\cleaned'])
        site = 'Gull';
    case 'North'
        cd([rootpath,'open-water-platform-data\north\cleaned'])
        site = 'North';
    case 'South'
        cd([rootpath,'open-water-platform-data\south\cleaned'])
        site = 'South';
end   

prompt = {'Choose the sonde dataset'};
answer = questdlg(prompt,'Sonde Selection','Sonde 1','Sonde 2','Cancel','Cancel');
switch answer
    case 'Sonde 1'
        load([site,'-bc-cleaned.mat'])
        dat = sonde1_cleaned;
        sondename = 'BC';
    case 'Sonde 2'
        load([site,'-erdc-cleaned.mat'])
        dat = sonde2_cleaned;
        sondename = 'ERDC';
end

cd([rootpath,'physical-data\wind-speed'])
load windSpeed.mat

cd([rootpath,'physical-data\par'])
load par.mat

cd([rootpath,'physical-data\baro-pressure'])
load baroPress.mat

% Retime cleaned Gull met data to same datetimes as sonde data 
newTimes = dat.datetime_utc(1):minutes(10):dat.datetime_utc(end);
metDat_rt = retime(metDat_cleaned,newTimes,'mean'); 

%====Gap fill Gull met station wind speed data using NOAA daily means======
% Retime NOAA data to same datetimes as sonde data 
noaaDat_rt = retime(noaaDat,newTimes,'previous');

% Find where there are gaps in Gull met station wind speed data
for i = 1:height(metDat_rt.datetime_utc)
    nid = find(isnan(metDat_rt.wspd));
    id(1:length(nid),1) = nid;
end

% figure(1),clf
% plot(metDat_rt.datetime_utc,metDat_rt.wspd,'.-')
% hold on
% plot(noaaDat_rt.date(id),noaaDat_rt.wspd_avg(id),'or','MarkerSize',4)
% ylabel('Wind Speed (m/s)')
% legend('Gull Met Station','NOAA Daily Mean')
% title('Gap-Filled Wind Speed Data')

% Replace missing Gull windspeed data with NOAA windspeed data
ind_nan = find(isnan(metDat_rt.wspd));
metDat_rt.wspd(ind_nan) = noaaDat_rt.wspd_avg(ind_nan);

%====Gap filling Gull met station air T data===============================
% Retime PAR data to same datetimes as sonde data 
parDat_rt = retime(parDat,newTimes,'previous');
% Cut off data once original dataset ends
endDate = parDat.datetime_utc(end);
ind_end = find(ismember(parDat_rt.datetime_utc,endDate));
parDat_rt(ind_end:end,{'Tair' 'light_lux' 'par'}) = {NaN};

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

% % See if the PAR or the BP sensor do a better job of matching the met station data
% figure(2),clf
% plot(metDat_rt.datetime_utc,metDat_rt.Tair,'.','MarkerSize',10)
% hold on
% plot(parDat_rt.datetime_utc,parDat_rt.Tair,'.r','MarkerSize',4)
% plot(bpDat_rt.datetime_utc,bpDat_rt.Tair,'.g','MarkerSize',4)
% ylabel('Air Temperature (^oC)')
% legend('Gull Met Station','PAR Dataset','Baro Pressure Dataset')
% title('Assessing Air Temperature Datasets')

% Replace missing Gull T_air data with HOBO Baro Pressure T_air data
ind_nan = find(isnan(metDat_rt.Tair));
metDat_rt.Tair(ind_nan) = bpDat_rt.Tair(ind_nan);

% figure(3),clf
% plot(metDat_rt.datetime_utc,metDat_rt.Tair,'.','MarkerSize',4)
% hold on
% plot(metDat_rt.datetime_utc(ind_nan),metDat_rt.Tair(ind_nan),'og','MarkerSize',6,'LineWidth',1)
% ylabel('Air Temperature (^oC)')
% legend('Gull Met Station','Baro Pressure Dataset','location','southeast')
% title('Gap-Filled Air Temperature Data')

%====Gap filling Gull met station atmos p data============================
% figure(4),clf
% plot(metDat_rt.datetime_utc,metDat_rt.patm,'.','MarkerSize',4)
% hold on
% plot(bpDat_rt.datetime_utc,bpDat_rt.patm,'.g','MarkerSize',4)
% ylabel('p_{atm} (hPa)')
% legend('Gull Met Station','Baro Pressure Dataset')
% title('Assessing Atmospheric Pressure Datasets')

% Replace missing Gull p_atm data with HOBO Baro Pressure p_atm data
ind_nan = find(isnan(metDat_rt.patm));
metDat_rt.patm(ind_nan) = bpDat_rt.patm(ind_nan);

% figure(5),clf
% plot(metDat_rt.datetime_utc,metDat_rt.patm,'.','MarkerSize',4)
% hold on
% plot(metDat_rt.datetime_utc(ind_nan),metDat_rt.patm(ind_nan),'og','MarkerSize',6,'LineWidth',1)
% ylabel('p_{atm} (hPa)')
% legend('Gull Met Station','Baro Pressure Dataset')
% title('Gap-Filled Atmospheric Pressure Data')

% Horizontally concatenate sonde and met data (already have common time vector)
dat = synchronize(dat,metDat_rt);
% dat = rmmissing(dat,'DataVariables',{'DO_conc','wspd'});

% Check how many days of data remain in synchronized data
dataPoints = daysact(dat.datetime_utc(1),dat.datetime_utc);
wholeDays = length(unique(floor(dataPoints)));

%====Export data in proper format to run Beck's R code=====================
DO_mgL = dat.DO_conc / 1000 * 31.999;
varNames = ["DateTimeStamp","Temp","Sal","DO_obs","ATemp","BP","WSpd","Tide"];
dat_R = table(dat.datetime_utc,dat.temperature,dat.salinity,DO_mgL,dat.Tair,dat.patm,dat.wspd,dat.depth,'VariableNames',varNames);
dat_R = rmmissing(dat_R);

option = questdlg('Save cleaned data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'\diel-method\owp-data\'])
        switch answer
            case 'Sonde 1'
                writetable(dat_R,[site,'-bc_obs.csv'])
                save([site,'-',sondename,'_obs.mat'],'dat')
            case 'Sonde 2'
                writetable(dat_R,[site,'-erdc_obs.csv'])
                save([site,'-',sondename,'_obs.mat'],'dat')
        end
        disp('Files saved!')
    
    case 'No'
        disp('Files not saved.')
end