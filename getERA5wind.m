%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getERA5wind.m
% This script plots and creates an output table (.mat) of wind speed data downloaded from ERA5.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 5/25/2024
% Last Updated: 8/26/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import wind data downloaded from ERA5=================================
% https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
cd([rootpath,'\physical-data\ERA5'])

fileName = 'all_2021.nc';
timeData1 = ncread(fileName,'time');
timeData1 = double(timeData1);
u10 = ncread(fileName,'u10');
v10 = ncread(fileName,'v10');
for i = 1:length(timeData1)
    u10_1(i,1) = u10(1,1,i);
    v10_1(i,1) = v10(1,1,i);
end

fileName = 'all_2022.nc';
timeData2 = ncread(fileName,'time');
timeData2 = double(timeData2);
u10 = ncread(fileName,'u10');
v10 = ncread(fileName,'v10');
for i = 1:length(timeData2)
    u10_2(i,1) = u10(1,1,i);
    v10_2(i,1) = v10(1,1,i);
end

fileName = 'all_2023.nc';
timeData3 = ncread(fileName,'time');
timeData3 = double(timeData3);
u10 = ncread(fileName,'u10');
v10 = ncread(fileName,'v10');
for i = 1:length(timeData3)
    u10_3(i,1) = u10(1,1,i);
    v10_3(i,1) = v10(1,1,i);
end

fileName = '1-3_2024.nc';
timeData4 = ncread(fileName,'time');
timeData4 = double(timeData4);
u10 = ncread(fileName,'u10');
v10 = ncread(fileName,'v10');
for i = 1:length(timeData4)
    u10_4(i,1) = u10(1,1,i);
    v10_4(i,1) = v10(1,1,i);
end

fileName = '4-6_2024.nc';
timeData5 = ncread(fileName,'time');
timeData5 = double(timeData5);
u10 = ncread(fileName,'u10');
v10 = ncread(fileName,'v10');
for i = 1:length(timeData5)
    u10_5(i,1) = u10(1,1,i);
    v10_5(i,1) = v10(1,1,i);
end

lat = ncread(fileName,'latitude');
lon = ncread(fileName,'longitude');

clearvars timeData u10 v10

u10 = [u10_1;u10_2;u10_3;u10_4;u10_5];
v10 = [v10_1;v10_2;v10_3;v10_4;v10_5];
timeData = [timeData1;timeData2;timeData3;timeData4;timeData5];

%====Convert time unit=====================================================
% (https://memg.ocean.dal.ca/grosse/Teaching/MM2019/Materials/Labs/Lab10/Lab_10.pdf)
timeUnit = ncreadatt(fileName,'time','units');

switch timeUnit(1:3)
    case 'sec'
        timeFac = 1/86400;
    case 'min'
        timeFac = 1/1440;
    case {'hou', 'hrs'}
        timeFac = 1/24;
    case 'day'
        timeFac = 1;
    otherwise
        error('Invalid time unit.');
end

% get the reference time used in the NetCDF file
i = strfind(timeUnit, 'since') + length('since');
refTimeStr = strtrim(timeUnit(i:end));
refTimeNum = datenum(refTimeStr);

% do time conversion
timeData = refTimeNum + timeData*timeFac;

% convert the day counter into a date vector
timeVec = datevec(timeData);
datetime = datetime(timeVec);
datetime.TimeZone = 'UTC';

%====Calculate wind speed from u- and v-components=========================
% (https://confluence.ecmwf.int/pages/viewpage.action?pageId=133262398)
wspd = sqrt(u10.^2 + v10.^2);

%====Create time table of ERA5 data========================================
era5Dat = timetable(datetime,u10,v10,wspd);
era5Dat.Properties.VariableUnits = {'m/s','m/s','m/s'};

%====Plot the ERA5 data====================================================
cd([rootpath,'figures\physical-data'])

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(era5Dat.datetime,era5Dat.wspd,'.','DisplayName','ERA5')
hold off
ylabel('Wind speed (m/s)')

%====Option to save data===================================================
option = questdlg('Save data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'physical-data\ERA5'])
        save('era5.mat','era5Dat')
        disp('Files saved!')
    case 'No'
        disp('Files not saved.')
end
