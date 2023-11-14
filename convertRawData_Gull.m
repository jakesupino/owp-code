%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convertRawData_gull.m
% This script reads in the raw .csv data from the Gull open-water platform
% for one deployment from both sondes (BC and ERDC).
%
% AUTHOR:
% Emily Chua
%
% DATE:
% 9/14/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===Read in sonde data=====================================================
clear all;close all;clc

site = 'gull';

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

cd([rootpath,'open-water-platform-data\',site,'\original\deployments'])

% Interactively select files (from same deployment)
[fileNames,dataPath] = uigetfile('*.csv','MultiSelect','on');
fileName1 = extractBefore(fileNames{1},'.csv');
fileName2 = extractBefore(fileNames{2},'.csv');

% Extract the deployment number from the filename and make a column to add to data table
depNum = extractBetween(fileName1,"dep","-");
depNum = str2double(depNum);

% Extract the site from the filepath for plotting later
depSite = extractBetween(dataPath,'platform-data\','\');
depSite = [upper(depSite{1}(1)),depSite{1}(2:end)];

% Name parameters based on order in .csv file
if depNum == 1
    paramNames1 = ["datetime_utc","depth","temperature","salinity","chla","nitrate","DO_conc","DO_sat"];
    paramNames2 = ["datetime_utc","depth","temperature","salinity","turbidity","pH","DO_conc","DO_sat"];
elseif depNum == 2
    paramNames1 = ["datetime_utc","depth","temperature","salinity","chla","nitrate","DO_conc"];
    paramNames2 = ["datetime_utc","depth","temperature","salinity","turbidity","pH","DO_conc"];
elseif depNum == 5
    paramNames1 = ["datetime_local","datetime_utc","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","chla",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
    paramNames2 = ["datetime_local","datetime_utc","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
elseif depNum == 6
    paramNames1 = ["datetime_local","datetime_utc","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","chla","pH","pH_raw","ORP",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
    paramNames2 = ["datetime_local","datetime_utc","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
elseif (depNum >= 7) && (depNum <= 15)
    paramNames1 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","chla",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
    paramNames2 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
end

% sonde 1 = BC sonde
sonde1 = readtable([dataPath,fileName1],'TextType','string');
sonde1.Properties.VariableNames = paramNames1;

% sonde 2 = ERDC sonde
sonde2 = readtable([dataPath,fileName2],'TextType','string');
sonde2.Properties.VariableNames = paramNames2;

%===Adjust properties======================================================
% Remove rows with missing datetime, NaT (headers; also, some files have log notes at end)
sonde1 = sonde1(~any(ismissing(sonde1(:,1)),2),:);
sonde2 = sonde2(~any(ismissing(sonde2(:,1)),2),:);

% Set timezones. If not present, create a column for datetime in UTC.
if (depNum >= 1) && (depNum <= 2)
    sonde1.datetime_utc.TimeZone = 'UTC';
    sonde2.datetime_utc.TimeZone = 'UTC';
elseif (depNum >= 5) && (depNum <= 6)
    sonde1.datetime_local.TimeZone = 'America/New_York';
    sonde1.datetime_utc.TimeZone = 'UTC';
    sonde2.datetime_local.TimeZone = 'America/New_York';
    sonde2.datetime_utc.TimeZone = 'UTC';
elseif (depNum >= 7) && (depNum <= 15)
    sonde1.datetime_local.TimeZone = 'America/New_York';
    datetime_utc1 = datetime(sonde1.datetime_local,'TimeZone','UTC');
    datetime_utc1 = table(datetime_utc1,'VariableNames',"datetime_utc");
    sonde1 = [datetime_utc1,sonde1];
    sonde2.datetime_local.TimeZone = 'America/New_York';
    datetime_utc2 = datetime(sonde2.datetime_local,'TimeZone','UTC');
    datetime_utc2 = table(datetime_utc2,'VariableNames',"datetime_utc");
    sonde2 = [datetime_utc2,sonde2];
end

% Create a column for the deployment number
depNumCol = array2table(repmat(depNum,height(sonde1),1),'VariableNames',"deployment");
sonde1 = [depNumCol,sonde1];

depNumCol = array2table(repmat(depNum,height(sonde2),1),'VariableNames',"deployment");
sonde2 = [depNumCol,sonde2];

% Create columns filled with NaNs for parameters not measured
if depNum == 1 || depNum == 2
    % sonde1 (BC)
    datetime_local = array2table(NaT(height(sonde1),1,'TimeZone','America/New_York'),'VariableNames',"datetime_local");
    actual_cond = array2table(NaN(height(sonde1),1),'VariableNames',"actual_cond");
    specific_cond = array2table(NaN(height(sonde1),1),'VariableNames',"specific_cond");
    resistivity = array2table(NaN(height(sonde1),1),'VariableNames',"resistivity");
    density = array2table(NaN(height(sonde1),1),'VariableNames',"density");
    barometric_p = array2table(NaN(height(sonde1),1),'VariableNames',"barometric_p");
    p = array2table(NaN(height(sonde1),1),'VariableNames',"p");
    TDS = array2table(NaN(height(sonde1),1),'VariableNames',"TDS");
    pO2 = array2table(NaN(height(sonde1),1),'VariableNames',"pO2");
    pH = array2table(NaN(height(sonde1),1),'VariableNames',"pH");
    pH_raw = array2table(NaN(height(sonde1),1),'VariableNames',"pH_raw");
    ORP = array2table(NaN(height(sonde1),1),'VariableNames',"ORP");
    external_voltage = array2table(NaN(height(sonde1),1),'VariableNames',"external_voltage");
    battery_capacity = array2table(NaN(height(sonde1),1),'VariableNames',"battery_capacity");
    sonde1 = [datetime_local actual_cond specific_cond resistivity density barometric_p p TDS pO2 pH pH_raw ORP external_voltage battery_capacity sonde1];

    % sonde2 (ERDC)
    datetime_local = array2table(NaT(height(sonde2),1,'TimeZone','America/New_York'),'VariableNames',"datetime_local");
    actual_cond = array2table(NaN(height(sonde2),1),'VariableNames',"actual_cond");
    specific_cond = array2table(NaN(height(sonde2),1),'VariableNames',"specific_cond");
    resistivity = array2table(NaN(height(sonde2),1),'VariableNames',"resistivity");
    density = array2table(NaN(height(sonde2),1),'VariableNames',"density");
    barometric_p = array2table(NaN(height(sonde2),1),'VariableNames',"barometric_p");
    p = array2table(NaN(height(sonde2),1),'VariableNames',"p");
    TDS = array2table(NaN(height(sonde2),1),'VariableNames',"TDS");
    pO2 = array2table(NaN(height(sonde2),1),'VariableNames',"pO2");
    pH_raw = array2table(NaN(height(sonde2),1),'VariableNames',"pH_raw");
    ORP = array2table(NaN(height(sonde2),1),'VariableNames',"ORP");
    external_voltage = array2table(NaN(height(sonde2),1),'VariableNames',"external_voltage");
    battery_capacity = array2table(NaN(height(sonde2),1),'VariableNames',"battery_capacity");
    sonde2 = [datetime_local actual_cond specific_cond resistivity density barometric_p p TDS pO2 pH_raw ORP external_voltage battery_capacity sonde2];
    
    if depNum == 2
        DO_sat = array2table(NaN(height(sonde1),1),'VariableNames',"DO_sat");
        sonde1 = [DO_sat sonde1];
        DO_sat = array2table(NaN(height(sonde2),1),'VariableNames',"DO_sat");
        sonde2 = [DO_sat sonde2];
    end

elseif (depNum >= 5) && (depNum <= 15)
    nitrate = array2table(NaN(height(sonde1),1),'VariableNames',"nitrate");
    sonde1 = [nitrate sonde1];
end

%===Convert units==========================================================

% Depth: Convert [ft] to [m]
if depNum == 11
    sonde1.depth = sonde1.depth/3.281;
elseif depNum == 9 || depNum == 13
    sonde2.depth = sonde2.depth/3.281;
elseif depNum == 14 || depNum == 15
    sonde1.depth = sonde1.depth/3.281;
    sonde2.depth = sonde2.depth/3.281;
end

% DO concentration: Convert [mg/L] to [umol/L]
sonde1.DO_conc = sonde1.DO_conc/31.999*10^3;
sonde2.DO_conc = sonde2.DO_conc/31.999*10^3;

% Pressure: Convert [mbar] to [psi]
if depNum == 9 || depNum == 12
    sonde1.p = sonde1.p/68.9476;
elseif depNum == 11
    sonde2.p = sonde2.p/68.9476;
end

%===Restructure tables=====================================================
%===so parameter columns are in same order between deployments=============
sonde1 = sonde1(:,{'deployment' 'datetime_utc' 'datetime_local' ...
    'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
    'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
    'pH' 'pH_raw' 'ORP' 'chla' 'nitrate' 'external_voltage' 'battery_capacity'});

sonde2 = sonde2(:,{'deployment' 'datetime_utc' 'datetime_local' ...
    'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
    'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
    'pH' 'pH_raw' 'ORP' 'turbidity' 'external_voltage' 'battery_capacity'});

paramUnits1 = ["","","",...
    "uS/cm","uS/cm","psu","ohm cm","g/cm3",...
    "degC","mmHg","psi","m","ppt","umol/L","%sat","torr",...
    "","mV","mV","mg/L","RFU","V","%"];

paramUnits2 = ["","","",...
    "uS/cm","uS/cm","psu","ohm cm","g/cm3",...
    "degC","mmHg","psi","m","ppt","umol/L","%sat","torr",...
    "","mV","mV","NTU","V","%"];

sonde1.Properties.VariableUnits = paramUnits1;

sonde2.Properties.VariableUnits = paramUnits2;

%====Save created tables in .mat files=====================================
option = questdlg(['Save .mat file in SMIIL\open-water-platform-data\',site,'\original\deployments?'],'Save File','Y','N','Y');

switch option
    case 'Y'
        saveFileName = extractBefore(fileName1,'-bc');
        save([saveFileName,'.mat'],"sonde1","sonde2")
        disp('File saved!')
    case 'N'
        disp('File not saved.')
end