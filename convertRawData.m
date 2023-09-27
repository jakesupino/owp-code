%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convertRawData.m
% This script reads in the raw .csv data from one open-water platform for one deployment from both
% sondes (BC and ERDC).
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% 9/14/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in sonde data

clear all;close all;clc

site = 'gull';  % CHANGE THIS

cd(['H:\My Drive\Postdoc\SMIIL\raw-data\open-water-platform-data\',site])

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
elseif (depNum >= 7) && (depNum <= 14) 
    paramNames1 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","chla",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
     paramNames2 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
end

%----For BC sonde----------------------------------------------------------
sonde1 = readtable([dataPath,fileName1],'TextType','string');
sonde1.Properties.VariableNames = paramNames1;

% Remove rows with missing datetime, NaT (headers; also, some files have extra lines at end)
sonde1 = sonde1(~any(ismissing(sonde1(:,1)),2),:);

% Create a column for datetime in UTC
if (depNum >= 1) && (depNum <= 2)
    % Do nothing (datetime already in UTC)
elseif (depNum >= 5) && (depNum <= 14)
    sonde1.datetime_local.TimeZone = 'America/New_York';
    datetime_utc = datetime(sonde1.datetime_local,'TimeZone','UTC');
    datetime_utc = table(datetime_utc,'VariableNames',"datetime_utc");
end

% Create a column for the deployment number
depNumCol = array2table(repmat(depNum,height(sonde1),1),'VariableNames',"deployment");
sonde1 = [depNumCol,sonde1];

if (depNum >= 1) && (depNum <= 2) || (depNum >= 5) && (depNum <= 6)
    sonde1.datetime_utc.TimeZone = 'UTC';
elseif (depNum >= 7) && (depNum <= 14)
    sonde1 = [datetime_utc,sonde1];
end

%----For ERDC sonde--------------------------------------------------------
sonde2 = readtable([dataPath,fileName2],'TextType','string');

sonde2.Properties.VariableNames = paramNames2;

% Remove rows with missing datetime, NaT (headers; also, some files have extra lines at end)
sonde2 = sonde2(~any(ismissing(sonde2(:,1)),2),:);

% Create a column for datetime in UTC
if (depNum >= 1) && (depNum <= 2) 
    % Do nothing (datetime already in UTC)
elseif (depNum >= 5) && (depNum <= 14)
    sonde2.datetime_local.TimeZone = 'America/New_York';
    datetime_utc = datetime(sonde2.datetime_local,'TimeZone','UTC');
    datetime_utc = table(datetime_utc,'VariableNames',"datetime_utc");
end

% Create a column for the deployment number
depNumCol = array2table(repmat(depNum,height(sonde2),1),'VariableNames',"deployment");
sonde2 = [depNumCol,sonde2];

if (depNum >= 1) && (depNum <= 2) || (depNum >= 5) && (depNum <= 6)
    sonde2.datetime_utc.TimeZone = 'UTC';
elseif (depNum >= 7) && (depNum <= 14)
    sonde2 = [datetime_utc,sonde2];
end

% Create columns filled with NaN for parameters not measured
if depNum == 1 || depNum == 2
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
        if depNum == 2
        DO_sat = array2table(NaN(height(sonde1),1),'VariableNames',"DO_sat");
        sonde1 = [DO_sat sonde1];
    end
    
    datetime_local = array2table(NaT(height(sonde2),1,'TimeZone','America/New_York'),'VariableNames',"datetime_local");
    actual_cond = array2table(NaN(height(sonde2),1),'VariableNames',"actual_cond");
    specific_cond = array2table(NaN(height(sonde2),1),'VariableNames',"specific_cond");
    resistivity = array2table(NaN(height(sonde2),1),'VariableNames',"resistivity");
    density = array2table(NaN(height(sonde2),1),'VariableNames',"density");
    barometric_p = array2table(NaN(height(sonde2),1),'VariableNames',"barometric_p");
    p = array2table(NaN(height(sonde2),1),'VariableNames',"p");
    TDS = array2table(NaN(height(sonde2),1),'VariableNames',"TDS");
    pO2 = array2table(NaN(height(sonde2),1),'VariableNames',"pO2");
    pH = array2table(NaN(height(sonde2),1),'VariableNames',"pH");
    pH_raw = array2table(NaN(height(sonde2),1),'VariableNames',"pH_raw");
    ORP = array2table(NaN(height(sonde2),1),'VariableNames',"ORP");
    external_voltage = array2table(NaN(height(sonde2),1),'VariableNames',"external_voltage");
    battery_capacity = array2table(NaN(height(sonde2),1),'VariableNames',"battery_capacity");
    sonde2 = [datetime_local actual_cond specific_cond resistivity density barometric_p p TDS pO2 pH_raw ORP external_voltage battery_capacity sonde2];
    if depNum == 2
        DO_sat = array2table(NaN(height(sonde2),1),'VariableNames',"DO_sat");
        sonde2 = [DO_sat sonde2];
    end
    
elseif depNum >= 5
    nitrate = array2table(NaN(height(sonde1),1),'VariableNames',"nitrate");
    sonde1 = [nitrate sonde1];
end

%----Convert units---------------------------------------------------------

% Convert [ft] to [m]
if depNum == 9 || depNum == 13
    sonde2.depth = sonde2.depth/3.281;
elseif depNum == 11
    sonde1.depth = sonde1.depth/3.281;
elseif depNum == 14
    sonde1.depth = sonde1.depth/3.281;
    sonde2.depth = sonde2.depth/3.281;
end

% Convert [mg/L] to [umol/L]
sonde1.DO_conc = sonde1.DO_conc/31.999*10^3;
sonde2.DO_conc = sonde2.DO_conc/31.999*10^3;

% Restructure table so parameter columns are in same order between deployments

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

%----Save created tables in .mat files-------------------------------------
cd(dataPath)

saveFileName = extractBefore(fileName1,'-bc');

save([saveFileName,'.mat'],"sonde1","sonde2")