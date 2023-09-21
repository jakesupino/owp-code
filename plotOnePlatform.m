%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotOnePlatform.m
% This script reads in and plots the data for one deployment from one sonde
% from one open-water platform.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% 9/14/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in sonde data

clear all;close all;clc

site = 'gull';  % CHANGE THIS

cd(['H:\My Drive\Postdoc\SMIIL\raw-data\open-water-platform-data\',site])

% Interactively select files (from same deployment)
[fileNames,dataPath] = uigetfile('*.csv','Select one or more files','MultiSelect','on');
fileName1 = extractBefore(fileNames{1},'.csv');
fileName2 = extractBefore(fileNames{2},'.csv');

% Extract the deployment number from the filename and make a column to add to data table
depNum = extractBetween(fileName1,"dep","-");
depNum = str2double(depNum);

% Extract the site from the filepath for plotting later
depSite = extractBetween(dataPath,'platform-data\','\');
depSite = [upper(depSite{1}(1)),depSite{1}(2:end)];

%----For BC sonde----------------------------------------------------------
sonde1 = readtable([dataPath,fileName1],'TextType','string');

if depNum == 5
    paramNames1 = ["datetime_local","datetime_utc","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","chla",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
elseif depNum == 6
    paramNames1 = ["datetime_local","datetime_utc","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","chla","pH","pH_raw","ORP",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
elseif 7 <= depNum <= 14
    paramNames1 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","chla",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
end

sonde1.Properties.VariableNames = paramNames1;

% Remove rows with missing datetime, NaT (headers; also some files have extra lines at end)
sonde1 = sonde1(~any(ismissing(sonde1(:,1)),2),:);

% Create a column for datetime in UTC
sonde1.datetime_local.TimeZone = 'America/New_York';
datetime_utc = datetime(sonde1.datetime_local,'TimeZone','UTC');
datetime_utc = table(datetime_utc,'VariableNames',"datetime_utc");

% Create a column for the deployment number
depNumVec = repmat(depNum,height(sonde1),1);
depNumTable = table(depNumVec,'VariableNames',"deployment");

if 5 <= depNum <= 6
    sonde1 = [depNumTable,sonde1];
elseif 7 <= depNum <= 14
    sonde1 = [depNumTable,datetime_utc,sonde1];
end

%----For ERDC sonde--------------------------------------------------------
sonde2 = readtable([dataPath,fileName2],'TextType','string');

if depNum == 5
    paramNames2 = ["datetime_local","datetime_utc","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
elseif depNum == 6
    paramNames2 = ["datetime_local","datetime_utc","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
elseif 7 <= depNum <= 14
    paramNames2 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
end

sonde2.Properties.VariableNames = paramNames2;

% Remove rows with missing datetime, NaT (some files have extra lines at end)
sonde2 = sonde2(~any(ismissing(sonde2(:,1)),2),:);

% Create a column for datetime in UTC
sonde2.datetime_local.TimeZone = 'America/New_York';
datetime_utc = datetime(sonde2.datetime_local,'TimeZone','UTC');
datetime_utc = table(datetime_utc,'VariableNames',"datetime_utc");

% Create a column for the deployment number
depNumVec = repmat(depNum,height(sonde2),1);
depNumTable = table(depNumVec,'VariableNames',"deployment");

if 5 <= depNum <= 6
    sonde2 = [depNumTable,sonde2];
elseif 7 <= depNum <= 14
    sonde2 = [depNumTable,datetime_utc,sonde2];
end

%% Convert units

% Convert [ft] to [m]
sonde1.depth = sonde1.depth/3.281;
sonde2.depth = sonde2.depth/3.281;

% Convert [mg/L] to [umol/L]
sonde1.DO_conc = sonde1.DO_conc/31.999*10^3;
sonde2.DO_conc = sonde2.DO_conc/31.999*10^3;

paramUnits1 = ["","","","uS/cm","uS/cm","psu","ohm cm",...
    "g/cm3","ppt","umol/L","%sat","torr","","mV","mV","RFU",...
    "degC","V","%","mmHg","psi","m"];

sonde1.Properties.VariableUnits = paramUnits1;

paramUnits2 = ["","","","uS/cm","uS/cm","psu","ohm cm",...
    "g/cm3","ppt","umol/L","%sat","torr","","mV","mV","NTU",...
    "degC","V","%","mmHg","psi","m"];

sonde2.Properties.VariableUnits = paramUnits2;

%% Plot data for both sondes together

red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde

figure,clf

tl = tiledlayout(4,2,'TileSpacing','Tight');
title(tl,[depSite,' Deployment #',num2str(depNum),' - Both Sondes'])
xlabel(tl,'Time (UTC)')

nexttile
plot(sonde1.datetime_utc,sonde1.depth,'color',red)
hold on
plot(sonde2.datetime_utc,sonde2.depth,'color',blue)
hold off
title('Depth')
ylabel('Depth (m)')
xlim('tight')

% Add a global legend to label sondes
lg = legend('BC','ERDC');
lg.Layout.Tile = 'north';

nexttile
plot(sonde1.datetime_utc,sonde1.pH,'color',red)
hold on
plot(sonde2.datetime_utc,sonde2.pH,'color',blue)
title('pH')
ylabel('pH')
xlim('tight')

nexttile
plot(sonde1.datetime_utc,sonde1.temperature,'color',red)
hold on
plot(sonde2.datetime_utc,sonde2.temperature,'color',blue)
title('Temperature')
ylabel('Temperature (^oC)')
xlim('tight')

nexttile
plot(sonde1.datetime_utc,sonde1.DO_conc,'color',red)
hold on
plot(sonde2.datetime_utc,sonde2.DO_conc,'color',blue)
hold off
title('DO')
ylabel('DO (\mumol/L)')
xlim('tight')

nexttile
plot(sonde1.datetime_utc,sonde1.salinity,'color',red)
hold on
plot(sonde2.datetime_utc,sonde2.salinity,'color',blue)
hold off
title('Salinity')
ylabel('Salinity (PSU)')
xlim('tight')

nexttile
plot(sonde1.datetime_utc,sonde1.ORP,'color',red)
hold on
plot(sonde2.datetime_utc,sonde2.ORP,'color',blue)
hold off
title('ORP')
ylabel('ORP (mV)')
xlim('tight')

nexttile
plot(sonde2.datetime_utc,sonde2.turbidity,'Color',blue)
title('Turbidity')
ylabel('Turbidity (NTU)')
xlim('tight')

nexttile
plot(sonde1.datetime_utc,sonde1.chla,'color',red)
title('Chl a')
ylabel('Chl a (RFU)')
xlim('tight')

cd(['H:\My Drive\Postdoc\SMIIL\figures\open-water-platform-figures\',site])

%% Save created tables in .mat files

cd(dataPath)

saveFileName = extractBefore(fileName1,'-bc');

save([saveFileName,'.mat'],"sonde1","sonde2")