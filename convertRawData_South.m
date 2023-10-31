%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convertRawData_South.m
% This script reads in the raw .csv data from the South open-water platform 
% for one deployment from both sondes (BC and ERDC).
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% 10/42023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===Read in sonde data=====================================================
clear all;close all;clc

site = 'south';

cd(['G:\My Drive\Postdoc\SMIIL\open-water-platform-data\',site,'\original\deployments'])

% Interactively select files (from same deployment)
[fileNames,dataPath] = uigetfile('*.csv','MultiSelect','on');

if numel(fileNames) == 2
    fileName1 = extractBefore(fileNames{1},'.csv');
    fileName2 = extractBefore(fileNames{2},'.csv');
    depNum = extractBetween(fileName1,"dep","-");
else
    sonde = extractBefore(extractAfter(fileNames,13),5);
    switch sonde
        case{'bc'}
        fileName1 = extractBefore(fileNames,'.csv');
        depNum = extractBetween(fileName1,"dep","-");
        case{'erdc'}
        fileName2 = extractBefore(fileNames,'.csv');
        depNum = extractBetween(fileName2,"dep","-");
    end
end
depNum = str2double(depNum);

% Extract the site from the filepath for plotting later
depSite = extractBetween(dataPath,'platform-data\','\');
depSite = [upper(depSite{1}(1)),depSite{1}(2:end)];

% Name parameters based on order in .csv file
switch depNum
    case{1, 2}
    paramNames1 = ["datetime_utc","depth","temperature","salinity","chla","nitrate","DO_conc"];
    paramNames2 = ["datetime_utc","depth","temperature","salinity","turbidity","pH","DO_conc"];
    
    case{4}
    paramNames1 = ["datetime_local","nitrate","nitrate_raw","specific_cond","salinity",...
        "DO_conc","DO_sat","pO2","chla","temperature","external_voltage","depth"];
    
    case{5}
    paramNames1 = ["datetime_local","datetime_utc","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","chla",...
        "temperature","barometric_p","p","depth"];
    paramNames2 = ["datetime_local","datetime_utc","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];

    case{6}
    paramNames1 = ["datetime_utc","depth","temperature","salinity",...
        "specific_cond","chla","pH","pH_raw","DO_conc","DO_sat"];

    case{7}
    paramNames1 = ["datetime_local","chla","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
    paramNames2 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
    
    case{9}
    paramNames1 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","chla",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
    paramNames2 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","chla",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
    
    case{8,10,11,12,13,14}
    paramNames1 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","chla",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
    paramNames2 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
        "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
        "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];

    case{15}
        paramNames2 = ["datetime_local","actual_cond","specific_cond","salinity","resistivity",...
            "density","TDS","DO_conc","DO_sat","pO2","pH","pH_raw","ORP","turbidity",...
            "temperature","external_voltage","battery_capacity","barometric_p","p","depth"];
end

if numel(fileNames) == 2
    % sonde 1 = BC sonde
    sonde1 = readtable([dataPath,fileName1],'TextType','string');
    sonde1.Properties.VariableNames = paramNames1;
    % sonde 2 = ERDC sonde
    sonde2 = readtable([dataPath,fileName2],'TextType','string');
    sonde2.Properties.VariableNames = paramNames2;
else
    switch sonde
        case{'bc'}
            sonde1 = readtable([dataPath,fileName1],'TextType','string');
            sonde1.Properties.VariableNames = paramNames1;
        case{'erdc'}
            sonde2 = readtable([dataPath,fileName2],'TextType','string');
            sonde2.Properties.VariableNames = paramNames2;
    end
end

%===Adjust properties======================================================
if numel(fileNames) == 2
    % Remove rows with missing datetime, NaT (headers; also, some files have log notes at end)
    sonde1 = sonde1(~any(ismissing(sonde1(:,1)),2),:);
    sonde2 = sonde2(~any(ismissing(sonde2(:,1)),2),:);
    % Create a column for the deployment number
    depNumCol = array2table(repmat(depNum,height(sonde1),1),'VariableNames',"deployment");
    sonde1 = [depNumCol,sonde1];
    depNumCol = array2table(repmat(depNum,height(sonde2),1),'VariableNames',"deployment");
    sonde2 = [depNumCol,sonde2];
else
    switch sonde
        case{'bc'}
            % Remove rows with missing datetime, NaT (headers; also, some files have log notes at end)
            sonde1 = sonde1(~any(ismissing(sonde1(:,1)),2),:);
            % Create a column for the deployment number
            depNumCol = array2table(repmat(depNum,height(sonde1),1),'VariableNames',"deployment");
            sonde1 = [depNumCol,sonde1];
        case{'erdc'}
            % Remove rows with missing datetime, NaT (headers; also, some files have log notes at end)
            sonde2 = sonde2(~any(ismissing(sonde2(:,1)),2),:);
            % Create a column for the deployment number
            depNumCol = array2table(repmat(depNum,height(sonde2),1),'VariableNames',"deployment");
            sonde2 = [depNumCol,sonde2];
    end
end

% Set timezones. If not present, create a column for datetime in UTC.
switch depNum
    case{1,2}
    sonde1.datetime_utc.TimeZone = 'UTC';
    sonde2.datetime_utc.TimeZone = 'UTC';

    case{4}
    sonde1.datetime_local.TimeZone = 'America/New_York';
    datetime_utc1 = datetime(sonde1.datetime_local,'TimeZone','UTC');
    datetime_utc1 = table(datetime_utc1,'VariableNames',"datetime_utc");
    sonde1 = [datetime_utc1,sonde1];

    case{5}
    sonde1.datetime_local.TimeZone = 'America/New_York';
    sonde2.datetime_local.TimeZone = 'America/New_York';
    sonde1.datetime_utc.TimeZone = 'UTC';
    datetime_utc2 = datetime(sonde2.datetime_local,'TimeZone','UTC');
    sonde2.datetime_utc = datetime_utc2;    % Had a datetime UTC column, but empty of values

    case{6}
    sonde1.datetime_utc.TimeZone = 'America/New_York';

    case{7,8,9,10,11,12,13,14}
    sonde1.datetime_local.TimeZone = 'America/New_York';
    datetime_utc1 = datetime(sonde1.datetime_local,'TimeZone','UTC');
    datetime_utc1 = table(datetime_utc1,'VariableNames',"datetime_utc");
    sonde1 = [datetime_utc1,sonde1];
    sonde2.datetime_local.TimeZone = 'America/New_York';
    datetime_utc2 = datetime(sonde2.datetime_local,'TimeZone','UTC');
    datetime_utc2 = table(datetime_utc2,'VariableNames',"datetime_utc");
    sonde2 = [datetime_utc2,sonde2];

    case{15}
    sonde2.datetime_local.TimeZone = 'America/New_York';
    datetime_utc2 = datetime(sonde2.datetime_local,'TimeZone','UTC');
    datetime_utc2 = table(datetime_utc2,'VariableNames',"datetime_utc");
    sonde2 = [datetime_utc2,sonde2];

end

% Create columns filled with NaNs for parameters not measured
switch depNum
    case {1,2}
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
    DO_sat = array2table(NaN(height(sonde1),1),'VariableNames',"DO_sat");
    pH = array2table(NaN(height(sonde1),1),'VariableNames',"pH");
    pH_raw = array2table(NaN(height(sonde1),1),'VariableNames',"pH_raw");
    ORP = array2table(NaN(height(sonde1),1),'VariableNames',"ORP");
    external_voltage = array2table(NaN(height(sonde1),1),'VariableNames',"external_voltage");
    battery_capacity = array2table(NaN(height(sonde1),1),'VariableNames',"battery_capacity");
    sonde1 = [datetime_local actual_cond specific_cond resistivity density barometric_p p TDS pO2 DO_sat pH pH_raw ORP external_voltage battery_capacity sonde1];
    
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
    DO_sat = array2table(NaN(height(sonde2),1),'VariableNames',"DO_sat");
    pH_raw = array2table(NaN(height(sonde2),1),'VariableNames',"pH_raw");
    ORP = array2table(NaN(height(sonde2),1),'VariableNames',"ORP");
    external_voltage = array2table(NaN(height(sonde2),1),'VariableNames',"external_voltage");
    battery_capacity = array2table(NaN(height(sonde2),1),'VariableNames',"battery_capacity");
    sonde2 = [datetime_local actual_cond specific_cond resistivity density barometric_p p TDS pO2 DO_sat pH_raw ORP external_voltage battery_capacity sonde2];

    case 4
    % sonde1 (BC)
    actual_cond = array2table(NaN(height(sonde1),1),'VariableNames',"actual_cond");
    resistivity = array2table(NaN(height(sonde1),1),'VariableNames',"resistivity");
    density = array2table(NaN(height(sonde1),1),'VariableNames',"density");
    barometric_p = array2table(NaN(height(sonde1),1),'VariableNames',"barometric_p");
    p = array2table(NaN(height(sonde1),1),'VariableNames',"p");
    TDS = array2table(NaN(height(sonde1),1),'VariableNames',"TDS");
    pH = array2table(NaN(height(sonde1),1),'VariableNames',"pH");
    pH_raw = array2table(NaN(height(sonde1),1),'VariableNames',"pH_raw");
    ORP = array2table(NaN(height(sonde1),1),'VariableNames',"ORP");
    battery_capacity = array2table(NaN(height(sonde1),1),'VariableNames',"battery_capacity");
    sonde1 = [actual_cond resistivity density barometric_p p TDS pH pH_raw ORP battery_capacity sonde1];
    
    case 5
    % sonde1 (BC)
    pH = array2table(NaN(height(sonde1),1),'VariableNames',"pH");
    pH_raw = array2table(NaN(height(sonde1),1),'VariableNames',"pH_raw");
    ORP = array2table(NaN(height(sonde1),1),'VariableNames',"ORP");
    external_voltage = array2table(NaN(height(sonde1),1),'VariableNames',"external_voltage");
    battery_capacity = array2table(NaN(height(sonde1),1),'VariableNames',"battery_capacity");
    nitrate = array2table(NaN(height(sonde1),1),'VariableNames',"nitrate");
    sonde1 = [pH pH_raw ORP external_voltage battery_capacity nitrate sonde1];
    
    case 6
    % sonde1 (BC)
    datetime_local = array2table(NaT(height(sonde1),1,'TimeZone','America/New_York'),'VariableNames',"datetime_local");
    actual_cond = array2table(NaN(height(sonde1),1),'VariableNames',"actual_cond");
    resistivity = array2table(NaN(height(sonde1),1),'VariableNames',"resistivity");
    density = array2table(NaN(height(sonde1),1),'VariableNames',"density");
    barometric_p = array2table(NaN(height(sonde1),1),'VariableNames',"barometric_p");
    p = array2table(NaN(height(sonde1),1),'VariableNames',"p");
    TDS = array2table(NaN(height(sonde1),1),'VariableNames',"TDS");
    pO2 = array2table(NaN(height(sonde1),1),'VariableNames',"pO2");
    ORP = array2table(NaN(height(sonde1),1),'VariableNames',"ORP");
    external_voltage = array2table(NaN(height(sonde1),1),'VariableNames',"external_voltage");
    battery_capacity = array2table(NaN(height(sonde1),1),'VariableNames',"battery_capacity");
    nitrate = array2table(NaN(height(sonde1),1),'VariableNames',"nitrate");
    sonde1 = [datetime_local actual_cond resistivity density barometric_p p TDS pO2 ORP external_voltage battery_capacity nitrate sonde1];
   
    case 9
    % sonde1 (BC)
    nitrate = array2table(NaN(height(sonde1),1),'VariableNames',"nitrate");
    sonde1 = [nitrate sonde1];
    % sonde2 (ERDC)
    turbidity = array2table(NaN(height(sonde2),1),'VariableNames',"turbidity");
    sonde2 = [turbidity sonde2];

    case {7,8,10,11,12,13,14}
    % sonde1 (BC)
    nitrate = array2table(NaN(height(sonde1),1),'VariableNames',"nitrate");
    sonde1 = [nitrate sonde1];
end

%===Convert units==========================================================

% Depth: Convert [ft] to [m]
switch depNum
    case 9
        sonde1.depth = sonde1.depth/3.281;
    case {10,11}
        sonde2.depth = sonde2.depth/3.281;
    case{12,13,14}
        sonde1.depth = sonde1.depth/3.281;
        sonde2.depth = sonde2.depth/3.281;
end

% DO concentration: Convert [mg/L] to [umol/L]
if numel(fileNames) == 2
    sonde1.DO_conc = sonde1.DO_conc/31.999*10^3;
    sonde2.DO_conc = sonde2.DO_conc/31.999*10^3;
else
    switch sonde
        case{'bc'}
            sonde1.DO_conc = sonde1.DO_conc/31.999*10^3;
        case{'erdc'}
            sonde2.DO_conc = sonde2.DO_conc/31.999*10^3;
    end
end

% Pressure: Convert [mbar] to [psi]
switch depNum
    case {10,11,12,14}
        sonde1.p = sonde1.p/68.9476;
    case 8
        sonde1.p = sonde1.p/68.9476;
        sonde2.p = sonde2.p/68.9476;
end

%===Restructure tables=====================================================
%===so parameter columns are in same order between deployments=============
if numel(fileNames) == 2
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
else
    switch sonde
        case{'bc'}
            sonde1 = sonde1(:,{'deployment' 'datetime_utc' 'datetime_local' ...
                'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
                'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
                'pH' 'pH_raw' 'ORP' 'chla' 'nitrate' 'external_voltage' 'battery_capacity'});
            paramUnits1 = ["","","",...
                "uS/cm","uS/cm","psu","ohm cm","g/cm3",...
                "degC","mmHg","psi","m","ppt","umol/L","%sat","torr",...
                "","mV","mV","mg/L","RFU","V","%"];
            sonde1.Properties.VariableUnits = paramUnits1;
        case{'erdc'}
            sonde2 = sonde2(:,{'deployment' 'datetime_utc' 'datetime_local' ...
                'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
                'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
                'pH' 'pH_raw' 'ORP' 'turbidity' 'external_voltage' 'battery_capacity'});
            paramUnits2 = ["","","",...
                "uS/cm","uS/cm","psu","ohm cm","g/cm3",...
                "degC","mmHg","psi","m","ppt","umol/L","%sat","torr",...
                "","mV","mV","NTU","V","%"];
            sonde2.Properties.VariableUnits = paramUnits2;
    end
end


%====Save created tables in .mat files=====================================
option = questdlg(['Save .mat file in SMIIL\open-water-platform-data\',site,'\original\deployments?'],'Save File','Y','N','Y');

switch option
    case 'Y'
        if numel(fileNames) == 2
            saveFileName = extractBefore(fileName1,'-bc');
            save([saveFileName,'.mat'],"sonde1","sonde2")
        else
            switch sonde
                case{'bc'}
                    saveFileName = extractBefore(fileName1,'-bc');
                    save([saveFileName,'.mat'],"sonde1")
                case{'erdc'}
                    saveFileName = extractBefore(fileName2,'-erdc');
                    save([saveFileName,'.mat'],"sonde2")
            end
        end
        disp('File saved!')
    case 'N'
        disp('File not saved.')
end