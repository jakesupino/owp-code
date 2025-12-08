%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compileDiscreteSamples.m
% This script pulls out the final QC'd AquaTroll salinity, temperature, and
% DO values for the open water platforms at the sample times corresponding to the
% discrete Winkler and DIC/TA samples.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 10/29/2024
% Last updated: 11/4/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
% Load the salinity inventory sheet, which now lives in the shared drive
cd('G:\Shared drives\SMIIL\Shared Data\')
filename = 'SMIIL DIC_TA Sample Salinities - OWP.xlsx';
dic_ta = readtable(filename);
dic_ta.Sampling_Time = datetime(dic_ta.Sampling_Time,'ConvertFrom','datenum','Format',"HH:mm:ss");
dic_ta.datetime_local = dic_ta.Date_Sampled + timeofday(dic_ta.Sampling_Time);
dic_ta.datetime_local.TimeZone = 'America/New_York';
dic_ta.datetime_utc = dic_ta.datetime_local;
dic_ta.datetime_utc.TimeZone = 'UTC';
dic_ta = table2timetable(dic_ta,'RowTimes','datetime_utc');

% Load the Winkler data sheet
filename = 'winklers_owp.csv';
winklers = readtable(filename);
winklers.Datetime_UTC.TimeZone = 'UTC';
winklers = table2timetable(winklers,'RowTimes','Datetime_UTC');

% Load North data and add to structure
cd([rootpath,'open-water-platform-data\north\cleaned\final-qc'])
load('north-cleaned.mat');

tbl_all.north = finalQC;

cd([rootpath,'open-water-platform-data\gull\cleaned\final-qc'])
load('gull-cleaned.mat');
tbl_all.gull = finalQC;

cd([rootpath,'open-water-platform-data\south\cleaned\final-qc'])
load('south-cleaned.mat');
tbl_all.south = finalQC;

clearvars finalQC

%
% dic_ta_new = removevars(dic_ta,{'Bottle_No','Bottle_No_Clean','datetime_local',...
%     'Date_Sampled','Year','Month','Day','Sampling_Time','Field_Temp','Final_Salinity',...
%     'Lab_Salinity','Latitude','Longitude','trip_ID'});
dic_ta_new = removevars(dic_ta,{'Bottle_No_Clean','datetime_local',...
    'Date_Sampled','Year','Month','Day','Sampling_Time','Field_Temp','Final_Salinity',...
    'Lab_Salinity','Latitude','Longitude','trip_ID'});
dic_ta_new.Properties.VariableNames("Location_ID") = "Platform";

winklers_new = removevars(winklers,{'Datetime_EST','S_lab','umol_L_std','x_err'});
winklers_new.Properties.VariableNames("umol_L_avg") = "Winkler_DOconc_avg";

mergedTT = outerjoin(dic_ta_new,winklers_new);

platformName = strings(height(mergedTT),1);
for i = 1:height(mergedTT)
    if mergedTT.Platform_dic_ta_new(i) == ""
        platformName(i) = mergedTT.Platform_winklers_new{i};
    elseif mergedTT.Platform_winklers_new(i) == ""
        platformName(i) = mergedTT.Platform_dic_ta_new{i};
    else
        platformName(i) = mergedTT.Platform_dic_ta_new{i};
    end
end

mergedTT.Platform = platformName;

% Create a new table with all the discrete sample data
varNames = {'datetime_utc','Platform','Year_UTC','Month_UTC','Day_UTC','Time_UTC','Latitude','Longitude','DIC_TA_bottle',...
    'Aquatroll_Press','Aquatroll_Depth','Aquatroll_Temp','Aquatroll_Temp_flag','Aquatroll_Sal','Aquatroll_Sal_flag',....
    'Aquatroll_DOconc','Aquatroll_DOconc_flag','Winkler_DOconc_avg','Winkler_DOconc_flag',...
    'Discrete_DIC','Discrete_DIC_flag','Discrete_TA','Discrete_TA_flag'};
varTypes = {'datetime','string','double','double','double','double','double','double','string',...
    'double','double','double','double','double','double',...
    'double','double','double','double',...
    'double','double','double','double'};
sz = [height(mergedTT) length(varNames)];
discreteSamples = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
discreteSamples.datetime_utc = mergedTT.datetime_utc;

% Fill it in
discreteSamples.Platform = mergedTT.Platform;
discreteSamples.DIC_TA_bottle = mergedTT.Bottle_No; % Add DIC/TA sample bottle numbers
discreteSamples.Year_UTC = year(mergedTT.datetime_utc);
discreteSamples.Month_UTC = month(mergedTT.datetime_utc);
discreteSamples.Day_UTC = day(mergedTT.datetime_utc);
discreteSamples.Time_UTC = timeofday(mergedTT.datetime_utc);

ind.north = find(strcmp("North",discreteSamples.Platform));
ind.gull = find(strcmp("Gull",discreteSamples.Platform));
ind.south = find(strcmp("South",discreteSamples.Platform));

% Define platform lat/lon
lat_north = 39.103164;
lon_north = 74.765492;
lat_gull = 39.072131;
lon_gull = 74.778139;
lat_south = 39.044011;
lon_south = 74.78825;

discreteSamples.Latitude(ind.north) = lat_north;
discreteSamples.Longitude(ind.north) = lon_north;
discreteSamples.Latitude(ind.gull) = lat_gull;
discreteSamples.Longitude(ind.gull) = lon_gull;
discreteSamples.Latitude(ind.south) = lat_south;
discreteSamples.Longitude(ind.south) = lon_south;

% Pull the the Aquatroll values for (pressure,) depth, T, S, and DOconc
fn = fieldnames(tbl_all);
for i = 1:numel(fn)
    % Find QC'd values for closest matching datetimes to discrete samples for each site
    nearest_ind = interp1(tbl_all.(fn{i}).datetime_utc,1:length(tbl_all.(fn{i}).datetime_utc),discreteSamples.datetime_utc(ind.(fn{i})),'nearest');

    % Write the QC'd values into the discrete samples table
    discreteSamples.Aquatroll_Depth(ind.(fn{i})) = tbl_all.(fn{i}).depth(nearest_ind);  
    discreteSamples.Aquatroll_Temp(ind.(fn{i})) = tbl_all.(fn{i}).temperature(nearest_ind);
    discreteSamples.Aquatroll_Sal(ind.(fn{i})) = tbl_all.(fn{i}).salinity(nearest_ind);
    discreteSamples.Aquatroll_DOconc(ind.(fn{i})) = tbl_all.(fn{i}).DOconc(nearest_ind);

end

discreteSamples.Winkler_DOconc_avg = mergedTT.Winkler_DOconc_avg;

discreteSamples.Aquatroll_Press(isnan(discreteSamples.Aquatroll_Press)) = -999;
discreteSamples.Aquatroll_Depth(isnan(discreteSamples.Aquatroll_Depth)) = -999;
discreteSamples.Aquatroll_Temp(isnan(discreteSamples.Aquatroll_Temp)) = -999;
discreteSamples.Aquatroll_Sal(isnan(discreteSamples.Aquatroll_Sal)) = -999;
discreteSamples.Aquatroll_DOconc(isnan(discreteSamples.Aquatroll_DOconc)) = -999;
discreteSamples.Winkler_DOconc_avg(isnan(discreteSamples.Winkler_DOconc_avg)) = -999;

discreteSamples = removevars(discreteSamples,"datetime_utc");

cd('G:\Shared drives\SMIIL\Shared Data\')

writetable(discreteSamples,'discreteSamples_OWP.xlsx')
