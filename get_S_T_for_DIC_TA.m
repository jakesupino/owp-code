%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_S_T_for_DIC_TA.m
% This script pulls out the final QC'd AquaTroll salinity and temperature
% values for the open water platforms at the sample times corresponding to the
% discrete DIC and TA samples.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 7/31/2024
% Last updated: 9/19/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';   

%====Import data===========================================================
% Load the salinity inventory sheet, which now lives in the shared drive
cd('G:\Shared drives\SMIIL\Shared Data\')
% filename = 'SMIIL DIC_TA Sample Salinities.xlsx';
filename = 'SMIIL DIC_TA Sample Salinities - OWP.xlsx';
sampleLog = readtable(filename);

% Load North data and add to structure
cd([rootpath,'open-water-platform-data\north\cleaned\dupcheck'])
load('north-bestGuess.mat');
tbl_all.north = [bestguess.salinity, bestguess.temperature];
% Remove any rows that contain a NaN
idx = find(isnan(tbl_all.north.salinity));
tbl_all.north(idx,:) = [];

% Load Gull data and add to structure
cd([rootpath,'open-water-platform-data\gull\cleaned\dupcheck'])
load('gull-bestGuess.mat');
tbl_all.gull = [bestguess.salinity, bestguess.temperature];
% Remove any rows that contain a NaN for salinity
idx = find(isnan(tbl_all.gull.salinity));
tbl_all.gull(idx,:) = [];

% Load South data and add to structure
cd([rootpath,'open-water-platform-data\south\cleaned\dupcheck'])
load('south-bestGuess.mat');
tbl_all.south = [bestguess.salinity, bestguess.temperature];
% Remove any rows that contain a NaN for salinity
idx = find(isnan(tbl_all.south.salinity));
tbl_all.south(idx,:) = [];

clearvars bestguess

% Create datetimes for all samples in local and UTC format
T  = datetime(sampleLog.Sampling_Time,'ConvertFrom','datenum','Format','hh:mm:ss');
datetime_local = sampleLog.Date_Sampled + timeofday(T);
datetime_local.TimeZone = 'America/New_York';
datetime_utc = datetime_local;
datetime_utc.TimeZone = 'UTC';

% Pull out rows corresponding to North, Gull, and South
ind.north = find(strcmp('North',sampleLog.Location_ID));
ind.gull = find(strcmp('Gull',sampleLog.Location_ID));
ind.south = find(strcmp('South',sampleLog.Location_ID));

fn = fieldnames(tbl_all);

for i = 1:numel(fn)
    % Find QC'd S and T values for closest matching datetimes to discrete samples for each site
    nearest_ind = interp1(tbl_all.(fn{i}).datetime_utc,1:length(tbl_all.(fn{i}).datetime_utc),datetime_utc(ind.(fn{i})),'nearest');

    % Overwrite the QC'd S and T value into the 'sampleLog' table
    sampleLog.Final_Salinity(ind.(fn{i})) = tbl_all.(fn{i}).salinity(nearest_ind);  % Take the AquaTroll best-guess S as the "Final Salinity"
    sampleLog.Field_Temp(ind.(fn{i})) = tbl_all.(fn{i}).temperature(nearest_ind);
end

% Write table to Excel spreadsheet
cd('G:\Shared drives\SMIIL\Shared Data\')
writetable(sampleLog,filename)