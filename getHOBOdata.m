%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getHOBOdata.m
% This script creates an output table (.mat) of data from the HOBO
% water-level logger at TWI (which measures atmospheric pressure and air temperature).
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: January 2024
% Last updated: 5/30/204
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
% HOBO water-level logger
folder = [rootpath,'physical-data\HOBO'];
cd(folder)
myFiles = dir(fullfile(folder,'*.csv'));

hoboDat = array2timetable(NaN(1,3),'RowTimes',datetime("today"));
hoboDat.Properties.VariableNames = ["sample","patm","Tair"];

varNames = ["sample","datetime_local","patm","Tair"];
varUnits = ["","","kPa","degC"];

for i = 1:length(myFiles)
    file = myFiles(i).name;
    fprintf(1,'Now reading %s\n',file);
    opt = detectImportOptions(file);
    opt.VariableNames = varNames;
    opt = setvaropts(opt,{'datetime_local'},'InputFormat','dd/MM/yyyy HH:mm:ss');
    dat = readtable(file);
    dat.Properties.VariableNames = varNames;
    dat.Properties.VariableUnits = varUnits;
    dat = table2timetable(dat);
    hoboDat = [hoboDat;dat];
end

hoboDat.datetime_local.TimeZone = 'America/New_York';
hoboDat(1,:) = [];    % Remove the first row from initializing the timetable
hoboDat = removevars(hoboDat,{'sample'});
datetime_utc = table(datetime(hoboDat.datetime_local,'TimeZone','UTC'),'VariableNames',"datetime_utc");
datetime_utc = table2timetable(datetime_utc);
hoboDat = [datetime_utc,hoboDat];

%====Option to save data===================================================
option = questdlg('Save data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'physical-data\HOBO'])
        save('HOBO.mat','hoboDat')
        disp('Files saved!')
    case 'No'
        disp('Files not saved.')
end