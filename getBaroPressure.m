%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getBaroPressure.m
% This script creates an output table (.mat) of data from the HOBO Water
% Level Logger (which measures atmospheric pressure and air temperature) located at TWI.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% January 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

%====Import atm p dataset measured at TWI==================================
varNames = ["sample","datetime_local","patm","Tair"];
varUnits = ["","","kPa","degC"];

folder = 'G:\My Drive\Postdoc\Work\SMIIL\physical-data\baro-pressure';
cd(folder)
myFiles = dir(fullfile(folder,'*.csv'));

bpDat = array2timetable(NaN(1,3),'RowTimes',datetime("today"));
bpDat.Properties.VariableNames = ["sample","patm","Tair"];

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
    bpDat = [bpDat;dat];
end

bpDat.datetime_local.TimeZone = 'America/New_York';
bpDat(1,:) = [];    % Remove the first row from initializing the timetable
bpDat = removevars(bpDat,{'sample'});
datetime_utc = table(datetime(bpDat.datetime_local,'TimeZone','UTC'),'VariableNames',"datetime_utc");
datetime_utc = table2timetable(datetime_utc);
bpDat = [datetime_utc,bpDat];
%====Sort timetable by date================================================
bpDat = sortrows(bpDat);

% Convert p_atm units from [kPa] to [hPa]
bpDat.patm = bpDat.patm .* 10;    % [hPa]
bpDat.Properties.VariableUnits = ["hPa","degC"];

%====Retime timetable======================================================
newTimeStep = minutes(10);
bpDat = retime(bpDat,'regular','mean','TimeStep',newTimeStep); % Calculate mean of values in each time bin
bpDat.datetime_utc = bpDat.datetime_utc + newTimeStep/2;

%====Plot air temperature and PAR data=====================================
t = tiledlayout(2,1);
ax1 = nexttile;
plot(bpDat.datetime_utc,bpDat.Tair,'.','markersize',4);ylabel('Air Temperature (^oC)')
ax2 = nexttile;
plot(bpDat.datetime_utc,bpDat.patm,'.','markersize',4);ylabel('Atmospheric Pressure (hPa)')
title(t,'HOBO Water Level Logger at TWI','fontsize',16)
t.TileSpacing = 'compact';

% Link the axes
linkaxes([ax1,ax2],'x')     

%====Save data=============================================================
cd('G:\My Drive\Postdoc\Work\SMIIL\physical-data\baro-pressure')
save('baroPress.mat','bpDat')