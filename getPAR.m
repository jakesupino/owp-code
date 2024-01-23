%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getPAR.m
% This script...
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% January 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

cd('G:\My Drive\Postdoc\Work\SMIIL\physical-data\PAR')

%====Import PAR dataset measured at TWI====================================
% lux and ppfd units - https://www.apogeeinstruments.com/conversion-ppfd-to-lux/

dat1 = readtable('PAR_SetAug2022.xlsx');
varNames = ["sample","datetime_local","Tair","light_lux","par"];
varUnits = ["","","degC","lx","umol m-2 s-1"];
dat1.Properties.VariableNames = varNames;
dat1.Properties.VariableUnits = varUnits;

dat2 = readtable('PAR_SepttoDec.xlsx');
dat2.Properties.VariableNames = varNames;
dat2.Properties.VariableUnits = varUnits;

dat3 = readtable('TWI 2023-10-19 19_05_10 EDT (Data EDT).xlsx');
dat3.Properties.VariableNames = varNames;
dat3.Properties.VariableUnits = varUnits;

parDat = [dat1;dat2;dat3];
parDat.datetime_local.TimeZone = 'America/New_York';
parDat = removevars(parDat,{'sample'});
datetime_utc = table(datetime(parDat.datetime_local,'TimeZone','UTC'),'VariableNames',"datetime_utc");
parDat = [datetime_utc,parDat];

% Retime timetable
parDat = table2timetable(parDat,'RowTimes','datetime_utc');
newTimeStep = minutes(10);
parDat = retime(parDat,'regular','mean','TimeStep',newTimeStep); % Calculate mean of values in each time bin

parDat.datetime_utc = parDat.datetime_utc + newTimeStep/2;
parDat.datetime_local = parDat.datetime_utc;
parDat.datetime_local.TimeZone = 'America/New_York';

% Plot air temperature and PAR data
t = tiledlayout(2,1);
ax1 = nexttile;
plot(parDat.datetime_local,parDat.Tair,'.','markersize',4);ylabel('Air Temperature (^oC)')
ax2 = nexttile;
plot(parDat.datetime_local,parDat.par,'.','markersize',4);ylabel('PAR (\mumol m^{-2} s^{-1})')
title(t,'HOBO sensor at TWI','fontsize',16)
t.TileSpacing = 'compact';

% Link the axes
linkaxes([ax1,ax2],'x')     

cd('G:\My Drive\Postdoc\Work\SMIIL\physical-data\PAR')
save('par.mat','parDat')