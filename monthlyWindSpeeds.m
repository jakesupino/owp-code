%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% monthlyWindSpeeds.m
% This script ...
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% January 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

cd('G:\My Drive\Postdoc\Work\SMIIL\physical-data\wind-speed')

dat1 = readtable('Gull_Met_210914_220625.csv');
varNames = ["datetime_utc","Tair","wspd","wdir","patm","rhumid","source"];
varUnits = ["","degC","m/s","","hPa","",""];
dat1.Properties.VariableNames = varNames;
dat1.Properties.VariableUnits = varUnits;

dat2 = readtable('Gull_Met_WL_220625_220824.xlsx');
varNames = ["datetime_utc","wspd","wdir","Tair","rhumid","patm"];
varUnits = ["","m/s","","degC","","hPa"];
dat2.Properties.VariableNames = varNames;
dat2.Properties.VariableUnits = varUnits;
dat2.datetime_utc = datetime(dat2.datetime_utc);

dat3 = readtable('Gull_Met_WL_220824_221008.xlsx');
dat3.Properties.VariableNames = varNames;
dat3.Properties.VariableUnits = varUnits;
dat3.datetime_utc = datetime(dat3.datetime_utc);

dat4 = readtable('Gull_Met_WL_221008_to_230207.xlsx');
dat4.Properties.VariableNames = varNames;
dat4.Properties.VariableUnits = varUnits;
dat4.datetime_utc = datetime(dat4.datetime_utc);

% Concatenate data tables
dat1 = dat1(:,{'datetime_utc','wspd','wdir','Tair','rhumid','patm'});   % Restructure table 1
metDat = [dat1;dat2;dat3;dat4];
metDat.datetime_utc.TimeZone = 'UTC';

metDat = table2timetable(metDat);

clear dat1 dat2 dat3 dat4 varNames varUnits

figure(1),clf
plot(metDat.datetime_utc,metDat.wspd,'.')
ylabel("Wind speed (m/s)")

save('metData.mat','metDat')

% Find hourly means by month
