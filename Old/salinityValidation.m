%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% salinityValidation.m
% This script uses the discrete DIC/TA sample salinities to validate the 
% salinity data.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 9/26/2024
% Last updated: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

% Import BC and ERDC salinity
cd([rootpath,'open-water-platform-data\',site,'\cleaned\movmed\'])
load([site,'-bc-cleaned.mat'])
load([site,'-erdc-cleaned.mat'])

%====Import best-guess salinity====================================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc\'])
load([site,'-cleaned.mat'])

%====Import discrete data==================================================
cd('G:\Shared drives\SMIIL\Shared Data')
discreteS = readtable('SMIIL DIC_TA Sample Salinities - OWP.xlsx');

varUnits = ["","","","","","","","","","psu","psu","psu","degC","",""];
discreteS.Properties.VariableUnits = varUnits;

discreteS.Sampling_Time = datetime(discreteS.Sampling_Time,'ConvertFrom','datenum');
myDatetime = discreteS.Date_Sampled + timeofday(discreteS.Sampling_Time);
myDatetime.TimeZone = 'America/New_York';
datetime_utc = myDatetime;
datetime_utc.TimeZone = 'UTC';
discreteS.datetime_utc = datetime_utc;
discreteS = movevars(discreteS,'datetime_utc','Before','trip_ID');

% Find indices of samples taken from this site
ind_discrete = find(ismember(discreteS.Location_ID,site));

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = finalQC.datetime_utc(end);
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 6;
circlesize = 6;

switch site
    case 'Gull'
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
            'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
            'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
            'Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
    case 'North'
        label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8',...
            'Deployment 9','Deployment 10','Deployment 11','Deployment 12',...
            'Deployment 13','Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
    case 'South'
        label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
            'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
            'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
            'Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
end

% Find indices of deployment changes
ind_dep = find(diff(finalQC.deployment) > 0);
switch site
    case 'Gull'
        dep = table([1;finalQC.deployment(ind_dep+1)],[1;ind_dep+1]);
    case 'North'
        dep = table([2;finalQC.deployment(ind_dep+1)],[1;ind_dep+1]);
    case 'South'
        dep = table([1;finalQC.deployment(ind_dep+1)],[1;ind_dep+1]);
end
dep.Properties.VariableNames = {'depNum','ind'};

cd([rootpath,'figures\open-water-platform\',site,'\validation\dic_ta\salinity'])

%====Plot time series datasets with discrete samples on top================
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(sonde1_cleaned.datetime_utc,sonde1_cleaned.salinity,'.','Color',red,'DisplayName','BC - Cleaned S')
hold on
plot(sonde2_cleaned.datetime_utc,sonde2_cleaned.salinity,'.','Color',blue,'DisplayName','ERDC - Cleaned S')
plot(finalQC.datetime_utc,finalQC.salinity,'.k','MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','Best-guess AquaTroll S')
plot(discreteS.datetime_utc(ind_discrete),discreteS.Lab_Salinity(ind_discrete),'o','Color',rgb('gold'),'MarkerSize',6,'LineWidth',2,'DisplayName','Lab S')
xline(finalQC.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
xlim([dt1 dt2])                 % Use same x limits
ylabel('Salinity (PSU)')
title([site,' - Salinity Validation'])
set(gca,'FontSize',fontsize)
legend('show','Location','best')
