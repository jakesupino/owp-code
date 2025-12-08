%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combineSondeData.m
% This script combines the data that have been QC'd up to this point for the 
% BC and ERDC sondes and plots the results.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 5/3/2024
% Last updated: 11/20/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

%====Import sonde data===================================================
% Best-guess data for depth, S, T, DO conc, and pH
cd([rootpath,'open-water-platform-data\',site,'\cleaned\dupcheck\'])
load([site,'-bestGuess.mat'])

% Turbidity and chlorophyll a
cd([rootpath,'open-water-platform-data\',site,'\cleaned\movmed\'])
load([site,'-bc-cleaned'])
load([site,'-erdc-cleaned'])
dat1 = sonde1_cleaned;
dat2 = sonde2_cleaned;
% Synchronize the BC and ERDC data to a common datetime vector
dat_syn = synchronize(dat1,dat2);

clearvars dat1 dat2 sonde1_cleaned sonde2_cleaned

%====Import re-calculated DO concentrations=============================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\corrected'])
load([site,'-DOCorr.mat'])

%====Import Winkler data================================================
cd('G:\Shared drives\SMIIL\Shared Data')
wink = readtable('winklers_owp.csv');
varNames = ["datetime_utc","datetime_local","platform","S_lab","DO_mean","DO_std","DO_%err"];
varUnits = ["","","","psu","umol/L","umol/L","%"];
wink.Properties.VariableNames = varNames;
wink.Properties.VariableUnits = varUnits;
wink.datetime_local.TimeZone = 'America/New_York';
wink.datetime_utc.TimeZone = 'UTC';
wink = table2timetable(wink);
% Find indices of samples taken from this site
ind_wink = find(ismember(wink.platform,site));

%====Import windspeed data==============================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

cd([rootpath,'physical-data\final-dataset'])
load('par.mat')

%====Create final QC data table============================================
var1 = bestguess.deployment;
var2 = bestguess.depth;
var3 = bestguess.salinity;
var4 = bestguess.temperature;
var5 = bestguess.DOconc;
var6 = bestguess.pH;
var7 = timetable(dat_syn.datetime_utc,dat_syn.chla,'VariableNames',{'chla'});
var8 = timetable(dat_syn.datetime_utc,dat_syn.turbidity,'VariableNames',{'turbidity'});

sondes_comb = synchronize(var1,var2,var3,var4,var5,var6,var7,var8);

%====Calculate DO %sat and add to final QC table===========================
DOeqbm = O2sol(sondes_comb.salinity,sondes_comb.temperature);    % [umol/kg]
% Convert units - density function from https://www.mbari.org/technology/matlab-scripts/oceanographic-calculations/
rho = density(sondes_comb.salinity,sondes_comb.temperature);    % [g/mL]
DOeqbm = DOeqbm.*rho;    % [umol/L]

DOsat = sondes_comb.DOconc./DOeqbm * 100;  % [%]

sondes_comb = addvars(sondes_comb,DOsat,'After','DOconc');
sondes_comb.Properties.VariableUnits = {'','m','psu','degC','umol/L','%sat','','RFU','NTU'};

%====Find the indices of deployment changes================================
ind_dep = find(diff(sondes_comb.deployment) > 0);
switch site
    case 'Gull'
        dep = table([1;sondes_comb.deployment(ind_dep+1)],[1;ind_dep+1]);
    case 'North'
        dep = table([2;sondes_comb.deployment(ind_dep+1)],[1;ind_dep+1]);
    case 'South'
        dep = table([1;sondes_comb.deployment(ind_dep+1)],[1;ind_dep+1]);
end
dep.Properties.VariableNames = {'depNum','ind'};

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

%====Global plotting settings===========================================
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 6;
circlesize = 6;

%====Save the cleaned data=================================================
option = questdlg('Save combined data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\',site,'\cleaned\combined'])
        save([site,'-combined.mat'],'sondes_comb')
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end

% %====Plot all parameters in linked figure==================================
% % fig1 = figure('units','inches','position',[5 9 14 10]);
% fig1 = figure(1);
% 
% t = tiledlayout(4,2,'TileSpacing','compact');
% 
% ax1 = nexttile;
% plot(sondes_comb.datetime_utc,sondes_comb.DOconc,'.k','MarkerSize',dotsize,'LineWidth',linewidth)
% hold on
% errorbar(wink.datetime_utc(ind_wink),wink.DO_mean(ind_wink),wink.DO_std(ind_wink),'.','Color',rgb('gold'),'MarkerSize',dotsize,'LineWidth',2,'DisplayName','Winkler')
% xline(sondes_comb.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('DO Conc (\mumol/L)')
% 
% ax2 = nexttile;
% plot(sondes_comb.datetime_utc,sondes_comb.salinity,'.k','MarkerSize',dotsize,'LineWidth',linewidth)
% xline(sondes_comb.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('Salinity (PSU)')
% 
% ax3 = nexttile;
% plot(sondes_comb.datetime_utc,sondes_comb.temperature,'.k','MarkerSize',dotsize,'LineWidth',linewidth)
% hold on
% xline(sondes_comb.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('T (^oC)')
% ylim([-5 35])
% 
% ax4 = nexttile;
% plot(sondes_comb.datetime_utc,sondes_comb.pH,'.k','MarkerSize',dotsize,'LineWidth',linewidth)
% xline(sondes_comb.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('pH')
% 
% ax5 = nexttile;
% plot(sondes_comb.datetime_utc,sondes_comb.turbidity,'.k','MarkerSize',dotsize,'LineWidth',linewidth)
% xline(sondes_comb.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('Turbidity (NTU)')
% 
% ax6 = nexttile;
% plot(era5Dat.datetime,era5Dat.wspd,'.k','MarkerSize',dotsize,'LineWidth',linewidth)
% xline(sondes_comb.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('Wind Speed (m/s)')
% 
% ax7 = nexttile;
% plot(sondes_comb.datetime_utc,sondes_comb.chla,'.k','MarkerSize',dotsize,'LineWidth',linewidth)
% xline(sondes_comb.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('Chl a (RFU)')
% 
% ax8 = nexttile;
% plot(parDat.datetime_utc,parDat.par,'.k','MarkerSize',dotsize,'LineWidth',linewidth)
% xline(sondes_comb.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('PAR (\mumol m^{-2} s^{-1})')
% 
% linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'x')
% xlim([sondes_comb.datetime_utc(1) sondes_comb.datetime_utc(end)])
% title(t,[site,' - Syncrhonized Sonde Data'],'FontSize',18)
