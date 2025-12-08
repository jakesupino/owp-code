%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pHValidation.m
% This script plots the discrete DIC/TA measurements on top of the AquaTroll
% pH data to help with visual assessment of whether the BC, ERDC, or mean
% is the "best guess" for this deployment (noted manually). Based on this
% assessment, the best guess calls are entered in "dupCheck_DO_pH.m".
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 8/15/2024
% Last updated: 11/4/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

%==========================================================================
%   Calculate pH from DIC and TA via CO2SYS
%==========================================================================

%----Load the discrete DIC/TA data-----------------------------------------
cd('G:\Shared drives\SMIIL\Shared Data\')
sample_log = readtable('SMIIL DIC_TA Sample Salinities - OWP.xlsx');

cd('G:\Shared drives\SMIIL\Shared Data\Processed')
load('DIC_TA_Output_SMIIL_OWP_KF');    % Downloaded from Lab_Sample_Processing/SMIIL GitHub

% Make sample log into timetable
date = datetime(sample_log.Date_Sampled);
timeOfDay = datetime(sample_log.Sampling_Time,'ConvertFrom','datenum');
datetime_local = date + timeofday(timeOfDay);
datetime_local.TimeZone = 'America/New_York';
datetime_utc = datetime_local;
datetime_utc.TimeZone = 'UTC';
sample_log.datetime_utc = datetime_utc;
sample_log = table2timetable(sample_log,'RowTimes','datetime_utc');

% Make data output into timetable
datetime_utc = datetime(DIC_TA_discrete.dn_UTC,'ConvertFrom','datenum');
datetime_utc.TimeZone = 'UTC';
DIC_TA_discrete.datetime_utc = datetime_utc;
DIC_TA_discrete = table2timetable(DIC_TA_discrete,"RowTimes","datetime_utc");
DIC_TA_discrete = sortrows(DIC_TA_discrete);

%----Calculate [H+] from discrete DIC/TA samples-----------------------------
par1type =    1; % Type "1" = "alkalinity"
par1     =    DIC_TA_discrete.TA; % Values of the first parameter
par2type =    2; % Type "2" = "DIC"
par2     =    DIC_TA_discrete.DIC; % Values of the second parameter
sal      =    sample_log.Final_Salinity; % Salinity of the sample
tempin   =    sample_log.Field_Temp; % In situ temperature
presin   =    0; % Lab pressure
tempout  =    0; % Doesn't matter here
presout  =    1; % In situ pressure
sil      =   50; % Concentration of silicate  in the sample (in umol/kg)
po4      =    2; % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c    =    4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% Do the calculation. See CO2SYS's help for syntax and output format
[carb_sys.results,headers,units] = CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);

carb_sys.results = array2table(carb_sys.results);
carb_sys.results.Properties.VariableNames = headers';

%----Estimate error--------------------------------------------------------
% For error in DIC, use the max of DIC_std and DIC_drift (ORIGINAL DATA FILE)
% for i = 1:length(DIC_TA_discrete.DIC)
%     eDIC(i,1) = max(DIC_TA_discrete.DIC_std(i),DIC_TA_discrete.DIC_drift(i)); 
% end

epar1     =    DIC_TA_discrete.TA_std; % Error for the first parameter
% epar2     =    eDIC; % Error for the second parameter
epar2     =    DIC_TA_discrete.DIC_std;
esal      =    0;
etemp     =    0;
esil      =    0;
epo4      =    0;
epk       =    0;
ebt       =    0;
er        =    0;
% Do the calculation. See CO2SYS error.m help for syntax and output format
[err, headers, units] = errors(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,...
                                    epar1,epar2,esal,etemp,esil,epo4,epk,ebt,er,...
                                    pHscale,k1k2c,kso4c);

carb_sys.errors = array2table(err);
carb_sys.errors.Properties.VariableNames = headers';
carb_sys.errors.Properties.VariableUnits = units';

% Make a table with just the discrete [H+] values and errors
Hfree = carb_sys.results.("Hfreeout")/10^6; % Convert Hfreeout from umol/kg to mol/kg
uHfree = carb_sys.errors.("u(Hout)")/10^9; % Convert u(Hout) from nmol/kg to to mol/kg
H_calc = table(DIC_TA_discrete.datetime_utc,DIC_TA_discrete.Platform,Hfree,uHfree);
H_calc.Properties.VariableNames = {'datetime_utc','platform','Hfree','error'};

% Extract the discrete values for the site being analyzed
ind_site = find(strcmp(site,H_calc.platform));

%==========================================================================
%   Compare discrete sample and AquaTroll [H+]
%==========================================================================
% Load the final QC'd AquaTroll data (best guess data)
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load('finalQC.mat');

% Load the cleaned BC and ERDC data for comparison
cd([rootpath,'open-water-platform-data\',site,'\cleaned\movmed\'])
load([site,'-bc-cleaned'])
load([site,'-erdc-cleaned'])
dat1 = sonde1_cleaned;
dat2 = sonde2_cleaned;
% Synchronize the BC and ERDC data to a common datetime vector
dat_syn = synchronize(dat1,dat2);

clearvars dat1 dat2 sonde1_cleaned sonde2_cleaned

% Convert sensor pH to [H+]
finalQC.Hfree = 10.^-finalQC.pH;
dat_syn.Hfree1 = 10.^-dat_syn.pH_dat1;
dat_syn.Hfree2 = 10.^-dat_syn.pH_dat2;

% Global plotting settings
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

cd([rootpath,'figures\open-water-platform\',site,'\validation\dic_ta'])

% Plot AquaTroll pH data with discrete [H+] on top
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(dat_syn.datetime_utc,dat_syn.Hfree1,'o','MarkerSize',8,'color',red,'DisplayName','BC')
hold on
plot(dat_syn.datetime_utc,dat_syn.Hfree2,'o','MarkerSize',8,'color',blue,'DisplayName','ERDC')
plot(finalQC.datetime_utc,finalQC.Hfree,'.k','MarkerSize',12,'DisplayName','"Best Guess"')
errorbar(H_calc.datetime_utc(ind_site),H_calc.Hfree(ind_site),H_calc.error(ind_site),'o','color',rgb('gold'),'MarkerSize',8,'LineWidth',2,'DisplayName','Discrete Sample')
xline(finalQC.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
xlim([dt1 dt2])    % Use same x limits
xlabel('UTC')
ylabel('[H+] (mol/kg)')
title(site,'FontSize',fontsize)
legend('show','location','best')
ylim([0 6E-8])
%%
% % Linear regression
% % Remove rows with missing pH values
% finalQC_trimmed = rmmissing(finalQC,'DataVariables',"pH");
% 
% % Find closest matching indices
% t1 = H_calc.datetime_utc;
% t2 = finalQC_trimmed.datetime_utc;
% ind_sonde = interp1(t2,1:length(t2),t1,'nearest','extrap');
% 
% ind_discrete = 1:1:height(H_calc);
% 
% 
% % If closest matching samples are more than 25 minutes off, do not include
% for i = 1:length(ind_sonde(~isnan(ind_sonde)))
%     dt = between(finalQC_trimmed.datetime_utc(ind_sonde(i)),H_calc.datetime_utc(i));
%     if abs(time(dt)) > minutes(25)
%         ind_sonde(i) = NaN;
%         ind_discrete(i) = NaN;
%     end
% end
% 
% ind_sonde = rmmissing(ind_sonde);
% ind_discrete = rmmissing(ind_discrete);
% 
% % Create the model
% x = finalQC_trimmed.Hfree(ind_sonde);
% y = H_calc.Hfree(ind_discrete);
% mdl = fitlm(x,y,'y~x1-1');  % Force intercept through zero; see Wilkinson notation
% 
% tbl = table(finalQC_trimmed.datetime_utc(ind_sonde),H_calc.datetime_utc(ind_discrete),x,y);
% tbl.Properties.VariableNames = {'AquaTroll datetime','Discrete datetime','AquaTroll [H+]','Discrete [H+]'};
% 
% % Create equation string for plot
% eqn = ['y = ',num2str(mdl.Coefficients.Estimate,3),'x'];
% R2 = num2str(mdl.Rsquared.Ordinary,2);
% 
% fig2 = figure(2);clf
% fig2.WindowState = 'maximized';
% h = plot(mdl,'marker','.','markersize',20);
% % delete(h([3 4]))    % Delete confidence bounds on plot
% hold on
% plot([min(min([x,y])) max(max([x,y]))], [min(min([x,y])) max(max([x,y]))],'--k')
% errorbar(x,y,H_calc.error(ind_discrete),'.b','MarkerSize',dotsize,'LineWidth',2)
% % x1=0:0.1:8.5;
% % y1=mdl.Coefficients.Estimate*x1;
% % plot(x1,y1,'-r')
% % plot([min(min([x1,y1])) max(max([x1,y1]))], [min(min([x1,y1])) max(max([x1,y1]))],'--k')
% legend('',[eqn,newline,'R^2 = ',R2],'95% confidence interval','','1:1 line','Location','southeast')
% xlabel('AquaTroll [H+] (mol/kg)')
% ylabel('Discrete [H+] (mol/kg)')
% title(site)
% daspect([1 1 1])
% ylim([min(min([x,y])) max(max([x,y]))])
% xlim([min(min([x,y])) max(max([x,y]))])

% % See Table 2, Rule 16 for propagation of uncertainty for y = A + log10(x)
% % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3387884/
% abs_uncert = carb_sys.errors.("u(Hout)")*10^-9; % Factor of 10^-9 to convert nmol/kg to mol/kg
% abs_val = 10.^(-carb_sys.results.pHout);        % [H+] in mol/kg
% rel_uncert = 0.4343 * abs_uncert ./ abs_val;  % Unitless
% pHerr = rel_uncert .* carb_sys.results.("pHout");    % pH units

%====Save the plots========================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\open-water-platform\',site,'\validation\dic_ta'])
        saveas(fig1,'pHValid_timeseries.png')
        saveas(fig1,'pHValid_timeseries.fig')
        % saveas(fig2,[site,'_linreg.png'])
        % saveas(fig2,[site,'_linreg.fig'])
        disp('Plots saved!')
    case 'No'
        disp('Plots not saved.')
end
