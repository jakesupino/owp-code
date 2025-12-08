%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% salinityLinReg.m
% This script performs linear regressions between the "best-guess" salinity 
% data (moving median) and the lab salinity values measured from the DIC/TA
% samples.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 10/2/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

% %====Import re-calculated DO concentrations===============================
% cd([rootpath,'open-water-platform-data\',site,'\cleaned\corrected'])
% load([site,'-DOCorr.mat'])

%====Import best-guess sonde data==========================================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc\'])
load('finalQC.mat')

%====Import lab S from discrete samples====================================
cd('G:\Shared drives\SMIIL\Shared Data')
discreteS = readtable('SMIIL DIC_TA Sample Salinities - OWP.xlsx');

varUnits = ["","","","","","","","","","psu","psu","degC","",""];
discreteS.Properties.VariableUnits = varUnits;

discreteS.Sampling_Time = datetime(discreteS.Sampling_Time,'ConvertFrom','datenum');
myDatetime = discreteS.Date_Sampled + timeofday(discreteS.Sampling_Time);
myDatetime.TimeZone = 'America/New_York';
datetime_utc = myDatetime;
datetime_utc.TimeZone = 'UTC';
discreteS.datetime_utc = datetime_utc;
discreteS = movevars(discreteS,'datetime_utc','Before','trip_ID');

% Find indices of samples taken from this site
ind_site = find(ismember(discreteS.Location_ID,site));

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

% cd([rootpath,'figures\open-water-platform\',site,'\validation\winklers'])

% %====Plot time series datasets with discrete samples on top==============
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(finalQC.datetime_utc,finalQC.salinity,'.k','MarkerSize',dotsize,'DisplayName','"Best Guess"')
hold on
plot(discreteS.datetime_utc(ind_site),discreteS.Lab_Salinity(ind_site),'o','Color',rgb('gold'),'MarkerSize',6,'LineWidth',2,'DisplayName','Lab S')
xline(finalQC.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Salinity (psu)')
title([site,' - "Best Guess" Salinity'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)
legend('show')
%%
%====Best-guess linear regression==========================================
Sbestguess = table(finalQC.datetime_utc,finalQC.salinity);
Sbestguess.Properties.VariableNames = {'datetime_utc','salinity'};
Sbestguess = rmmissing(Sbestguess);

t1 = discreteS.datetime_utc(ind_site);
t2 = Sbestguess.datetime_utc;
ind_sonde = interp1(t2,1:length(t2),t1,'nearest','extrap');

% If closest matching samples are more than 25 minutes off, do not include
for i = 1:length(ind_sonde(~isnan(ind_sonde)))
    dt = between(Sbestguess.datetime_utc(ind_sonde(i)),discreteS.datetime_utc(ind_site(i)));
    if abs(time(dt)) > minutes(25)
        ind_sonde(i) = NaN;
    end
end

ind_matching = table(ind_site,ind_sonde);
ind_matching = rmmissing(ind_matching);

% "Best guess" linear model
x = Sbestguess.salinity(ind_matching.ind_sonde);
y = discreteS.Lab_Salinity(ind_matching.ind_site);
mdl = fitlm(x,y,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl.Rsquared.Ordinary,2);

fig2 = figure(2);clf
fig2.WindowState = 'maximized';
h = plot(mdl,'marker','.','markersize',20);
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([0 ylim], [0 ylim],'--k')
% errorbar(x,y,wink.DO_std(ind_matching.ind_wink),'.b','MarkerSize',dotsize,'LineWidth',2)
legend('Data',[eqn,newline,'R^2 = ',R2],'1:1 line','Location','southeast')
xlabel('"Best-Guess" Salinity (PSU)','interpreter','tex')                                                                                                                                     
ylabel('Lab Salinity (PSU)','interpreter','tex')    
title([site,' Best Guess Salinity Validation'])
pbaspect([1 1 1]); % Make relative lengths of axes equal
xlim([min(min([x y])) max(max([x y]))])
ylim([min(min([x y])) max(max([x y]))])

%====Save regression model results=========================================
option = questdlg('Save regression model?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\',site,'\cleaned\validation'])
        save([site,'-discreteS_mdl.mat'],'mdl')
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end

%====Save the plots========================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\open-water-platform\',site,'\validation\dic_ta\salinity'])
        saveas(fig1,'S_bestguess-timeseries.fig')
        saveas(fig1,'S_bestguess-timeseries.png')
        saveas(fig2,'S_bestguess-linreg.fig')
        saveas(fig2,'S_bestguess-linreg.png')   
        disp('Plots saved!')
    case 'No'
        disp('Plots not saved.')
end
