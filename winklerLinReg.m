%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% winklerLinReg.m
% This script uses the Winkler samples to validate the DO concentration
% data. Linear regressions are performed for the recalculated BC and ERDC 
% data, as well as the overall "best-guess" data.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 5/9/2024
% Last updated: 5/23/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

%====Import re-calculated DO concentrations================================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\corrected'])
load([site,'-DOCorr.mat'])

%====Import best-guess DO concentration====================================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc\'])
load('finalQC.mat')

%====Import Winkler data===================================================
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

cd([rootpath,'figures\open-water-platform\',site,'\validation\winklers'])

%====Plot time series datasets with Winklers on top========================
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(DOcorr.bc.datetime_utc,DOcorr.bc.DOconc,'.','Color',rgb('darkred'),'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','Recalculated BC')
hold on
plot(DOcorr.erdc.datetime_utc,DOcorr.erdc.DOconc,'.','Color',rgb('darkblue'),'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','Recalculated ERDC')
plot(finalQC.datetime_utc,finalQC.DOconc,'.k','MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','"Best Guess"')
errorbar(wink.datetime_utc(ind_wink),wink.DO_mean(ind_wink),wink.DO_std(ind_wink),'.','Color',rgb('gold'),'MarkerSize',12,'LineWidth',2,'DisplayName','Winkler')
xline(finalQC.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
xlim([dt1 dt2])                 % Use same x limits
ylim([0 450])
ylabel('DO Conc (\mumol/L)')
title([site,' - "Best Guess" DO Concentration'])
set(gca,'FontSize',fontsize)
legend('show','Location','best')

%====BC linear regression==================================================
DObc = table(DOcorr.bc.datetime_utc,DOcorr.bc.DOconc);
DObc.Properties.VariableNames = {'datetime_utc','DOconc'};
DObc = rmmissing(DObc);

% Find closest matching indices
t1 = wink.datetime_utc(ind_wink);
t2 = DObc.datetime_utc;
ind_sonde = interp1(t2,1:length(t2),t1,'nearest','extrap');

% If closest matching samples are more than 25 minutes off, do not include
for i = 1:length(ind_sonde(~isnan(ind_sonde)))
    dt = between(DObc.datetime_utc(ind_sonde(i)),wink.datetime_utc(ind_wink(i)));
    if abs(time(dt)) > minutes(25)
        ind_sonde(i) = NaN;
    end
end

ind_matching = table(ind_wink,ind_sonde);
ind_matching = rmmissing(ind_matching);

% BC linear model
x = DObc.DOconc(ind_matching.ind_sonde);
y = wink.DO_mean(ind_matching.ind_wink);
mdl_bc = fitlm(x,y,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl_bc.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl_bc.Rsquared.Ordinary,2);

fig2 = figure(2);clf
fig2.WindowState = 'maximized';
h = plot(mdl_bc,'marker','.');
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([0 ylim], [0 ylim],'--k')
errorbar(x,y,wink.DO_std(ind_matching.ind_wink),'.b','MarkerSize',dotsize,'LineWidth',2)
legend('',[eqn,newline,'R^2 = ',R2],'1:1 line','Location','southeast')
xlabel('Re-calculated DO Conc (\mumol/L)','interpreter','tex')          
ylabel('Winkler DO Conc (\mumol/L)','interpreter','tex')  
title([site,' BC DO Validation'])
pbaspect([1 1 1]); % Make relative lengths of axes equal

%====ERDC linear regression================================================
DOerdc = table(DOcorr.erdc.datetime_utc,DOcorr.erdc.DOconc);
DOerdc.Properties.VariableNames = {'datetime_utc','DOconc'};
DOerdc = rmmissing(DOerdc);

% Find closest matching indices
t1 = wink.datetime_utc(ind_wink);
t2 = DOerdc.datetime_utc;
ind_sonde = interp1(t2,1:length(t2),t1,'nearest','extrap');

% If closest matching samples are more than 25 minutes off, do not include
for i = 1:length(ind_sonde(~isnan(ind_sonde)))
    dt = between(DOerdc.datetime_utc(ind_sonde(i)),wink.datetime_utc(ind_wink(i)));
    if abs(time(dt)) > minutes(25)
        ind_sonde(i) = NaN;
    end
end

ind_matching = table(ind_wink,ind_sonde);
ind_matching = rmmissing(ind_matching);

% ERDC linear model
x = DOerdc.DOconc(ind_matching.ind_sonde);
y = wink.DO_mean(ind_matching.ind_wink);
mdl_erdc = fitlm(x,y,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl_erdc.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl_erdc.Rsquared.Ordinary,2);

fig3 = figure(3);clf
fig3.WindowState = 'maximized';
h = plot(mdl_erdc,'marker','.');
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([0 ylim], [0 ylim],'--k')
errorbar(x,y,wink.DO_std(ind_matching.ind_wink),'.b','MarkerSize',dotsize,'LineWidth',2)
legend('',[eqn,newline,'R^2 = ',R2],'1:1 line','Location','southeast')
xlabel('Re-calculated DO Conc (\mumol/L)','interpreter','tex')                                                                                                                                     
ylabel('Winkler DO Conc (\mumol/L)','interpreter','tex')  
title([site,' ERDC DO Validation'])
pbaspect([1 1 1]); % Make relative lengths of axes equal

%====Best-guess linear regression==========================================
DObestguess = table(finalQC.datetime_utc,finalQC.DOconc);
DObestguess.Properties.VariableNames = {'datetime_utc','DOconc'};
DObestguess = rmmissing(DObestguess);

t1 = wink.datetime_utc(ind_wink);
t2 = DObestguess.datetime_utc;
ind_sonde = interp1(t2,1:length(t2),t1,'nearest','extrap');

% If closest matching samples are more than 25 minutes off, do not include
for i = 1:length(ind_sonde(~isnan(ind_sonde)))
    dt = between(DObestguess.datetime_utc(ind_sonde(i)),wink.datetime_utc(ind_wink(i)));
    if abs(time(dt)) > minutes(25)
        ind_sonde(i) = NaN;
    end
end

ind_matching = table(ind_wink,ind_sonde);
ind_matching = rmmissing(ind_matching);

% "Best guess" linear model
x = DObestguess.DOconc(ind_matching.ind_sonde);
y = wink.DO_mean(ind_matching.ind_wink);
mdl = fitlm(x,y,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl.Rsquared.Ordinary,2);

fig4 = figure(4);clf
fig4.WindowState = 'maximized';
h = plot(mdl,'marker','.');
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([0 ylim], [0 ylim],'--k')
errorbar(x,y,wink.DO_std(ind_matching.ind_wink),'.b','MarkerSize',dotsize,'LineWidth',2)
legend('Data',[eqn,newline,'R^2 = ',R2],'1:1 line','Location','southeast')
xlabel('"Best-Guess" DO Conc (\mumol/L)','interpreter','tex')                                                                                                                                     
ylabel('Winkler DO Conc (\mumol/L)','interpreter','tex')    
title([site,' Best Guess DO Validation'])
pbaspect([1 1 1]); % Make relative lengths of axes equal

%====Save regression model results=========================================
option = questdlg('Save regression model?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\',site,'\cleaned\validation'])
        save([site,'-winkler_mdl.mat'],'mdl')
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end

%====Save the plots========================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\open-water-platform\',site,'\validation\winklers'])
        saveas(fig2,[site,'_bc-linreg.png'])
        saveas(fig2,[site,'_bc-linreg.fig'])
        saveas(fig3,[site,'_erdc-linreg.png'])
        saveas(fig3,[site,'_erdc-linreg.fig'])  
        saveas(fig4,[site,'_bestguess-linreg.png'])
        saveas(fig4,[site,'_bestguess-linreg.png'])   
        disp('Plots saved!')
    case 'No'
        disp('Plots not saved.')
end

%% Gain correction to South Deployment 15 data
% Calculate gain corrections from beginning and end Winklers
% t0 = t2(ind_sonde(8));
% tn = t2(ind_sonde(9));
% gain_t0 = y(5)/x(5);
% gain_tn = y(6)/x(6);
% 
% tvec = t2(ind_sonde(8):ind_sonde(9));
% dt = seconds(tvec - t0);
% 
% slope = (gain_tn - gain_t0) / seconds(tn-t0);  % [umol L-1 s-1]
% 
% gain_ti = gain_t0 + slope*dt; % [umol L-1]
% 
% DOadj = DObestguess.DOconc(ind_sonde(8):ind_sonde(9)) .* gain_ti;
% 
% fig5 = figure(5);clf
% fig5.WindowState = 'maximized';
% plot(finalQC.datetime_utc,finalQC.DOconc,'.k','MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','"Best Guess"')
% hold on
% plot(DObestguess.datetime_utc(ind_sonde(8):ind_sonde(9)),DOadj,'.r','MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','Corrected DO')
% errorbar(wink.datetime_utc(ind_wink),wink.DO_mean(ind_wink),wink.DO_std(ind_wink),'.','Color',rgb('gold'),'MarkerSize',12,'LineWidth',2,'DisplayName','Winkler')
% xline(finalQC.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
% xlim([dt1 dt2])                 % Use same x limits
% % ylim([0 450])
% ylabel('DO Conc (\mumol/L)')
% title([site,' - "Best Guess" DO Concentration'])
% set(gca,'FontSize',fontsize)
% legend('show','Location','best')
