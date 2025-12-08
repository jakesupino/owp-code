%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dupCheck_DO_pH.m
% This script compares the QC'd results from the BC and ERDC sondes for the
% recalculated DO concentration and pH and produces the "best guess" time 
% series data for each of these parameters.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 4/26/2024
% Last updated: 10/24/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

%====Import cleaned BC and ERDC data=======================================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\movmed\'])
load([site,'-bc-cleaned'])
load([site,'-erdc-cleaned'])
dat1 = sonde1_cleaned;
dat2 = sonde2_cleaned;
% Synchronize the BC and ERDC data to a common datetime vector
dat_syn = synchronize(dat1,dat2);

clearvars dat1 dat2 sonde1_cleaned sonde2_cleaned

%====Import best-guess data for depth, S, and T============================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\dupcheck\'])
load([site,'-bestGuess.mat'])

%====Import re-calculated DO concentrations================================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\corrected'])
load([site,'-DOCorr.mat'])

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

%====Find the indices of deployment changes================================
switch site
    case 'South'    % South BC Dep. 15 is missing, so need to manually add in dep number
        ind_d15 = find(dat_syn.deployment_dat2 == 15,1);
        ind_d16 = find(dat_syn.deployment_dat1 == 16,1);
        dat_syn.deployment_dat1(ind_d15:ind_d16) = 15;
end

dep1 = fillmissing(dat_syn.deployment_dat1,'previous'); % Replace NaNs in BC deployment column with previous deployment #
dat_syn.deployment_dat1 = dep1;
ind_dep = find(diff(dat_syn.deployment_dat1) > 0);

switch site
    case 'Gull'
        dep = table([1;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
    case 'North'
        dep = table([2;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
    case 'South'
        dep = table([1;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
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

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = dat_syn.datetime_utc(end);
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 12;
circlesize = 4;

% % Compare DO and pH time series
% fig = figure(1);clf
% fig.WindowState = 'maximized';
% t = tiledlayout(2,1,'TileSpacing','tight');
% 
% ax1 = nexttile;
% plot(DOcorr.bc.datetime_utc,DOcorr.bc.DOconc,'.-','Color',rgb('darkred'),'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','BC - Re-calculated')
% hold on
% plot(DOcorr.erdc.datetime_utc,DOcorr.erdc.DOconc,'.-','Color',rgb('darkblue'),'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','ERDC - Re-calculated')
% xline([bestguess.deployment.datetime_utc(1); bestguess.deployment.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
% ylabel('DO Conc (\mumol/L)')
% legend('show')
% 
% ax2 = nexttile;
% plot(dat_syn.datetime_utc,dat_syn.pH_dat1,'.-','Color',red,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','BC')
% hold on
% plot(dat_syn.datetime_utc,dat_syn.pH_dat2,'.-','Color',blue,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','ERDC')
% xline([bestguess.deployment.datetime_utc(1); bestguess.deployment.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
% ylabel('pH')
% legend('show')
% 
% linkaxes([ax1 ax2],'x')

%====Create "best-guess" DO conc===========================================
% Find mean of DO time series 
DOcorr_syn = synchronize(DOcorr.bc,DOcorr.erdc);
DOcorr_syn.Properties.VariableNames = {'bc','erdc'};
mean_DOcorr = mean(DOcorr_syn(:,1:2),2);
DOcorr_syn.mean = mean_DOcorr.mean;
dt_utc = DOcorr_syn.datetime_utc;

% Find indices of deployment changes in synchronized, corrected DO data
t1 = dat_syn.datetime_utc(dep.ind);
t2 = DOcorr_syn.datetime_utc;
ind_depDO = interp1(t2,1:length(t2),t1,'nearest');
ind_depDO = [1;ind_depDO(2:end)];

cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\DOconc'])
% Check agreement between BC and ERDC sondes
DO_bc = DOcorr_syn.bc;
DO_erdc = DOcorr_syn.erdc;
% threshold = 20; % umol/L
% ind_diverge = find(abs(DO_bc - DO_erdc) > threshold);
% pc_exceed = height(ind_diverge)/height(DOcorr_syn)*100; % Percentage of synchronized data that exceed threshold
% txt = ['Differs by >',num2str(threshold),' \mumol/L (',num2str(pc_exceed,'%.1f'),'%)'];

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(dt_utc,DOcorr_syn.mean,'-','Color',rgb('magenta'),'DisplayName','Mean Value')
hold on
plot(dt_utc,DO_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC (recalculated)');
plot(dt_utc,DO_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC (recalculated)')
% plot(dt_utc(ind_diverge),DO_bc(ind_diverge),'o','color',rgb('goldenrod'),'MarkerSize',circlesize,'DisplayName',txt)
errorbar(wink.datetime_utc(ind_wink),wink.DO_mean(ind_wink),wink.DO_std(ind_wink),'o','MarkerSize',6,'LineWidth',2,'color',rgb('goldenrod'),'DisplayName','Winkler Sample')
xline(dt_utc(ind_depDO),'--',label,'HandleVisibility','off')
ylabel('DO concentration (\mumol/L)')
legend('show','location','best')
title([site,' - After Initial QC and Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% Choose "better" deployment
switch site
    case 'Gull'
        d1 = DOcorr_syn.erdc(ind_depDO(1):ind_depDO(2));       % Dep 1
        d2 = DOcorr_syn.mean(ind_depDO(2)+1:ind_depDO(3));     % Dep 2
        d5 = DOcorr_syn.mean(ind_depDO(3)+1:ind_depDO(4));     % Dep 5
        d6 = DOcorr_syn.bc(ind_depDO(4)+1:ind_depDO(5));       % Dep 6
        d7 = DOcorr_syn.bc(ind_depDO(5)+1:ind_depDO(6));       % Dep 7
        d8 = DOcorr_syn.erdc(ind_depDO(6)+1:ind_depDO(7));     % Dep 8
        d9 = DOcorr_syn.mean(ind_depDO(7)+1:ind_depDO(8));     % Dep 9
        d10 = DOcorr_syn.bc(ind_depDO(8)+1:ind_depDO(9));    % Dep 10
        d11 = DOcorr_syn.mean(ind_depDO(9)+1:ind_depDO(10));   % Dep 11
        d12 = DOcorr_syn.mean(ind_depDO(10)+1:ind_depDO(11));  % Dep 12
        d13 = DOcorr_syn.erdc(ind_depDO(11)+1:ind_depDO(12));  % Dep 13
        d14 = DOcorr_syn.erdc(ind_depDO(12)+1:ind_depDO(13));  % Dep 14
        d15 = DOcorr_syn.mean(ind_depDO(13)+1:ind_depDO(14));  % Dep 15
        d16 = DOcorr_syn.mean(ind_depDO(14)+1:ind_depDO(15));  % Dep 16
        d17 = DOcorr_syn.bc(ind_depDO(15)+1:ind_depDO(16));  % Dep 17
        d18 = DOcorr_syn.bc(ind_depDO(16)+1:end);            % Dep 18

        DO_bestguess = [d1;d2;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];
        
        % Calculate mean absolute difference for Deployments 11 & 12
        diff_DOconc = diff([DOcorr_syn.bc(ind_depDO(9):ind_depDO(11)),DOcorr_syn.erdc(ind_depDO(9):ind_depDO(11))],1,2);
        two_sigma.DOconc = mean(abs(diff_DOconc),'omitmissing');

    case 'North'
        d2 = DOcorr_syn.erdc(ind_depDO(1):ind_depDO(2));       % Dep 2
        d6 = DOcorr_syn.bc(ind_depDO(2)+1:ind_depDO(3));       % Dep 6
        d7 = DOcorr_syn.mean(ind_depDO(3)+1:ind_depDO(4));     % Dep 7
        d8 = DOcorr_syn.mean(ind_depDO(4)+1:ind_depDO(5));     % Dep 8
        d9 = DOcorr_syn.mean(ind_depDO(5)+1:ind_depDO(6));     % Dep 9
        d10 = DOcorr_syn.mean(ind_depDO(6)+1:ind_depDO(7));    % Dep 10
        d11 = DOcorr_syn.bc(ind_depDO(7)+1:ind_depDO(8));      % Dep 11
        d12 = DOcorr_syn.mean(ind_depDO(8)+1:ind_depDO(9));    % Dep 12
        d13 = DOcorr_syn.erdc(ind_depDO(9)+1:ind_depDO(10));   % Dep 13
        d14 = DOcorr_syn.erdc(ind_depDO(10)+1:ind_depDO(11));  % Dep 14
        d15 = DOcorr_syn.mean(ind_depDO(11)+1:ind_depDO(12));  % Dep 15
        d16 = DOcorr_syn.erdc(ind_depDO(12)+1:ind_depDO(13));  % Dep 16
        d17 = DOcorr_syn.erdc(ind_depDO(13)+1:end);            % Dep 17

        DO_bestguess = [d2;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17];
        
        % Calculate mean absolute difference for Deployments 10 & 11
        diff_DOconc = diff([DOcorr_syn.bc(ind_depDO(6):ind_depDO(8)),DOcorr_syn.erdc(ind_depDO(6):ind_depDO(8))],1,2);
        two_sigma.DOconc = mean(abs(diff_DOconc),'omitmissing');

    case 'South'
        d1 = DOcorr_syn.mean(ind_depDO(1):ind_depDO(2));       % Dep 1
        d2 = DOcorr_syn.mean(ind_depDO(2)+1:ind_depDO(3));     % Dep 2
        d4 = DOcorr_syn.bc(ind_depDO(3)+1:ind_depDO(4));       % Dep 4
        d5 = DOcorr_syn.bc(ind_depDO(4)+1:ind_depDO(5));       % Dep 5
        d6 = DOcorr_syn.bc(ind_depDO(5)+1:ind_depDO(6));       % Dep 6
        d7 = DOcorr_syn.erdc(ind_depDO(6)+1:ind_depDO(7));     % Dep 7
        d8 = DOcorr_syn.erdc(ind_depDO(7)+1:ind_depDO(8));     % Dep 8
        d9 = DOcorr_syn.erdc(ind_depDO(8)+1:ind_depDO(9));     % Dep 9
        d10 = DOcorr_syn.bc(ind_depDO(9)+1:ind_depDO(10));     % Dep 10
        d11 = DOcorr_syn.erdc(ind_depDO(10)+1:ind_depDO(11));  % Dep 11
        d12 = DOcorr_syn.erdc(ind_depDO(11)+1:ind_depDO(12));  % Dep 12
        d13 = DOcorr_syn.erdc(ind_depDO(12)+1:ind_depDO(13));  % Dep 13
        d14 = DOcorr_syn.erdc(ind_depDO(13)+1:ind_depDO(14));  % Dep 14
        d15 = NaN(abs(ind_depDO(14) - ind_depDO(15)),1);       % Dep 15
        d16 = DOcorr_syn.mean(ind_depDO(15)+1:ind_depDO(16));  % Dep 16
        d17 = DOcorr_syn.mean(ind_depDO(16)+1:ind_depDO(17));  % Dep 17
        d18 = DOcorr_syn.erdc(ind_depDO(17)+1:end);            % Dep 18

        DO_bestguess = [d1;d2;d4;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];

        % Remove "obvious" outliers
        DO_bestguess(28535:ind_depDO(6)) = NaN;
        
        % Calculate mean absolute difference for Deployments 10 & 11
        diff_DOconc = diff([DOcorr_syn.bc(ind_depDO(9):ind_depDO(11)),DOcorr_syn.erdc(ind_depDO(9):ind_depDO(11))],1,2);
        two_sigma.DOconc = mean(abs(diff_DOconc),'omitmissing');
end

clearvars d1 d2 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 d17

DO_bestguess = table(DOcorr_syn.datetime_utc,DO_bestguess);
DO_bestguess.Properties.VariableNames = {'datetime_utc','DOconc'};
DO_bestguess = table2timetable(DO_bestguess);
DO_bestguess = rmmissing(DO_bestguess);

cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\DOconc'])
fig2 = figure(2);clf
fig2.WindowState = 'maximized';
plot(DOcorr_syn.datetime_utc,DOcorr_syn.bc,'.','Color',red,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','BC (recalculated)')
hold on
plot(DOcorr_syn.datetime_utc,DOcorr_syn.erdc,'.','Color',blue,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','ERDC (recalculated)')
plot(DOcorr_syn.datetime_utc,DOcorr_syn.mean,'.','Color',rgb('magenta'),'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','Mean')
plot(DO_bestguess.datetime_utc,DO_bestguess.DOconc,'.k','MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','"Best Guess"')
errorbar(wink.datetime_utc(ind_wink),wink.DO_mean(ind_wink),wink.DO_std(ind_wink),'o','MarkerSize',6,'LineWidth',2,'color',rgb('goldenrod'),'DisplayName','Winkler Sample')
xline(dt_utc(ind_depDO),'--',label,'HandleVisibility','off')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylim([0 450])
ylabel('DO Conc (\mumol/L)')
title([site,' - "Best-Guess" DO Concentration'])
set(gca,'FontSize',fontsize)
legend('show','location','best')
%%
%====Create "best-guess" pH================================================
% Find mean of pH time series 
bc = table(dat_syn.datetime_utc,dat_syn.pH_dat1);
bc.Properties.VariableNames = {'datetime_utc','pH'};
bc = table2timetable(bc);
bc = rmmissing(bc);

erdc = table(dat_syn.datetime_utc,dat_syn.pH_dat2);
erdc.Properties.VariableNames = {'datetime_utc','pH'};
erdc = table2timetable(erdc);
erdc = rmmissing(erdc);

pH_syn = synchronize(bc,erdc);
pH_syn.Properties.VariableNames = {'bc','erdc'};
mean_pH = mean(pH_syn(:,1:2),2);
pH_syn.mean = mean_pH.mean;
dt_utc = pH_syn.datetime_utc;

% Find indices of deployment changes in synchronized pH data
t1 = dat_syn.datetime_utc(dep.ind);
t2 = pH_syn.datetime_utc;
ind_deppH = interp1(t2,1:length(t2),t1,'nearest');
ind_deppH = [1;ind_deppH(2:end)];

cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\pH'])
% Check agreement between BC and ERDC sondes
pH_bc = pH_syn.bc;
pH_erdc = pH_syn.erdc;
% threshold = 0.3; 
% ind_diverge = find(abs(pH_bc - pH_erdc) > threshold);
% pc_exceed = height(ind_diverge)/height(pH_syn)*100; % Percentage of synchronized data that exceed threshold
% txt = ['Differs by >',num2str(threshold),' (',num2str(pc_exceed,'%.1f'),'%)'];

fig3 = figure(3);clf
fig3.WindowState = 'maximized';
plot(dt_utc,pH_syn.mean,'-','Color',rgb('magenta'),'DisplayName','Mean Value')
hold on
plot(dt_utc,pH_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC');
plot(dt_utc,pH_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
% plot(dt_utc(ind_diverge),pH_bc(ind_diverge),'o','color',rgb('goldenrod'),'MarkerSize',circlesize,'DisplayName',txt)
xline(dt_utc(ind_deppH),'--',label,'HandleVisibility','off')
ylabel('pH')
legend('show','location','best')
title([site,' - After Initial QC and Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% Choose "better" deployment
switch site
    case 'Gull'
        % Dep 17 -- offset adjustment
        delta_d17a = pH_syn.erdc(ind_deppH(15)) - pH_syn.bc(ind_deppH(15)-1);
        delta_d17b = pH_syn.erdc(ind_deppH(16)-1) - pH_syn.bc(95389);
        delta_d17 = mean([delta_d17a delta_d17b]);

        d1 = pH_syn.erdc(ind_deppH(1):ind_deppH(2));       % Dep 1
        d2 = pH_syn.erdc(ind_deppH(2)+1:ind_deppH(3));     % Dep 2
        d5 = pH_syn.bc(ind_deppH(3)+1:ind_deppH(4));       % Dep 5
        d6 = pH_syn.mean(ind_deppH(4)+1:ind_deppH(5));     % Dep 6
        d7 = pH_syn.mean(ind_deppH(5)+1:ind_deppH(6));     % Dep 7
        d8 = pH_syn.erdc(ind_deppH(6)+1:ind_deppH(7));     % Dep 8
        d9 = pH_syn.erdc(ind_deppH(7)+1:ind_deppH(8));     % Dep 9
        d10 = pH_syn.mean(ind_deppH(8)+1:ind_deppH(9));      % Dep 10
        d11 = pH_syn.mean(ind_deppH(9)+1:ind_deppH(10));   % Dep 11
        d12 = pH_syn.mean(ind_deppH(10)+1:ind_deppH(11));  % Dep 12
        d13 = [pH_syn.erdc(ind_deppH(11)+1:69205); pH_syn.bc(69206:ind_deppH(12))];  % Dep 13
        d14 = pH_syn.erdc(ind_deppH(12)+1:ind_deppH(13));  % Dep 14
        d15 = pH_syn.bc(ind_deppH(13)+1:ind_deppH(14));  % Dep 15
        d16 = pH_syn.bc(ind_deppH(14)+1:ind_deppH(15));    % Dep 16
        d17 = pH_syn.erdc(ind_deppH(15)+1:ind_deppH(16)) - delta_d17;  % Dep 17
        d18 = pH_syn.bc(ind_deppH(16)+1:end);            % Dep 18

        pH_bestguess = [d1;d2;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];

        % Calculate mean absolute difference for Deployments 11 & 12
        diff_pH = diff([pH_syn.bc(ind_deppH(9):ind_deppH(11)),pH_syn.erdc(ind_deppH(9):ind_deppH(11))],1,2);
        two_sigma.pH = mean(abs(diff_pH),'omitmissing');

    case 'North'
        d2 = pH_syn.erdc(ind_deppH(1):ind_deppH(2));       % Dep 2
        d6 = pH_syn.bc(ind_deppH(2)+1:ind_deppH(3));       % Dep 6
        d7 = pH_syn.erdc(ind_deppH(3)+1:ind_deppH(4));     % Dep 7
        % d8 = pH_syn.mean(ind_deppH(4)+1:ind_deppH(5));     % Dep 8
        d8 = pH_syn.erdc(ind_deppH(4)+1:ind_deppH(5));     % Dep 8
        % d9 = pH_syn.mean(ind_deppH(5)+1:ind_deppH(6));     % Dep 9
        d9 = pH_syn.erdc(ind_deppH(5)+1:ind_deppH(6));     % Dep 9
        d10 = pH_syn.mean(ind_deppH(6)+1:ind_deppH(7));    % Dep 10
        d11 = pH_syn.bc(ind_deppH(7)+1:ind_deppH(8));      % Dep 11
        d12 = pH_syn.mean(ind_deppH(8)+1:ind_deppH(9));    % Dep 12
        d13 = pH_syn.erdc(ind_deppH(9)+1:ind_deppH(10));   % Dep 13
        d14 = pH_syn.mean(ind_deppH(10)+1:ind_deppH(11));    % Dep 14
        % d15 = pH_syn.mean(ind_deppH(11)+1:ind_deppH(12));  % Dep 15
        d15 = pH_syn.bc(ind_deppH(11)+1:ind_deppH(12));  % Dep 15
        d16 = pH_syn.bc(ind_deppH(12)+1:ind_deppH(13));    % Dep 16
        % d17 = pH_syn.mean(ind_deppH(13)+1:ind_deppH(14);            % Dep 17
        d17 = pH_syn.bc(ind_deppH(13)+1:ind_deppH(14));            % Dep 17
        d18 = [pH_syn.bc(ind_deppH(14)+1:98936); pH_syn.erdc(98937:end)];        % Dep 18
        
        pH_bestguess = [d2;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];
        
        % Remove "obvious" outliers
        pH_bestguess(811) = NaN;
        pH_bestguess(6422:6576) = NaN;

        % Calculate mean absolute difference for Deployments 9 & 10
        diff_pH = diff([pH_syn.bc(ind_deppH(5):ind_deppH(7)),pH_syn.erdc(ind_deppH(5):ind_deppH(7))],1,2);
        two_sigma.pH = mean(abs(diff_pH),'omitmissing');

    case 'South'
        % % Funky deployments
        % % Dep 13, 14, and 15 -- offset adjustment
        % ind1 = ind_deppH(12)-1;
        % ind2 = find(~isnan(pH_syn.erdc(ind_deppH(12):end)),1);
        % ind2 = ind1 + ind2;
        % delta = pH_syn.erdc(ind1-1) - pH_syn.erdc(ind2);

        % Dep 13 -- offset adjustment
        delta_d13a = pH_syn.erdc(58943) - pH_syn.erdc(ind_deppH(12)-1);
        delta_d13b = pH_syn.erdc(ind_deppH(13)-1) - pH_syn.erdc(ind_deppH(13)+1);
        delta_d13 = mean([delta_d13a delta_d13b]);
        
        d1 = pH_syn.erdc(ind_deppH(1):ind_deppH(2));       % Dep 1
        d2 = pH_syn.erdc(ind_deppH(2)+1:ind_deppH(3));     % Dep 2
        d4 = NaN;                                          % Dep 4
        d5 = pH_syn.erdc(ind_deppH(4)+1:ind_deppH(5));     % Dep 5
        d6 = pH_syn.bc(ind_deppH(5)+1:ind_deppH(6));       % Dep 6
        d7 = pH_syn.bc(ind_deppH(6)+1:ind_deppH(7));       % Dep 7
        d8 = pH_syn.erdc(ind_deppH(7)+1:ind_deppH(8));     % Dep 8
        d9 = pH_syn.mean(ind_deppH(8)+1:ind_deppH(9));     % Dep 9
        d10 = pH_syn.bc(ind_deppH(9)+1:ind_deppH(10));     % Dep 10
        d11 = pH_syn.erdc(ind_deppH(10)+1:ind_deppH(11));  % Dep 11
        d12 = pH_syn.erdc(ind_deppH(11)+1:ind_deppH(12));  % Dep 12
        d13 = pH_syn.erdc(ind_deppH(12)+1:ind_deppH(13)) - delta_d13;  % Dep 13
        d14 = pH_syn.erdc(ind_deppH(13)+1:ind_deppH(14));  % Dep 14
        d15 = pH_syn.erdc(ind_deppH(14)+1:ind_deppH(15));  % Dep 15
        d16 = pH_syn.bc(ind_deppH(15)+1:ind_deppH(16));  % Dep 16
        d17 = pH_syn.erdc(ind_deppH(16)+1:ind_deppH(17));  % Dep 17
        d18 = pH_syn.bc(ind_deppH(17)+1:end);            % Dep 18

        pH_bestguess = [d1;d2;d4;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];

        % Calculate mean absolute difference for Deployments 16 & 17
        diff_pH = diff([pH_syn.bc(ind_deppH(15):ind_deppH(17)),pH_syn.erdc(ind_deppH(15):ind_deppH(17))],1,2);
        two_sigma.pH = mean(abs(diff_pH),'omitmissing');
end

clearvars d1 d2 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 d17 d18

pH_bestguess = table(pH_syn.datetime_utc,pH_bestguess);
pH_bestguess.Properties.VariableNames = {'datetime_utc','pH'};
pH_bestguess = table2timetable(pH_bestguess);
pH_bestguess = rmmissing(pH_bestguess);

cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\pH'])
fig4 = figure(4);clf
fig4.WindowState = 'maximized';
plot(pH_syn.datetime_utc,pH_syn.bc,'.','Color',red,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','BC')
hold on
plot(pH_syn.datetime_utc,pH_syn.erdc,'.','Color',blue,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','ERDC')
plot(pH_syn.datetime_utc,pH_syn.mean,'.','Color',rgb('magenta'),'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','Mean')
plot(pH_bestguess.datetime_utc,pH_bestguess.pH,'.k','MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','"Best Guess"')
xline(dt_utc(ind_deppH),'--',label,'HandleVisibility','off')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylabel('pH')
title([site,' - "Best Guess" pH'])
set(gca,'FontSize',fontsize)
legend('show','location','best')

%====Save "best guess" data================================================
% Add "best-guess" DO and pH into structure
bestguess.DOconc = DO_bestguess;
bestguess.pH = pH_bestguess;

option = questdlg('Save best guess data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\',site,'\cleaned\dupcheck'])
        save([site,'-bestGuess.mat'],'bestguess','two_sigma')
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end

%====Save the plots=======================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\DOconc'])
        saveas(fig1,'duplicateComparison.png')
        saveas(fig1,'duplicateComparison.fig')
        saveas(fig2,'bestGuessComparison.png')
        saveas(fig2,'bestGuessComparison.fig')
        cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\pH'])
        saveas(fig3,'duplicateComparison.png')
        saveas(fig3,'duplicateComparison.fig')
        saveas(fig4,'bestGuessComparison.png')
        saveas(fig4,'bestGuessComparison.fig')
        
        disp('Plots saved!')
    
    case 'No'
        disp('Plots not saved.')
end