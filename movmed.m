%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% movmed.m
% This script applies the moving median test to the cleaned sonde data
% (from initialQC.m).
%
% Code that requires manual input is commented with "INPUTS".
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 
% Last updated: 9/30/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

prompt = {'Chose which site to QC'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

cd([rootpath,'open-water-platform-data\',site,'\cleaned\initial-qc'])

prompt = {'Choose which sonde to QC'};
sonde = questdlg(prompt,'Sonde Selection','BC','ERDC','Cancel','BC');

switch sonde
    case 'BC'
        load([site,'-bc-cleaned.mat'])
        dat = sonde1_cleaned;
    case 'ERDC'
        load([site,'-erdc-cleaned.mat'])
        dat = sonde2_cleaned;
end

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = dat.datetime_utc(end);
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 12;
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
        switch sonde
            case 'BC'
                label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
                    'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
                    'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
                    'Deployment 14','Deployment 16','Deployment 17','Deployment 18'};
            case 'ERDC'
                label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 7',...
                    'Deployment 8','Deployment 9','Deployment 10','Deployment 11',...
                    'Deployment 12','Deployment 13','Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
        end
end

% Find indices of deployment changes
ind_dep = find(diff(dat.deployment) > 0);

% Double Moving Median Test
% window = 144; % 24 h
window = 150;   % 25 h
C = 1; % Consistency constant (should be 1.4826 for normally distributed values)
k = 3;

%====DEPTH=================================================================
% Find numbers below and above the moving median
mmed = movmedian(dat.depth,window,"omitmissing");
ind_low = find(dat.depth <= mmed);
ind_high = find(dat.depth >= mmed);
% Calculate moving MAD for values below the median
mmed_low = movmedian(dat.depth(ind_low),window,"omitmissing");
mmad_low = C*movmad(dat.depth(ind_low),window,"omitmissing");
% Calculate moving MAD for values above the median
mmed_high = movmedian(dat.depth(ind_high),window,"omitmissing");
mmad_high = C*movmad(dat.depth(ind_high),window,"omitmissing");

% Lower outlier threshold
lower = mmed_low - k*mmad_low;
ind_bad_low  = ind_low(dat.depth(ind_low) < lower);
% Upper outlier threshold
upper = mmed_high + k*mmad_high;
ind_bad_high = ind_high(dat.depth(ind_high) > upper);

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(dat.datetime_utc,dat.depth,'.k','DisplayName','Initially Cleaned Data')
hold on
plot(dat.datetime_utc,mmed,'.','DisplayName','Moving Median')
% plot(dat.datetime_utc(ind_low),lower,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Lower threshold')
% plot(dat.datetime_utc(ind_high),upper,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Upper threshold')
plot(dat.datetime_utc(ind_bad_low),dat.depth(ind_bad_low)','o','MarkerSize',circlesize,'DisplayName','Exceeds lower threshold')
plot(dat.datetime_utc(ind_bad_high),dat.depth(ind_bad_high)','o','MarkerSize',circlesize,'DisplayName','Exceeds upper threshold')
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend('show','location','best')
ylabel('Depth (m)')
title([site,' ',sonde,' - Flagged Points from Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites

% Compare distributions of data before/after
ind_bad = sort([ind_bad_low;ind_bad_high]);
depth_new = dat.depth;
depth_new(ind_bad) = [];

fig2 = figure(2);clf
fig2.WindowState = 'maximized';
histogram(dat.depth,'DisplayName','Before Moving Median Test')
hold on
histogram(depth_new,'DisplayName','After Moving Median Test')
legend('show')
xlabel('Depth (m)')
ylabel('Frequency')
title([site,' ',sonde,' - Data Histogram'])

depth_flags.ind_movmed = ind_bad;

%====TEMPERATURE===========================================================
% Find numbers below and above the moving median
mmed = movmedian(dat.temperature,window,"omitmissing");
ind_low = find(dat.temperature <= mmed);
ind_high = find(dat.temperature >= mmed);
% Calculate moving MAD for values below the median
mmed_low = movmedian(dat.temperature(ind_low),window,"omitmissing");
mmad_low = C*movmad(dat.temperature(ind_low),window,"omitmissing");
% Calculate moving MAD for values above the median
mmed_high = movmedian(dat.temperature(ind_high),window,"omitmissing");
mmad_high = C*movmad(dat.temperature(ind_high),window,"omitmissing");

% Lower outlier threshold
lower = mmed_low - k*mmad_low;
ind_bad_low  = ind_low(dat.temperature(ind_low) < lower);
% Upper outlier threshold
upper = mmed_high + k*mmad_high;
ind_bad_high = ind_high(dat.temperature(ind_high) > upper);

fig3 = figure(3);clf
fig3.WindowState = 'maximized';
plot(dat.datetime_utc,dat.temperature,'.k','DisplayName','Initially Cleaned Data')
hold on
plot(dat.datetime_utc,mmed,'.','DisplayName','Moving Median')
% plot(dat.datetime_utc(ind_low),lower,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Lower threshold')
% plot(dat.datetime_utc(ind_high),upper,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Upper threshold')
plot(dat.datetime_utc(ind_bad_low),dat.temperature(ind_bad_low)','o','MarkerSize',circlesize,'DisplayName','Exceeds lower threshold')
plot(dat.datetime_utc(ind_bad_high),dat.temperature(ind_bad_high)','o','MarkerSize',circlesize,'DisplayName','Exceeds upper threshold')
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend('show','location','best')
ylabel('Temperature (^oC)')
title([site,' ',sonde,' - Flagged Points from Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites

% Compare distributions of data before/after
ind_bad = sort([ind_bad_low;ind_bad_high]);
T_new = dat.temperature;
T_new(ind_bad) = [];

fig4 = figure(4);clf
fig4.WindowState = 'maximized';
histogram(dat.temperature,'DisplayName','Before Moving Median Test')
hold on
histogram(T_new,'DisplayName','After Moving Median Test')
legend('show')
xlabel('Temperature (^oC)')
ylabel('Frequency')
title([site,' ',sonde,' - Data Histogram'])

T_flags.ind_movmed = ind_bad;

%====DO CONCENTRATION======================================================
% Find numbers below and above the moving median
mmed = movmedian(dat.DO_conc,window,"omitmissing");
ind_low = find(dat.DO_conc <= mmed);
ind_high = find(dat.DO_conc >= mmed);
% Calculate moving MAD for values below the median
mmed_low = movmedian(dat.DO_conc(ind_low),window,"omitmissing");
mmad_low = C*movmad(dat.DO_conc(ind_low),window,"omitmissing");
% Calculate moving MAD for values above the median
mmed_high = movmedian(dat.DO_conc(ind_high),window,"omitmissing");
mmad_high = C*movmad(dat.DO_conc(ind_high),window,"omitmissing");

% Lower outlier threshold
lower = mmed_low - k*mmad_low;
ind_bad_low  = ind_low(dat.DO_conc(ind_low) < lower);
% Upper outlier threshold
upper = mmed_high + k*mmad_high;
ind_bad_high = ind_high(dat.DO_conc(ind_high) > upper);

fig5 = figure(5);clf
fig5.WindowState = 'maximized';
plot(dat.datetime_utc,dat.DO_conc,'.k','DisplayName','Initially Cleaned Data')
hold on
plot(dat.datetime_utc,mmed,'.','DisplayName','Moving Median')
plot(dat.datetime_utc(ind_bad_low),dat.DO_conc(ind_bad_low)','o','MarkerSize',circlesize,'DisplayName','Exceeds lower threshold')
plot(dat.datetime_utc(ind_bad_high),dat.DO_conc(ind_bad_high)','o','MarkerSize',circlesize,'DisplayName','Exceeds upper threshold')
% plot(dat.datetime_utc(ind_low),lower,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Lower threshold')
% plot(dat.datetime_utc(ind_high),upper,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Upper threshold')
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
hold off
xlabel('UTC')
ylabel('DO Concentration (\mumol/L)')
legend('show','location','best')
title([site,' ',sonde,' - Flagged Points from Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites

% Compare distributions of data before/after
ind_bad = sort([ind_bad_low;ind_bad_high]);
DO_new = dat.DO_conc;
DO_new(ind_bad) = [];

fig6 = figure(6);clf
fig6.WindowState = 'maximized';
histogram(dat.DO_conc,'DisplayName','Before Moving Median Test')
hold on
histogram(DO_new,'DisplayName','After Moving Median Test')
legend('show')
xlabel('DO Concentration (\mumol/L)')
ylabel('Frequency')
title([site,' ',sonde,' - Data Histogram'])

DO_flags.ind_movmed = ind_bad;

%====SALINITY==============================================================
% Find numbers below and above the moving median
mmed = movmedian(dat.salinity,window,"omitmissing");
ind_low = find(dat.salinity <= mmed);
ind_high = find(dat.salinity >= mmed);
% Calculate moving MAD for values below the median
mmed_low = movmedian(dat.salinity(ind_low),window,"omitmissing");
mmad_low = C*movmad(dat.salinity(ind_low),window,"omitmissing");
% Calculate moving MAD for values above the median
mmed_high = movmedian(dat.salinity(ind_high),window,"omitmissing");
mmad_high = C*movmad(dat.salinity(ind_high),window,"omitmissing");

% Lower outlier threshold
lower = mmed_low - k*mmad_low;
ind_bad_low  = ind_low(dat.salinity(ind_low) < lower);
% Upper outlier threshold
upper = mmed_high + k*mmad_high;
ind_bad_high = ind_high(dat.salinity(ind_high) > upper);

fig7 = figure(7);clf
fig7.WindowState = 'maximized';
plot(dat.datetime_utc,dat.salinity,'.k','DisplayName','Initially Cleaned Data')
hold on
plot(dat.datetime_utc,mmed,'.','DisplayName','Moving Median')
plot(dat.datetime_utc(ind_bad_low),dat.salinity(ind_bad_low)','o','MarkerSize',circlesize,'DisplayName','Exceeds lower threshold')
plot(dat.datetime_utc(ind_bad_high),dat.salinity(ind_bad_high)','o','MarkerSize',circlesize,'DisplayName','Exceeds upper threshold')
% plot(dat.datetime_utc(ind_low),lower,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Lower threshold')
% plot(dat.datetime_utc(ind_high),upper,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Upper threshold')
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend('show','location','best')
ylabel('Salinity (psu)')
title([site,' ',sonde,' - Flagged Points from Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites

% Compare distributions of data before/after
ind_bad = sort([ind_bad_low;ind_bad_high]);
S_new = dat.salinity;
S_new(ind_bad) = [];

fig8 = figure(8);clf
fig8.WindowState = 'maximized';
histogram(dat.salinity,'DisplayName','Before Moving Median Test')
hold on
histogram(S_new,'DisplayName','After Moving Median Test')
legend('show')
xlabel('Salinity (psu)')
ylabel('Frequency')
title([site,' ',sonde,' - Data Histogram'])

S_flags.ind_movmed = ind_bad;

%====pH====================================================================
% Find numbers below and above the moving median
mmed = movmedian(dat.pH,window,"omitmissing");
ind_low = find(dat.pH <= mmed);
ind_high = find(dat.pH >= mmed);
% Calculate moving MAD for values below the median
mmed_low = movmedian(dat.pH(ind_low),window,"omitmissing");
mmad_low = C*movmad(dat.pH(ind_low),window,"omitmissing");
% Calculate moving MAD for values above the median
mmed_high = movmedian(dat.pH(ind_high),window,"omitmissing");
mmad_high = C*movmad(dat.pH(ind_high),window,"omitmissing");

% Lower outlier threshold
lower = mmed_low - k*mmad_low;
ind_bad_low  = ind_low(dat.pH(ind_low) < lower);
% Upper outlier threshold
upper = mmed_high + k*mmad_high;
ind_bad_high = ind_high(dat.pH(ind_high) > upper);

fig9 = figure(9);clf
fig9.WindowState = 'maximized';
plot(dat.datetime_utc,dat.pH,'.k','DisplayName','Initially Cleaned Data')
hold on
plot(dat.datetime_utc,mmed,'.','DisplayName','Moving Median')
plot(dat.datetime_utc(ind_bad_low),dat.pH(ind_bad_low)','o','MarkerSize',circlesize,'DisplayName','Exceeds lower threshold')
plot(dat.datetime_utc(ind_bad_high),dat.pH(ind_bad_high)','o','MarkerSize',circlesize,'DisplayName','Exceeds upper threshold')
% plot(dat.datetime_utc(ind_low),lower,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Lower threshold')
% plot(dat.datetime_utc(ind_high),upper,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Upper threshold')
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend('show','location','best')
ylabel('pH')
title([site,' ',sonde,' - Flagged Points from Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites

% Compare distributions of data before/after
ind_bad = sort([ind_bad_low;ind_bad_high]);
pH_new = dat.pH;
pH_new(ind_bad) = [];

fig10 = figure(10);clf
fig10.WindowState = 'maximized';
histogram(dat.pH,'DisplayName','Before Moving Median Test')
hold on
histogram(pH_new,'DisplayName','After Moving Median Test')
legend('show')
xlabel('pH')
ylabel('Frequency')
title([site,' ',sonde,' - Data Histogram'])

pH_flags.ind_movmed = ind_bad;

switch sonde
    case 'BC'
        %====Chl a=================================================================
        % Find numbers below and above the moving median
        mmed = movmedian(dat.chla,window,"omitmissing");
        ind_low = find(dat.chla <= mmed);
        ind_high = find(dat.chla >= mmed);
        % Calculate moving MAD for values below the median
        mmed_low = movmedian(dat.chla(ind_low),window,"omitmissing");
        mmad_low = C*movmad(dat.chla(ind_low),window,"omitmissing");
        % Calculate moving MAD for values above the median
        mmed_high = movmedian(dat.chla(ind_high),window,"omitmissing");
        mmad_high = C*movmad(dat.chla(ind_high),window,"omitmissing");

        % Lower outlier threshold
        lower = mmed_low - k*mmad_low;
        ind_bad_low  = ind_low(dat.chla(ind_low) < lower);
        % Upper outlier threshold
        upper = mmed_high + k*mmad_high;
        ind_bad_high = ind_high(dat.chla(ind_high) > upper);

        fig11 = figure(11);clf
        fig11.WindowState = 'maximized';
        plot(dat.datetime_utc,dat.chla,'.k','DisplayName','Initially Cleaned Data')
        hold on
        plot(dat.datetime_utc,mmed,'.','DisplayName','Moving Median')
        plot(dat.datetime_utc(ind_bad_low),dat.chla(ind_bad_low)','o','MarkerSize',circlesize,'DisplayName','Exceeds lower threshold')
        plot(dat.datetime_utc(ind_bad_high),dat.chla(ind_bad_high)','o','MarkerSize',circlesize,'DisplayName','Exceeds upper threshold')
        % plot(dat.datetime_utc(ind_low),lower,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Lower threshold')
        % plot(dat.datetime_utc(ind_high),upper,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Upper threshold')
        xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
        legend('show','location','best')
        ylabel('Chl a (RFU)')
        title([site,' ',sonde,' - Flagged Points from Moving Median Test'])
        xlim([dt1 dt2])                 % Use same x limits for comparing sites

        % Compare distributions of data before/after
        ind_bad = sort([ind_bad_low;ind_bad_high]);
        chla_new = dat.chla;
        chla_new(ind_bad) = [];

        fig12 = figure(12);clf
        fig12.WindowState = 'maximized';
        histogram(dat.chla,'DisplayName','Before Moving Median Test')
        hold on
        histogram(chla_new,'DisplayName','After Moving Median Test')
        legend('show')
        xlabel('Chl a (RFU)')
        ylabel('Frequency')
        title([site,' ',sonde,' - Data Histogram'])

        chla_flags.ind_movmed = ind_bad;

    case 'ERDC'
        %====Turbidity=====================================================
        % Find numbers below and above the moving median
        mmed = movmedian(dat.turbidity,window,"omitmissing");
        ind_low = find(dat.turbidity <= mmed);
        ind_high = find(dat.turbidity >= mmed);
        % Calculate moving MAD for values below the median
        mmed_low = movmedian(dat.turbidity(ind_low),window,"omitmissing");
        mmad_low = C*movmad(dat.turbidity(ind_low),window,"omitmissing");
        % Calculate moving MAD for values above the median
        mmed_high = movmedian(dat.turbidity(ind_high),window,"omitmissing");
        mmad_high = C*movmad(dat.turbidity(ind_high),window,"omitmissing");

        % Lower outlier threshold
        lower = mmed_low - k*mmad_low;
        ind_bad_low  = ind_low(dat.turbidity(ind_low) < lower);
        % Upper outlier threshold
        upper = mmed_high + k*mmad_high;
        ind_bad_high = ind_high(dat.turbidity(ind_high) > upper);

        fig11 = figure(13);clf
        fig11.WindowState = 'maximized';
        plot(dat.datetime_utc,dat.turbidity,'.k','DisplayName','Initially Cleaned Data')
        hold on
        plot(dat.datetime_utc,mmed,'.','DisplayName','Moving Median')
        plot(dat.datetime_utc(ind_bad_low),dat.turbidity(ind_bad_low)','o','MarkerSize',circlesize,'DisplayName','Exceeds lower threshold')
        plot(dat.datetime_utc(ind_bad_high),dat.turbidity(ind_bad_high)','o','MarkerSize',circlesize,'DisplayName','Exceeds upper threshold')
        % plot(dat.datetime_utc(ind_low),lower,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Lower threshold')
        % plot(dat.datetime_utc(ind_high),upper,'.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','Upper threshold')
        xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
        legend('show','location','best')
        ylabel('Turbidity (NTU)')
        title([site,' ',sonde,' - Flagged Points from Moving Median Test'])
        xlim([dt1 dt2])                 % Use same x limits for comparing sites

        % Compare distributions of data before/after
        ind_bad = sort([ind_bad_low;ind_bad_high]);
        turbidity_new = dat.turbidity;
        turbidity_new(ind_bad) = [];

        fig12 = figure(12);clf
        fig12.WindowState = 'maximized';
        histogram(dat.turbidity,'DisplayName','Before Moving Median Test')
        hold on
        histogram(turbidity_new,'DisplayName','After Moving Median Test')
        legend('show')
        xlabel('Turbidity (NTU)')
        ylabel('Frequency')
        title([site,' ',sonde,' - Data Histogram'])

        turbidity_flags.ind_movmed = ind_bad;
end

% figure,clf
% % plot(dat.datetime_utc,dat.turbidity,'.k','DisplayName','Initially Cleaned Data')
% % hold on
% plot(dat.datetime_utc,mmed,'-','DisplayName','Moving Median (24-h)')
% xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
% legend('show','location','best')
% ylabel('Turbidity (NTU)')
% title([site,' ',sonde])
% xlim([dt1 dt2])                 % Use same x limits for comparing sites


%====Clean data============================================================
option = questdlg('Remove points that failed moving median test (exclude depth and T)?');
switch option
    case 'Yes'
        % dat.depth(depth_flags.ind_movmed) = NaN;  % Movmed test looks like it cuts off tops & bottoms of diel cycles
        % dat.temperature(T_flags.ind_movmed) = NaN;  % Movmed test looks like it cuts off tops & bottoms of diel cycles
        dat.DO_conc(DO_flags.ind_movmed) = NaN;
        dat.salinity(S_flags.ind_movmed) = NaN;
        dat.pH(pH_flags.ind_movmed) = NaN;
        switch sonde
            case 'BC'
                dat.chla(chla_flags.ind_movmed) = NaN;
            case 'ERDC'
                dat.turbidity(turbidity_flags.ind_movmed) = NaN;
        end
        disp('Data cleaned!')
    case 'No'
        disp('Data not cleaned.')
end

%====Plot the cleaned data after movmed test===============================
cd([rootpath,'\figures\open-water-platform\',site,'\data-qc\',sonde,'\movmed'])

% DO Concentration
fig13 = figure(13);clf
fig13.WindowState = 'maximized';
plot(dat.datetime_utc,dat.DO_conc,'.k','MarkerSize',dotsize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
ylabel('DO Concentration (\mumol/L)')
title([site,' ',sonde,' - After Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% Salinity
fig14 = figure(14);clf
fig14.WindowState = 'maximized';
plot(dat.datetime_utc,dat.salinity,'.k','MarkerSize',dotsize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
ylabel('Salinity (psu)')
title([site,' ',sonde,' - After Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
% ylim([-.5 50])
set(gca,'FontSize',fontsize)

% pH
fig15 = figure(15);clf
fig15.WindowState = 'maximized';
plot(dat.datetime_utc,dat.pH,'.k','MarkerSize',dotsize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
ylabel('pH')
title([site,' ',sonde,' - After Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

switch sonde
    case 'BC'
        % Chl a
        fig16 = figure(16);clf
        fig16.WindowState = 'maximized';
        plot(dat.datetime_utc,dat.chla,'.k','MarkerSize',dotsize);
        xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
        ylabel('Chl a (RFU)')
        title([site,' ',sonde,' - After Moving Median Test'])
        xlim([dt1 dt2])                 % Use same x limits for comparing sites
        set(gca,'FontSize',fontsize)

    case 'ERDC'
        % Turbidity
        fig16 = figure(16);clf
        fig16.WindowState = 'maximized';
        plot(dat.datetime_utc,dat.turbidity,'.k','MarkerSize',dotsize);
        xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
        ylabel('Turbidity (NTU)')
        title([site,' ',sonde,' - After Moving Median Test'])
        xlim([dt1 dt2])                 % Use same x limits for comparing sites
        set(gca,'FontSize',fontsize)
end

%====Save results from moving median test==================================
option = questdlg('Save cleaned data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\',site,'\cleaned\movmed'])
        switch sonde
            case 'BC'
                sonde1_cleaned = dat;
                save([site,'-bc-cleaned.mat'],'sonde1_cleaned')
            case 'ERDC'
                sonde2_cleaned = dat;
                save([site,'-erdc-cleaned.mat'],'sonde2_cleaned')
        end
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end

%====Save the flagged and cleaned plots====================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        switch sonde
            case 'BC'
                saveFilePath = ['figures\open-water-platform\',site,'\data-qc\bc\movmed'];
            case 'ERDC'
                saveFilePath = ['figures\open-water-platform\',site,'\data-qc\erdc\movmed'];
        end
        cd([rootpath,saveFilePath])
        saveas(fig1,'depth_movmed.png')
        saveas(fig1,'depth_movmed.fig')
        saveas(fig2,'depth_histogram.png')
        saveas(fig2,'depth_histogram.fig')  
        saveas(fig3,'T_movmed.png')
        saveas(fig3,'T_movmed.fig')
        saveas(fig4,'T_histogram.png')
        saveas(fig4,'T_histogram.fig')
        saveas(fig5,'DO_movmed.png')
        saveas(fig5,'DO_movmed.fig')
        saveas(fig6,'DO_histogram.png')
        saveas(fig6,'DO_histogram.fig')
        saveas(fig7,'S_movmed.png')
        saveas(fig7,'S_movmed.fig')
        saveas(fig8,'S_histogram.png')
        saveas(fig8,'S_histogram.fig')
        saveas(fig9,'pH_movmed.png')
        saveas(fig9,'pH_movmed.fig')
        saveas(fig10,'pH_histogram.png')
        saveas(fig10,'pH_histogram.fig')
        switch sonde
            case 'BC'
                saveas(fig11,'chla_movmed.png')
                saveas(fig11,'chla_movmed.fig')
                saveas(fig12,'chla_histogram.png')
                saveas(fig12,'chla_histogram.fig')
            case 'ERDC'
                saveas(fig11,'turbidity_movmed.png')
                saveas(fig11,'turbidity_movmed.fig')
                saveas(fig12,'turbidity_histogram.png')
                saveas(fig12,'turbidity_histogram.fig')
        end
        saveas(fig13,'DOconc_cleaned.png')
        saveas(fig13,'DOconc_cleaned.fig')
        saveas(fig14,'S_cleaned.png')
        saveas(fig14,'S_cleaned.fig')
        saveas(fig15,'pH_cleaned.png')
        saveas(fig15,'pH_cleaned.fig')
        switch sonde
            case 'BC'
                saveas(fig16,'chla_cleaned.png')
                saveas(fig16,'chla_cleaned.fig')
            case 'ERDC'
                saveas(fig16,'turbidity_cleaned.png')
                saveas(fig16,'turbidity_cleaned.fig')
        end
        
        disp('Plots saved!')
    
    case 'No'
        disp('Plots not saved.')
end


%% Test
% x = [1, 4, 4, 4, 5, 5, 5, 5, 7, 7, 8, 10, 16, 30];
% % x = [4, 10, 15, 18, 19, 20, 501, 502, 503, 504, 3000];
% % C = 1.4826;
% C = 1;
% k = 3;
% 
% med=median(x);
% ind_low = find(x<=med);
% ind_high = find(x>=med);
% mad_low = C*median(abs(x(ind_low) - med));
% mad_high = C*median(abs(x(ind_high) - med));
% lower = med - k*mad_low;
% upper = med + k*mad_high;
% x_low = x(find(x<lower))
% x_high = x(find(x>upper))
% 
% figure,clf
% plot(x,'.')
% hold on
% plot(repmat(lower,length(x),1),'-')
% plot(repmat(upper,length(x),1),'-')