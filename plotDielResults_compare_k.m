%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotDielResults_compare_k.m
% This script plots the metabolic rates calculated using different k parameterizations
% in dielAnalysis_compare_k.m on daily, monthly, seasonal, and annual time
% scales.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 9/5/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

% Load the sonde data
cd([rootpath,'\diel-method\owp-data\final-qc'])
load([site,'_obs.mat'])   

% Load R detiding results
cd([rootpath,'diel-method\R-results\final-qc\',site])
wtreg_res = readtable('wtreg_res.csv');
wtreg_res.Properties.DimensionNames{1} = 'datetime_utc';
wtreg_res.DateTimeStamp.TimeZone = "UTC";
wtreg_res = table2timetable(wtreg_res);

% Retime weighted regression data to same datetimes as sonde/physical data
wtreg_res.solar_period = [];
newTimes = dat.datetime_utc(1):minutes(10):dat.datetime_utc(end);
wtreg_res_rt = retime(wtreg_res,newTimes,'mean');

% Load R ecometab results
metab_dtd = table2timetable(readtable('metab_dtd.csv'));
metab_dtd.Date.TimeZone = 'UTC';

% Load the MATLAB metabolism results from different parameterizations
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\ro_hunt_theb_const'])
R1a = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\ro_hunt_theb_vary'])
R1b = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\ro_hunt_wann_const'])
R2a = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\ro_hunt_wann_vary'])
R2b = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\wann_const'])
W2a = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\wann_vary'])
W2b = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\emer_const'])
Ea = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\emer_vary'])
Eb = load('diel_res.mat');
             
cd([rootpath,'figures\diel-analysis\sensitivity-analysis'])

% Convert wtreg_res_rt from 10-min intervals to daily intervals
wtreg_res_daily = groupsummary(wtreg_res,"DateTimeStamp","day","mean");
wtreg_res_daily = rmmissing(wtreg_res_daily);
wtreg_res_daily.day_DateTimeStamp = datetime(string(wtreg_res_daily.day_DateTimeStamp),'InputFormat','dd-MMM-yyyy','TimeZone','UTC');

%==========================================================================
%   Plots of metabolic rates using different k parameterizations/implementations
%==========================================================================
% Set bounds for zoomed-in plots
dl = datetime('12-Jun-2023','TimeZone','UTC');
dr = datetime('01-Jul-2023','TimeZone','UTC');

%% ===Compare different implementations of Ro & Hunt 2006==================
% fig = figure;clf
% fig.WindowState = 'maximized';
% 
% yyaxis left
% plot(wtreg_res_daily.day_DateTimeStamp,wtreg_res_daily.mean_WSpd,':','color',rgb('sandybrown'),'HandleVisibility','off')
% ylabel('Wind speed (m/s)')
% ax = gca;
% ax.YColor = rgb('sandybrown');
% 
% yyaxis right
% plot(metab_dtd.Date,metab_dtd.Pg,'b.:','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% hold on
% plot(metab_dtd.Date,metab_dtd.Rt,'k.:','MarkerSize',12,'LineWidth',1,'DisplayName','ecometab - R&H-06')
% plot(metab_dtd.Date,metab_dtd.NEM,'r.:','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% 
% plot(R1a.diel_dtd.date,R1a.diel_dtd.GPP,'b.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% plot(R1a.diel_dtd.date,R1a.diel_dtd.ER,'k.-','MarkerSize',12,'LineWidth',1,'DisplayName','MATLAB - R&H-06 (T-08 implementation w/ const H)')
% plot(R1a.diel_dtd.date,R1a.diel_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% 
% plot(R2b.diel_dtd.date,R1b.diel_dtd.GPP,'b.--','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% plot(R2b.diel_dtd.date,R1b.diel_dtd.ER,'k.--','MarkerSize',12,'LineWidth',1,'DisplayName','MATLAB - R&H-06 (T-08 implementation w/ time-varying H)')
% plot(R2b.diel_dtd.date,R1b.diel_dtd.NEM,'r.--','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% 
% plot(R2a.diel_dtd.date,R2a.diel_dtd.GPP,'b.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% plot(R2a.diel_dtd.date,R2a.diel_dtd.ER,'k.-','MarkerSize',12,'LineWidth',1,'DisplayName','MATLAB - R&H-06 (Sc from W-14; w/ const H)')
% plot(R2a.diel_dtd.date,R2a.diel_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% 
% xlabel('UTC')
% ylabel('Metabolic rate (mmol O_2 m^{-2} d^{-1})','FontSize',14)
% set(gca,'FontSize',14,'LineWidth',2,'YColor','k')
% legend('show','location','southwest')
% title([site,': Comparison of different implementations of Ro & Hunt 2006'])
% xlim([dl dr])
% 
%% ====Compare MATLAB results with different parameterizations (daily time scales)
% fig = figure;clf
% fig.WindowState = 'maximized';
% 
% yyaxis left
% plot(wtreg_res_daily.day_DateTimeStamp,wtreg_res_daily.mean_WSpd,':','color',rgb('sandybrown'),'HandleVisibility','off')
% ylabel('Wind speed (m/s)')
% ax = gca;
% ax.YColor = rgb('sandybrown');
% 
% yyaxis right
% plot(R2a.diel_dtd.date,R2a.diel_dtd.GPP,'b.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% hold on
% plot(R2a.diel_dtd.date,R2a.diel_dtd.ER,'k.-','MarkerSize',12,'LineWidth',1,'DisplayName','R&H-06 (Sc from W-14) - const H')
% plot(R2a.diel_dtd.date,R2a.diel_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% 
% plot(R2b.diel_dtd.date,R2b.diel_dtd.GPP,'b.--','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% plot(R2b.diel_dtd.date,R2b.diel_dtd.ER,'k.--','MarkerSize',12,'LineWidth',1,'DisplayName','R&H-06 (Sc from W-14) - H(t)')
% plot(R2b.diel_dtd.date,R2b.diel_dtd.NEM,'r.--','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
% 
% plot(W2a.diel_dtd.date,W2a.diel_dtd.GPP,'b*-','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')
% plot(W2a.diel_dtd.date,W2a.diel_dtd.ER,'k*-','MarkerSize',8,'LineWidth',1,'DisplayName','W-14 - const H')
% plot(W2a.diel_dtd.date,W2a.diel_dtd.NEM,'r*-','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')
% 
% plot(W2b.diel_dtd.date,W2b.diel_dtd.GPP,'b*--','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')
% plot(W2b.diel_dtd.date,W2b.diel_dtd.ER,'k*--','MarkerSize',8,'LineWidth',1,'DisplayName','W-14 - H(t)')
% plot(W2b.diel_dtd.date,W2b.diel_dtd.NEM,'r*--','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')
% 
% plot(Ea.diel_dtd.date,Ea.diel_dtd.GPP,'b^-','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')
% plot(Ea.diel_dtd.date,Ea.diel_dtd.ER,'k^-','MarkerSize',6,'LineWidth',1,'DisplayName','E-19 - const H')
% plot(Ea.diel_dtd.date,Ea.diel_dtd.NEM,'r^-','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')
% 
% plot(Eb.diel_dtd.date,Eb.diel_dtd.GPP,'b^--','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')
% plot(Eb.diel_dtd.date,Eb.diel_dtd.ER,'k^--','MarkerSize',6,'LineWidth',1,'DisplayName','E-19 - H(t)')
% plot(Eb.diel_dtd.date,Eb.diel_dtd.NEM,'r^--','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')
% 
% xlabel('UTC')
% ylabel('Metabolic rate (mmol O_2 m^{-2} d^{-1})','FontSize',14)
% set(gca,'FontSize',14,'LineWidth',2,'YColor','k')
% legend('show','location','southwest')
% title([site,': Comparison of different k parameterizations (daily time scales)'])
% % xlim([dl dr])
% 
%% ====Compare different parameterizations on MONTHLY time scales==========
% Calculate monthly means 
wtreg_res_monthly = groupsummary(wtreg_res,"DateTimeStamp","month","mean");
R2a_monthly = groupsummary(R2a.diel_dtd,"date","month","mean");
R2b_monthly = groupsummary(R2b.diel_dtd,"date","month","mean");

W2a_monthly = groupsummary(W2a.diel_dtd,"date","month","mean");
W2b_monthly = groupsummary(W2b.diel_dtd,"date","month","mean");

Ea_monthly = groupsummary(Ea.diel_dtd,"date","month","mean");
Eb_monthly = groupsummary(Eb.diel_dtd,"date","month","mean");

fig = figure;clf
fig.WindowState = 'maximized';

plot(R2a_monthly.month_date,R2a_monthly.mean_GPP,'b.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
hold on
plot(R2a_monthly.month_date,R2a_monthly.mean_ER,'k.-','MarkerSize',12,'LineWidth',1,'DisplayName','R&H-06 (Sc from W-14) - const H')
plot(R2a_monthly.month_date,R2a_monthly.mean_NEM,'r.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')

plot(R2b_monthly.month_date,R2b_monthly.mean_GPP,'b.--','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
plot(R2b_monthly.month_date,R2b_monthly.mean_ER,'k.--','MarkerSize',12,'LineWidth',1,'DisplayName','R&H-06 (Sc from W-14) - H(t)')
plot(R2b_monthly.month_date,R2b_monthly.mean_NEM,'r.--','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')

plot(W2a_monthly.month_date,W2a_monthly.mean_GPP,'b*-','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')
plot(W2a_monthly.month_date,W2a_monthly.mean_ER,'k*-','MarkerSize',8,'LineWidth',1,'DisplayName','W-14 - const H')
plot(W2a_monthly.month_date,W2a_monthly.mean_NEM,'r*-','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')

plot(W2b_monthly.month_date,W2b_monthly.mean_GPP,'b*--','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')
plot(W2b_monthly.month_date,W2b_monthly.mean_ER,'k*--','MarkerSize',8,'LineWidth',1,'DisplayName','W-14 - H(t)')
plot(W2b_monthly.month_date,W2b_monthly.mean_NEM,'r*--','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')

plot(Ea_monthly.month_date,Ea_monthly.mean_GPP,'b^-','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')
plot(Ea_monthly.month_date,Ea_monthly.mean_ER,'k^-','MarkerSize',6,'LineWidth',1,'DisplayName','E-19 - const H')
plot(Ea_monthly.month_date,Ea_monthly.mean_NEM,'r^-','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')

plot(Eb_monthly.month_date,Eb_monthly.mean_GPP,'b^--','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')
plot(Eb_monthly.month_date,Eb_monthly.mean_ER,'k^--','MarkerSize',6,'LineWidth',1,'DisplayName','E-19 - H(t)')
plot(Eb_monthly.month_date,Eb_monthly.mean_NEM,'r^--','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')

xlabel('UTC')
ylabel('Metabolic rate (mmol O_2 m^{-2} d^{-1})','FontSize',14)
set(gca,'FontSize',14,'LineWidth',2,'YColor','k')
legend('show','location','southwest')
title([site,': Comparison of different k parameterizations (monthly time scales)'])

%% ====Compare different parameterizations on SEASONAL time scales=========
params = struct('R2a',R2a,'W2a',W2a,'Ea',Ea,'R2b',R2b,'W2b',W2b,'Eb',Eb);
paramNames = {'R2a','W2a','Ea','R2b','W2b','Eb'};

% Calculate seasonal means
% Get month numbers (same time vectors for all params, so doesn't matter which one is used)
mo = month(R2b.diel_dtd.date);
% Get year numbers
yr = year(R2b.diel_dtd.date);

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);
isSpring = ismember(mo,[3 4 5]);
isSummer = ismember(mo,[6 7 8]);
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

for m = 1:length(paramNames)
    % Add a column for month
    params.(paramNames{m}).diel_dtd.month = mo;

    % Add a column for season number
    params.(paramNames{m}).diel_dtd.season(indWinter) = 1;
    params.(paramNames{m}).diel_dtd.season(indSpring) = 2;
    params.(paramNames{m}).diel_dtd.season(indSummer) = 3;
    params.(paramNames{m}).diel_dtd.season(indFall) = 4;
    params.(paramNames{m}).diel_dtd.season = categorical(params.(paramNames{m}).diel_dtd.season);

    % Add a column for year
    % If statement so December is included in following winter
    for i = 1:height(R2b.diel_dtd)
        if params.(paramNames{m}).diel_dtd.month(i) == 12
            params.(paramNames{m}).diel_dtd.year(i) = categorical(yr(i) + 1);
        else
            params.(paramNames{m}).diel_dtd.year(i) = categorical(yr(i));
        end
    end

    ind_2021 = find(params.(paramNames{m}).diel_dtd.year == '2021');
    ind_2022 = find(params.(paramNames{m}).diel_dtd.year == '2022');
    ind_2023 = find(params.(paramNames{m}).diel_dtd.year == '2023');
    ind_2024 = find(params.(paramNames{m}).diel_dtd.year == '2024');

    seasonal_means_2021 = groupsummary(params.(paramNames{m}).diel_dtd(ind_2021,:),"season","mean",["GPP" "ER" "NEM"]);
    seasonal_means_2022 = groupsummary(params.(paramNames{m}).diel_dtd(ind_2022,:),"season","mean",["GPP" "ER" "NEM"]);
    seasonal_means_2023 = groupsummary(params.(paramNames{m}).diel_dtd(ind_2023,:),"season","mean",["GPP" "ER" "NEM"]);
    seasonal_means_2024 = groupsummary(params.(paramNames{m}).diel_dtd(ind_2024,:),"season","mean",["GPP" "ER" "NEM"]);
    
    seasonal_means_2021.year = repelem(categorical(2021), height(seasonal_means_2021))';
    seasonal_means_2022.year = repelem(categorical(2022), height(seasonal_means_2022))';
    seasonal_means_2023.year = repelem(categorical(2023), height(seasonal_means_2023))';
    seasonal_means_2024.year = repelem(categorical(2024), height(seasonal_means_2024))';

    % Compile seasonal means into one table
    seasonal_means.(paramNames{m}) = [seasonal_means_2021; seasonal_means_2022; seasonal_means_2023; seasonal_means_2024];
end

% xaxisLbl = {'Fall 2021',...
%     'Winter 2022','Spring 2022','Summer 2022','Fall 2022',...
%     'Winter 2023','Spring 2023','Summer 2023','Fall 2023',...
%     'Winter 2024','Spring 2024','Summer 2024'};

xaxisLbl = {'F ''21',...
    'W ''22','Sp ''22','Su ''22','F ''22',...
    'W ''23','Sp ''23','Su ''23','F ''23',...
    'W ''24','Sp ''24','Su ''24'};

% Plot the results (plot style from Palevsky et al. 2016, Fig. 2)
fig = figure;clf
fig.WindowState = 'maximized';
plot(seasonal_means.R2a.mean_GPP,'b.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
hold on
plot(seasonal_means.R2a.mean_ER,'k.-','MarkerSize',12,'LineWidth',1,'DisplayName','Ro & Hunt (2006) - const H')
plot(seasonal_means.R2a.mean_NEM,'r.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')

plot(seasonal_means.R2b.mean_GPP,'b.--','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
plot(seasonal_means.R2b.mean_ER,'k.--','MarkerSize',12,'LineWidth',1,'DisplayName','Ro & Hunt (2006) - H(t)')
plot(seasonal_means.R2b.mean_NEM,'r.--','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')

plot(seasonal_means.W2a.mean_GPP,'b*-','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')
plot(seasonal_means.W2a.mean_ER,'k*-','MarkerSize',8,'LineWidth',1,'DisplayName','Wanninkhof (2014) - const H')
plot(seasonal_means.W2a.mean_NEM,'r*-','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')

plot(seasonal_means.W2b.mean_GPP,'b*--','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')
plot(seasonal_means.W2b.mean_ER,'k*--','MarkerSize',8,'LineWidth',1,'DisplayName','Wanninkhof (2014) - H(t)')
plot(seasonal_means.W2b.mean_NEM,'r*--','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')

plot(seasonal_means.Ea.mean_GPP,'b^-','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')
plot(seasonal_means.Ea.mean_ER,'k^-','MarkerSize',6,'LineWidth',1,'DisplayName','Emerson et al. (2019) - const H')
plot(seasonal_means.Ea.mean_NEM,'r^-','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')

plot(seasonal_means.Eb.mean_GPP,'b^--','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')
plot(seasonal_means.Eb.mean_ER,'k^--','MarkerSize',6,'LineWidth',1,'DisplayName','Emerson et al. (2019) - H(t)')
plot(seasonal_means.Eb.mean_NEM,'r^--','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')

set(gca,'XTickLabel',xaxisLbl)
ylabel('Metabolic rate (mmol O_2 m^{-1} d^{-1})')
legend('show','location','southwest')
title([site,': Comparison of different k parameterizations (seasonal time scales)'])

%% ====Compare different parameterizations on ANNUAL time scales===========
wtreg_res_annual = groupsummary(wtreg_res,"DateTimeStamp","year","mean");

R2a_yearly = groupsummary(R2a.diel_dtd,"date","year","mean");
R2b_yearly = groupsummary(R2b.diel_dtd,"date","year","mean");

W2a_yearly = groupsummary(W2a.diel_dtd,"date","year","mean");
W2b_yearly = groupsummary(W2b.diel_dtd,"date","year","mean");

Ea_yearly = groupsummary(Ea.diel_dtd,"date","year","mean");
Eb_yearly = groupsummary(Eb.diel_dtd,"date","year","mean");

fig = figure;clf
fig.WindowState = 'maximized';

plot(R2a_yearly.year_date,R2a_yearly.mean_GPP,'b.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
hold on
plot(R2a_yearly.year_date,R2a_yearly.mean_ER,'k.-','MarkerSize',12,'LineWidth',1,'DisplayName','R&H-06 (Sc from W-14) - const H')
plot(R2a_yearly.year_date,R2a_yearly.mean_NEM,'r.-','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')

plot(R2b_yearly.year_date,R2b_yearly.mean_GPP,'b.--','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')
plot(R2b_yearly.year_date,R2b_yearly.mean_ER,'k.--','MarkerSize',12,'LineWidth',1,'DisplayName','R&H-06 (Sc from W-14) - H(t)')
plot(R2b_yearly.year_date,R2b_yearly.mean_NEM,'r.--','MarkerSize',12,'LineWidth',1,'HandleVisibility','off')

plot(W2a_yearly.year_date,W2a_yearly.mean_GPP,'b*-','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')
plot(W2a_yearly.year_date,W2a_yearly.mean_ER,'k*-','MarkerSize',8,'LineWidth',1,'DisplayName','W-14 - const H')
plot(W2a_yearly.year_date,W2a_yearly.mean_NEM,'r*-','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')

plot(W2b_yearly.year_date,W2b_yearly.mean_GPP,'b*--','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')
plot(W2b_yearly.year_date,W2b_yearly.mean_ER,'k*--','MarkerSize',8,'LineWidth',1,'DisplayName','W-14 - H(t)')
plot(W2b_yearly.year_date,W2b_yearly.mean_NEM,'r*--','MarkerSize',8,'LineWidth',1,'HandleVisibility','off')

plot(Ea_yearly.year_date,Ea_yearly.mean_GPP,'b^-','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')
plot(Ea_yearly.year_date,Ea_yearly.mean_ER,'k^-','MarkerSize',6,'LineWidth',1,'DisplayName','E-19 - const H')
plot(Ea_yearly.year_date,Ea_yearly.mean_NEM,'r^-','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')

plot(Eb_yearly.year_date,Eb_yearly.mean_GPP,'b^--','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')
plot(Eb_yearly.year_date,Eb_yearly.mean_ER,'k^--','MarkerSize',6,'LineWidth',1,'DisplayName','E-19 - H(t)')
plot(Eb_yearly.year_date,Eb_yearly.mean_NEM,'r^--','MarkerSize',6,'LineWidth',1,'HandleVisibility','off')

xlabel('UTC')
ylabel('Metabolic rate (mmol O_2 m^{-2} d^{-1})','FontSize',14)
set(gca,'FontSize',14,'LineWidth',2,'YColor','k')
legend('show','location','best')
title([site,': Comparison of different k parameterizations (yearly time scales)'])