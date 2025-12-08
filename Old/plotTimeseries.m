%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotTimeseries.m
% This script...
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 8/6/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import final QC'd data & diel analysis results========================
site = 'gull';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
gull_params = finalQC;
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
gull_metab = diel_dtd;

site = 'north';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
north_params = finalQC;
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
north_metab = diel_dtd;

site = 'south';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
south_params = finalQC;
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
south_metab = diel_dtd;

clearvars finalQC diel_dtd diel_obs

%====Import windspeed and PAR data=========================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

cd([rootpath,'physical-data\final-dataset'])
load('par.mat')

%====Calculate daily means for each site===================================
gull_dailyAvg = retime(gull_params,'daily','mean');
north_dailyAvg = retime(north_params,'daily','mean');
south_dailyAvg = retime(south_params,'daily','mean');

wspd_dailyAvg = retime(era5Dat,'daily','mean');
par_dailyAvg = retime(parDat,'daily','mean');

cd([rootpath,'figures\stats-analyses'])

%====Create tables containing all data for each site=======================
% Gull
dt2 = dateshift(gull_metab.daystart_dt,'start','day');
gull_metab.daystart_dt = dt2;
daily_gull = synchronize(gull_dailyAvg,wspd_dailyAvg,par_dailyAvg,gull_metab);
daily_gull = removevars(daily_gull,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
daily_gull(daily_gull.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily_gull(daily_gull.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% North
dt2 = dateshift(north_metab.daystart_dt,'start','day');
north_metab.daystart_dt = dt2;
daily_north = synchronize(north_dailyAvg,wspd_dailyAvg,par_dailyAvg,north_metab);
daily_north = removevars(daily_north,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
daily_north(daily_north.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily_north(daily_north.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% South
dt2 = dateshift(south_metab.daystart_dt,'start','day');
south_metab.daystart_dt = dt2;
daily_south = synchronize(south_dailyAvg,wspd_dailyAvg,par_dailyAvg,south_metab);
daily_south = removevars(daily_south,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
daily_south(daily_south.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily_south(daily_south.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];
%%
%====Plots of daily means==================================================
daily.gull = daily_gull;
daily.north = daily_north;
daily.south = daily_south;

fig1 = figure(1);clf
t1 = tiledlayout(4,1,'TileSpacing','tight');

ax1 = nexttile;
yyaxis left
plot(daily.north.datetime_utc,daily.north.temperature)
hold on
plot(daily.gull.datetime_utc,daily.gull.temperature)
plot(daily.south.datetime_utc,daily.south.temperature)
set(gca,'XTickLabel',[])
ylabel('T (^oC)')

yyaxis right
plot(daily.north.datetime_utc,daily.north.pH)
hold on
plot(daily.gull.datetime_utc,daily.gull.pH)
plot(daily.south.datetime_utc,daily.south.pH)
ylim([7.4 8.6])
ylabel('pH')

ax2 = nexttile;
yyaxis left
plot(daily.north.datetime_utc,daily.north.DOsat)
hold on
plot(daily.gull.datetime_utc,daily.gull.DOsat)
plot(daily.south.datetime_utc,daily.south.DOsat)
set(gca,'XTickLabel',[])
ylabel('DO (%sat)')

yyaxis right
plot(daily.north.datetime_utc,daily.north.salinity)
hold on
plot(daily.gull.datetime_utc,daily.gull.salinity)
plot(daily.south.datetime_utc,daily.south.salinity)
ylabel('S (psu)')

ax3 = nexttile;
yyaxis left
plot(daily.north.datetime_utc,daily.north.chla,'DisplayName','North')
hold on
plot(daily.gull.datetime_utc,daily.gull.chla,'DisplayName','Gull')
plot(daily.south.datetime_utc,daily.south.chla,'DisplayName','South')
set(gca,'XTickLabel',[])
ylabel('Chl a (RFU)')
legend('show','Location','northwest')

yyaxis right
plot(daily.north.datetime_utc,daily.north.turbidity,'HandleVisibility','off')
hold on
plot(daily.gull.datetime_utc,daily.gull.turbidity,'HandleVisibility','off')
plot(daily.south.datetime_utc,daily.south.turbidity,'HandleVisibility','off')
ylabel('Turbidity (NTU)')

ax4 = nexttile;
yyaxis left
plot(daily.north.datetime_utc,daily.north.GPP,'-k')
hold on
plot(daily.gull.datetime_utc,daily.gull.GPP,'--k')
plot(daily.south.datetime_utc,daily.south.GPP,':k')
plot(daily.north.datetime_utc,daily.north.ER,'-k')
plot(daily.gull.datetime_utc,daily.gull.ER,'--k')
plot(daily.south.datetime_utc,daily.south.ER,':k')
ylabel('ER and GPP')

yyaxis right
plot(daily.north.datetime_utc,daily.north.NEM)
hold on
plot(daily.gull.datetime_utc,daily.gull.NEM)
plot(daily.south.datetime_utc,daily.south.NEM)
ylabel('NEM')

title(t1,'Daily Means','FontSize',18,'FontWeight','bold')

%%
%====Plots of monthly means================================================
monthly_north = groupsummary(daily_north,'datetime_utc','month','mean');
monthly_gull = groupsummary(daily_gull,'datetime_utc','month','mean');
monthly_south = groupsummary(daily_south,'datetime_utc','month','mean');

% HERE
monthly_north.month_datetime_utc = datetime(char(monthly_north.month_datetime_utc),'Format','MMM-uuuu');
monthly_north = table2timetable(monthly_north);
monthly_south.month_datetime_utc = datetime(char(monthly_south.month_datetime_utc),'Format','MMM-uuuu');
monthly_south = table2timetable(monthly_south);
monthly_gull.month_datetime_utc = datetime(char(monthly_gull.month_datetime_utc),'Format','MMM-uuuu');
monthly_gull = table2timetable(monthly_gull);

monthly.north = monthly_north;
monthly.gull = monthly_gull;
monthly.south = monthly_south;

fig2 = figure(2);clf
t2 = tiledlayout(4,1,'TileSpacing','tight');

ax1 = nexttile;
yyaxis left
plot(monthly.north.month_datetime_utc,monthly.north.mean_temperature)
hold on
plot(monthly.gull.month_datetime_utc,monthly.gull.mean_temperature)
plot(monthly.south.month_datetime_utc,monthly.south.mean_temperature)
set(gca,'XTickLabel',[])
ylabel('T (^oC)')

yyaxis right
plot(monthly.north.month_datetime_utc,monthly.north.mean_pH)
hold on
plot(monthly.gull.month_datetime_utc,monthly.gull.mean_pH)
plot(monthly.south.month_datetime_utc,monthly.south.mean_pH)
ylim([7.4 8.6])
ylabel('pH')

ax2 = nexttile;
yyaxis left
plot(monthly.north.month_datetime_utc,monthly.north.mean_DOsat)
hold on
plot(monthly.gull.month_datetime_utc,monthly.gull.mean_DOsat)
plot(monthly.south.month_datetime_utc,monthly.south.mean_DOsat)
set(gca,'XTickLabel',[])
ylabel('DO (%sat)')

yyaxis right
plot(monthly.north.month_datetime_utc,monthly.north.mean_salinity)
hold on
plot(monthly.gull.month_datetime_utc,monthly.gull.mean_salinity)
plot(monthly.south.month_datetime_utc,monthly.south.mean_salinity)
ylim([28 38])
ylabel('S (psu)')

ax3 = nexttile;
yyaxis left
plot(monthly.north.month_datetime_utc,monthly.north.mean_chla,'DisplayName','North')
hold on
plot(monthly.gull.month_datetime_utc,monthly.gull.mean_chla,'DisplayName','Gull')
plot(monthly.south.month_datetime_utc,monthly.south.mean_chla,'DisplayName','South')
set(gca,'XTickLabel',[])
ylabel('Chl a (RFU)')
legend('show','Location','northwest')

yyaxis right
plot(monthly.north.month_datetime_utc,monthly.north.mean_turbidity,'HandleVisibility','off')
hold on
plot(monthly.gull.month_datetime_utc,monthly.gull.mean_turbidity,'HandleVisibility','off')
plot(monthly.south.month_datetime_utc,monthly.south.mean_turbidity,'HandleVisibility','off')
ylabel('Turbidity (NTU)')

ax4 = nexttile;
yyaxis left
plot(monthly.north.month_datetime_utc,monthly.north.mean_GPP,'-k')
hold on
plot(monthly.gull.month_datetime_utc,monthly.gull.mean_GPP,'--k')
plot(monthly.south.month_datetime_utc,monthly.south.mean_GPP,':k')
plot(monthly.north.month_datetime_utc,monthly.north.mean_ER,'-k')
plot(monthly.gull.month_datetime_utc,monthly.gull.mean_ER,'--k')
plot(monthly.south.month_datetime_utc,monthly.south.mean_ER,':k')
ylim([-300 300])
ylabel('ER and GPP')

yyaxis right
plot(monthly.north.month_datetime_utc,monthly.north.mean_NEM)
hold on
plot(monthly.gull.month_datetime_utc,monthly.gull.mean_NEM)
plot(monthly.south.month_datetime_utc,monthly.south.mean_NEM)
ylabel('NEM')

labels = string(ax4.XAxis.TickLabels); % extract
labels(2:3:end) = NaN; % remove every second one
labels(3:3:end) = NaN; % remove every third one
ax4.XAxis.TickLabels = labels; % set
title(t2,'Monthly Means','FontSize',18,'FontWeight','bold')

%%
%====Plots of seasonal means================================================
% Create vector of dates to use for extracting seasonal months
startdate = datetime(2021,06,01,'Format',"MMM-uuuu");
enddate = datetime(2024,06,04,'Format',"MMM-uuuu");
dates = startdate:calmonths(1):enddate;
winterMonth = month(dates) == 12 | month(dates) == 1 | month(dates) == 2;
springMonth = month(dates) == 3 | month(dates) == 4 | month(dates) == 5;
summerMonth = month(dates) == 6 | month(dates) == 7 | month(dates) == 8;
fallMonth = month(dates) == 9 | month(dates) == 10 | month(dates) == 11;

fn = fieldnames(monthly);
for j = 1:numel(fn)
    % Extract seasonal data from monthly timetables
    winter = monthly.(fn{j})(winterMonth,:);
    spring = monthly.(fn{j})(springMonth,:);
    summer = monthly.(fn{j})(summerMonth,:);
    fall = monthly.(fn{j})(fallMonth,:);

    % Create grouping variable for each season each year
    winterYr = year(dates(winterMonth))';
    winterYr(month(dates(winterMonth)) == 12) = winterYr(month(dates(winterMonth))==12)+1;
    springYr = year(dates(springMonth))';
    summerYr = year(dates(summerMonth))';
    fallYr = year(dates(fallMonth))';

    % Find winter means
    G = findgroups(winterYr);
    winterYrMean = NaN(length(unique(winterYr)),13);
    for i = 1:width(winterYrMean)
        winterYrMean(:,i) = splitapply(@nanmean,winter(:,i+1),G);
    end
    winterYrMean = array2table(winterYrMean);
    winterYrMean.Properties.VariableNames = {'depth','salinity','temperature','DOconc','DOsat','pH','chla','turbidity','wspd','par','GPP','ER','NEM'};
    winterYrMean.datetime = datetime(unique(winterYr),1,1); % Set dummy "winter month" as January
    winterYrMean = table2timetable(winterYrMean);

    % Find spring means
    G = findgroups(springYr);
    springYrMean = NaN(length(unique(springYr)),13);
    for i = 1:width(springYrMean)
        springYrMean(:,i) = splitapply(@nanmean,spring(:,i+1),G);
    end
    springYrMean = array2table(springYrMean);
    springYrMean.Properties.VariableNames = {'depth','salinity','temperature','DOconc','DOsat','pH','chla','turbidity','wspd','par','GPP','ER','NEM'};
    springYrMean.datetime = datetime(unique(springYr),4,1); % Set dummy "spring month" as April
    springYrMean = table2timetable(springYrMean);

    % Find summer means
    G = findgroups(summerYr);
    summerYrMean = NaN(length(unique(summerYr)),13);
    for i = 1:width(summerYrMean)
        summerYrMean(:,i) = splitapply(@nanmean,summer(:,i+1),G);
    end
    summerYrMean = array2table(summerYrMean);
    summerYrMean.Properties.VariableNames = {'depth','salinity','temperature','DOconc','DOsat','pH','chla','turbidity','wspd','par','GPP','ER','NEM'};
    summerYrMean.datetime = datetime(unique(summerYr),7,1); % Set dummy "summer month" as July
    summerYrMean = table2timetable(summerYrMean);

    % Find fall means
    G = findgroups(fallYr);
    fallYrMean = NaN(length(unique(fallYr)),13);
    for i = 1:width(fallYrMean)
        fallYrMean(:,i) = splitapply(@nanmean,fall(:,i+1),G);
    end
    fallYrMean = array2table(fallYrMean);
    fallYrMean.Properties.VariableNames = {'depth','salinity','temperature','DOconc','DOsat','pH','chla','turbidity','wspd','par','GPP','ER','NEM'};
    fallYrMean.datetime = datetime(unique(fallYr),10,1); % Set dummy "fall month" as October
    fallYrMean = table2timetable(fallYrMean);

    % Merge individual season tables
    seasonal.(fn{j}) = sortrows([winterYrMean; springYrMean; summerYrMean; fallYrMean]);
end

fig3 = figure(3);clf
t3 = tiledlayout(4,2,'TileSpacing','tight');
xaxisLbl = {'Summer 2021','Fall 2021','Winter 2022','Spring 2022','Summer 2022','Fall 2022','Winter 2023','Spring 2023','Summer 2023','Fall 2023','Winter 2024','Spring 2024','Summer 2024'};

nexttile
T_mat = [seasonal.north.temperature, seasonal.gull.temperature, seasonal.south.temperature];
bar(seasonal.north.datetime,T_mat)
set(gca,'XTickLabel',[])
ylabel('T (^oC)')

nexttile
pH_mat = [seasonal.north.pH, seasonal.gull.pH, seasonal.south.pH];
bar(seasonal.north.datetime,pH_mat)
ylim([7.4 8.6])
set(gca,'XTickLabel',[])
ylabel('pH')

nexttile
DOsat_mat = [seasonal.north.DOsat, seasonal.gull.DOsat, seasonal.south.DOsat];
bar(seasonal.north.datetime,DOsat_mat)
ylim([0 120])
set(gca,'XTickLabel',[])
ylabel('DO (%sat)')

nexttile
S_mat = [seasonal.north.salinity, seasonal.gull.salinity, seasonal.south.salinity];
bar(seasonal.north.datetime,S_mat)
ylim([26 36])
set(gca,'XTickLabel',[])
ylabel('S (psu)')

nexttile
chla_mat = [seasonal.north.chla, seasonal.gull.chla, seasonal.south.chla];
bar(seasonal.north.datetime,chla_mat)
set(gca,'XTickLabel',[])
ylim([0 18])
ylabel('Chl a (RFU)')
legend({'North','Gull','South'},'Location','northwest')

nexttile
turb_mat = [seasonal.north.turbidity, seasonal.gull.turbidity, seasonal.south.turbidity];
bar(seasonal.north.datetime,turb_mat)
set(gca,'XTickLabel',[])
ylabel('Turbidity (NTU)')

nexttile
GPP_mat = [seasonal.north.GPP, seasonal.gull.GPP, seasonal.south.GPP];
bar(seasonal.north.datetime,GPP_mat)
hold on
ER_mat = [seasonal.north.ER, seasonal.gull.ER, seasonal.south.ER];
b = bar(seasonal.north.datetime,ER_mat);
ylim([-200 200])
ylabel('ER and GPP')
set(gca,'XTickLabel',xaxisLbl)
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];

nexttile
NEM_mat = [seasonal.north.NEM, seasonal.gull.NEM, seasonal.south.NEM];
bar(seasonal.north.datetime,NEM_mat)
ylim([-25 10])
ylabel('NEM')
set(gca,'XTickLabel',xaxisLbl)


title(t3,'Seasonal Means','FontSize',18,'FontWeight','bold')