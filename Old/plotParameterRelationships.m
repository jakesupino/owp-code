%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotParameterRelationships.m
% This script...
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 7/22/2024
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

% Find daily means
gull_dailyAvg = retime(gull_params,'daily','mean');
north_dailyAvg = retime(north_params,'daily','mean');
south_dailyAvg = retime(south_params,'daily','mean');

wspd_dailyAvg = retime(era5Dat,'daily','mean');
par_dailyAvg = retime(parDat,'daily','mean');

cd([rootpath,'figures\stats-analyses'])

%====Check distributions===================================================
% histogram(north_metab.NEM)

%====Create data table for all parameters/sites============================
% Gull
dt2 = dateshift(gull_metab.daystart_dt,'start','day');
gull_metab.daystart_dt = dt2;
TT_gull = synchronize(gull_dailyAvg,wspd_dailyAvg,par_dailyAvg,gull_metab);
TT_gull = removevars(TT_gull,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
TT_gull(TT_gull.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
TT_gull(TT_gull.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Get month numbers
mo = month(TT_gull.datetime_utc);  

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);   
isSpring = ismember(mo,[3 4 5]);    
isSummer = ismember(mo,[6 7 8]);    
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

TT_gull.season(indWinter) = 1;
TT_gull.season(indSpring) = 2;
TT_gull.season(indSummer) = 3;
TT_gull.season(indFall) = 4;
TT_gull.season = categorical(TT_gull.season);

TT_gull.site = repelem(2,height(TT_gull))';
TT_gull.site = categorical(TT_gull.site);

% North
dt2 = dateshift(north_metab.daystart_dt,'start','day');
north_metab.daystart_dt = dt2;
TT_north = synchronize(north_dailyAvg,wspd_dailyAvg,par_dailyAvg,north_metab);
TT_north = removevars(TT_north,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
TT_north(TT_north.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
TT_north(TT_north.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Get month numbers
mo = month(TT_north.datetime_utc);  

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);   
isSpring = ismember(mo,[3 4 5]);    
isSummer = ismember(mo,[6 7 8]);    
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

TT_north.season(indWinter) = 1;
TT_north.season(indSpring) = 2;
TT_north.season(indSummer) = 3;
TT_north.season(indFall) = 4;
TT_north.season = categorical(TT_north.season);

TT_north.site = repelem(1,height(TT_north))';
TT_north.site = categorical(TT_north.site);

% South
dt2 = dateshift(south_metab.daystart_dt,'start','day');
south_metab.daystart_dt = dt2;
TT_south = synchronize(south_dailyAvg,wspd_dailyAvg,par_dailyAvg,south_metab);
TT_south = removevars(TT_south,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
TT_south(TT_south.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
TT_south(TT_south.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Get month numbers
mo = month(TT_south.datetime_utc);  

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);   
isSpring = ismember(mo,[3 4 5]);    
isSummer = ismember(mo,[6 7 8]);    
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

TT_south.season(indWinter) = 1;
TT_south.season(indSpring) = 2;
TT_south.season(indSummer) = 3;
TT_south.season(indFall) = 4;
TT_south.season = categorical(TT_south.season);

TT_south.site = repelem(3,height(TT_south))';
TT_south.site = categorical(TT_south.site);

tbl_gull = timetable2table(TT_gull);
tbl_gull.datetime_utc = [];

tbl_north = timetable2table(TT_north);
tbl_north.datetime_utc = [];

tbl_south = timetable2table(TT_south);
tbl_south.datetime_utc = [];

data = [tbl_north;tbl_gull;tbl_south];

%% ====Plot metabolic rate vs. parameter (by season)=======================
% Make a for loop to make 3 (2x2) figs for N, G, S
% Add lin reg lines
% for i = 1:length(tbl)
tbl_all.north = tbl_north;
tbl_all.gull = tbl_gull;
tbl_all.south = tbl_south;

fn = fieldnames(tbl_all);
siteName = {'North','Gull','South'};

for i = 1:numel(fn)
    fig = figure;clf
    t = tiledlayout(2,2,'TileSpacing','compact');

    ax1 = nexttile;
    plot(tbl_all.(fn{i}).temperature,tbl_all.(fn{i}).ER,'.','color',rgb('lightgrey'))
    hold on
    plot(tbl_all.(fn{i}).temperature(indWinter),tbl_all.(fn{i}).ER(indWinter),'.k')
    title('Winter')

    ax2 = nexttile;
    plot(tbl_all.(fn{i}).temperature,tbl_all.(fn{i}).ER,'.','color',rgb('lightgrey'))
    hold on
    plot(tbl_all.(fn{i}).temperature(indSpring),tbl_all.(fn{i}).ER(indSpring),'.k')
    title('Spring')

    ax3 = nexttile;
    plot(tbl_all.(fn{i}).temperature,tbl_all.(fn{i}).ER,'.','color',rgb('lightgrey'))
    hold on
    plot(tbl_all.(fn{i}).temperature(indSummer),tbl_all.(fn{i}).ER(indSummer),'.k')
    title('Summer')

    ax4 = nexttile;
    plot(tbl_all.(fn{i}).temperature,tbl_all.(fn{i}).ER,'.','color',rgb('lightgrey'))
    hold on
    plot(tbl_all.(fn{i}).temperature(indFall),tbl_all.(fn{i}).ER(indFall),'.k')
    title('Fall')

    title(t,siteName(i),'fontsize',16)
    xlabel(t,'Temperature (^oC)','fontsize',14)
    ylabel(t,'ER (mmol O_2 m^{-2} d^{-1})','fontsize',14)

    fig.Position = [335,105,660,582];
end

% plot(tbl_south.DOsat(indSummer),tbl_south.pH(indSummer),'.')


%% ====Plot metabolic rate vs. parameter (all data)========================
% First attempt
site = questdlg('Site selection','Choose the site','North','Gull','South','Gull');

switch site
    case 'Gull'
        tbl = tbl_gull;
    case 'North'
        tbl = tbl_north;
    case 'South'
        tbl = tbl_south;
end

% Temperature
fig1 = figure(1);clf
t1 = tiledlayout(3,1,'TileSpacing','tight');
fig1.Position = [310,41,828,740];

ax1 = nexttile;
mdl = fitlm(tbl,'GPP~temperature');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.temperature,tbl.GPP,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Temperature (^oC)')
ylabel('GPP (mmol O_2 m^{-2} d^{-1})')
title(site)

ax2 = nexttile;
mdl = fitlm(tbl,'ER~temperature');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.temperature,tbl.ER,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Temperature (^oC)')
ylabel('ER (mmol O_2 m^{-2} d^{-1})')

ax3 = nexttile;
mdl = fitlm(tbl,'NEM~temperature');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.temperature,tbl.NEM,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Temperature (^oC)')
ylabel('NEM (mmol O_2 m^{-2} d^{-1})')

% Salinity
fig2 = figure(2);clf
t2 = tiledlayout(3,1,'TileSpacing','tight');
fig2.Position = [310,41,828,740];

ax1 = nexttile;
mdl = fitlm(tbl,'GPP~salinity');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.salinity,tbl.GPP,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Salinity')
ylabel('GPP (mmol O_2 m^{-2} d^{-1})')
title(site)

ax2 = nexttile;
mdl = fitlm(tbl,'ER~salinity');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.salinity,tbl.ER,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Salinity')
ylabel('ER (mmol O_2 m^{-2} d^{-1})')

ax3 = nexttile;
mdl = fitlm(tbl,'NEM~salinity');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.salinity,tbl.NEM,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Salinity')
ylabel('NEM (mmol O_2 m^{-2} d^{-1})')

% pH
fig3 = figure(3);clf
t3 = tiledlayout(3,1,'TileSpacing','tight');
fig3.Position = [310,41,828,740];

ax1 = nexttile;
mdl = fitlm(tbl,'GPP~pH');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.pH,tbl.GPP,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('pH')
ylabel('GPP (mmol O_2 m^{-2} d^{-1})')
title(site)

ax2 = nexttile;
mdl = fitlm(tbl,'ER~pH');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.pH,tbl.ER,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('pH')
ylabel('ER (mmol O_2 m^{-2} d^{-1})')

ax3 = nexttile;
mdl = fitlm(tbl,'NEM~pH');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.pH,tbl.NEM,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('pH')
ylabel('NEM (mmol O_2 m^{-2} d^{-1})')

% Chla
fig4 = figure(4);clf
t4 = tiledlayout(3,1,'TileSpacing','tight');
fig4.Position = [310,41,828,740];

ax1 = nexttile;
mdl = fitlm(tbl,'GPP~chla');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.chla,tbl.GPP,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Chl a (RFU)')
ylabel('GPP (mmol O_2 m^{-2} d^{-1})')
title(site)

ax2 = nexttile;
mdl = fitlm(tbl,'ER~chla');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.chla,tbl.ER,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Chl a (RFU)')
ylabel('ER (mmol O_2 m^{-2} d^{-1})')

ax3 = nexttile;
mdl = fitlm(tbl,'NEM~chla');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.chla,tbl.NEM,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Chl a (RFU)')
ylabel('NEM (mmol O_2 m^{-2} d^{-1})')

% Turbidity
fig5 = figure(5);clf
t5 = tiledlayout(3,1,'TileSpacing','tight');
fig5.Position = [310,41,828,740];

ax1 = nexttile;
mdl = fitlm(tbl,'GPP~turbidity');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.turbidity,tbl.GPP,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Turbidity (NTU)')
ylabel('GPP (mmol O_2 m^{-2} d^{-1})')
title(site)

ax2 = nexttile;
mdl = fitlm(tbl,'ER~turbidity');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.turbidity,tbl.ER,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Turbidity (NTU)')
ylabel('ER (mmol O_2 m^{-2} d^{-1})')

ax3 = nexttile;
mdl = fitlm(tbl,'NEM~turbidity');
coefs = mdl.Coefficients.Estimate;  % [intercept; slope]
plot(tbl.turbidity,tbl.NEM,'.')
hold on
refline(coefs(2),coefs(1))
xlabel('Turbidity (NTU)')
ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
