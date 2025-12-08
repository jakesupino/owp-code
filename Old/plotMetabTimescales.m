%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotMetabTimescales.m
% This script...
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 11/6/24
% Last updated: 11/11/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import Monte Carlo diel analysis results from all 3 OWPs==============
fn = {'north','gull','south'};

for i = 1:length(fn)
    cd([rootpath,'diel-method\uncertainty-analysis\',(fn{i})])
    load('MonteCarloResults.mat')
    daily.(fn{i}) = diel_dtd_MC;

    % Remove times with anomalous GPP and ER values
    anomER = find(daily.(fn{i}).ER_avg > 0);
    anomGPP = find(daily.(fn{i}).GPP_avg < 0);
    daily.(fn{i})([anomER;anomGPP],:) = [];

    monthly.(fn{i}) = retime(daily.(fn{i}),"monthly","mean");
    yearly.(fn{i}) = retime(daily.(fn{i}),"yearly","mean");
end

%====Calculate means across sites on different time scales=================
% Daily means across sites
tbl_syn = synchronize(daily.north,daily.gull,daily.south);
daily_all = timetable(tbl_syn.date,tbl_syn.NEM_avg_1,tbl_syn.NEM_avg_2,tbl_syn.NEM_avg_3,tbl_syn.NEM_sd_1,tbl_syn.NEM_sd_2,tbl_syn.NEM_sd_3);
daily_all.Properties.VariableNames = {'NEM_avg_N' 'NEM_avg_G' 'NEM_avg_S' 'NEM_sd_N' 'NEM_sd_G' 'NEM_sd_S'};
site_avg = mean([daily_all.NEM_avg_N daily_all.NEM_avg_G daily_all.NEM_avg_S],2,'omitmissing');
daily_mean = timetable(daily_all.Time,site_avg,'VariableNames',"NEM");

% Monthly means across sites
monthly_all = retime(daily_all,"monthly","mean");
site_avg = mean([monthly_all.NEM_avg_N monthly_all.NEM_avg_G monthly_all.NEM_avg_S],2,'omitmissing');
monthly_mean = timetable(monthly_all.Time,site_avg,'VariableNames',"NEM");

% Quasi-annual means across sites (Sep-June of each year)
y1=mean(monthly_mean.NEM(1:10),'omitmissing');          % Sep 2021 - Jun 2022 (November is missing)
y2=mean(monthly_mean.NEM([13:14 16:21]),'omitmissing'); % Sep 2022 - Jun 2023 (omitting November)
y3=mean(monthly_mean.NEM([25:26 28:33]),'omitmissing'); % Sep 2023 - Jun 2024 (omitting November)

cd([rootpath,'figures\diel-analysis\monte-carlo\'])

fig = figure;clf
t = tiledlayout(3,1,'TileSpacing','compact');

nexttile
bar(daily_mean.Time,daily_mean.NEM,'k')
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'Daily','Location','SouthEast')
legend boxoff

nexttile
bar(monthly_mean.Time,monthly_mean.NEM,'k')
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'Monthly','Location','SouthEast')
legend boxoff

nexttile
bar({'2021','2022','2023'},[y1 y2 y3],'k')
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'Quasi-annual','Location','SouthEast')
legend boxoff

ylabel(t,'NEM (mmol O_2 m^{-2} d^{-1})','FontSize',14)

set(gcf,'Position',[-1000 -200 400 300])
fig.Units               = 'centimeters';
fig.Position(3)         = 18;
fig.Position(4)         = 22;

%% Old way
% Daily means across sites
% daily_NG = outerjoin(daily.north,daily.gull);
% daily_mean = table(daily_NG.date,mean([daily_NG.NEM_avg_left daily_NG.NEM_avg_right],2,'omitmissing'));
% daily_mean.Properties.VariableNames = {'date','NEM_avg'};
% daily_mean = table2timetable(daily_mean);
% 
% daily_NGS = outerjoin(daily_mean,daily.south);
% daily_mean = table(daily_NGS.date,mean([daily_NGS.NEM_avg_daily_mean daily_NGS.NEM_avg_right],2,'omitmissing'));
% daily_mean.Properties.VariableNames = {'date','NEM_avg'};
% daily_mean_old = table2timetable(daily_mean);
% 
% monthly_mean_old = retime(daily_mean_old,"monthly","mean");
% 
% y1=mean(monthly_mean.NEM(1:10),'omitmissing');   % Sep 2021 - Jun 2022
% y2=mean(monthly_mean.NEM(13:22),'omitmissing');  % Sep 2022 - Jun 2023
% y3=mean(monthly_mean.NEM(25:34),'omitmissing');  % Sep 2023 - Jun 2024

% % Yearly means across sites, using ALL available months per year
% yearly_all = retime(daily_all,"yearly","mean");
% site_avg = mean([yearly_all.NEM_avg_N yearly_all.NEM_avg_G yearly_all.NEM_avg_S],2,'omitmissing');
% yearly_mean = timetable(yearly_all.Time,site_avg,'VariableNames',"NEM");

% %% Attempt to find mean and sd using the daily and monthly all sites table 
% % (DON'T THINK THIS WORKS B/C SOME SITES HAVE A LOT OF MISSING MONTHS DURING 2021-2022, GIVING MORE WEIGHT TO MONTHS THAT HAVE MORE MEASUREMENTS
% y1_mean = mean([daily_all.NEM_avg_N(1:200) daily_all.NEM_avg_G(1:200),daily_all.NEM_avg_S(1:200)],"all",'omitmissing');
% y1_sd = std([daily_all.NEM_avg_N(1:200) daily_all.NEM_avg_G(1:200),daily_all.NEM_avg_S(1:200)],0,"all",'omitmissing');
% 
% y2_mean = mean([daily_all.NEM_avg_N(306:560) daily_all.NEM_avg_G(306:560),daily_all.NEM_avg_S(306:560)],"all",'omitmissing');
% y2_sd = std([daily_all.NEM_avg_N(306:560) daily_all.NEM_avg_G(306:560),daily_all.NEM_avg_S(306:560)],0,"all",'omitmissing');
% 
% y3_mean = mean([daily_all.NEM_avg_N(667:921) daily_all.NEM_avg_G(667:921),daily_all.NEM_avg_S(667:921)],"all",'omitmissing');
% y3_sd = std([daily_all.NEM_avg_N(667:921) daily_all.NEM_avg_G(667:921),daily_all.NEM_avg_S(667:921)],0,"all",'omitmissing');

% mean([monthly_all.NEM_avg_N(1:10) monthly_all.NEM_avg_G(1:10) monthly_all.NEM_avg_S(1:10)],'omitmissing')
% std([monthly_all.NEM_avg_N(1:10) monthly_all.NEM_avg_G(1:10) monthly_all.NEM_avg_S(1:10)],'omitmissing')
