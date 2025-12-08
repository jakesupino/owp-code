%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fillGaps.m
% This script calculates a smoothed climatological day-of-year mean for 
% each parameter over the whole time series and uses it to fill in missing 
% days of data.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 3/6/2025
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%==========================================================================
% Import data and create table of daily means for each site
%==========================================================================
site = {'north','gull','south'};

for i = 1:length(site)
    cd([rootpath,'open-water-platform-data\',site{i},'\cleaned\final-qc'])
    load('finalQC.mat');
    params.(site{i}) = finalQC;
    cd([rootpath,'diel-method\uncertainty-analysis\',site{i}])
    load('MonteCarloResults')
    metab.(site{i}) = diel_dtd_MC;
end
clearvars finalQC diel_dtd diel_obs

% Find daily tidal range for each site
gull_range = groupsummary(params.gull,"datetime_utc","day","range","depth");
gull_range(find(gull_range.GroupCount < 144),:) = []; % Remove days that don't have all depth measurements (full day = 144 points)
gull_range.day_datetime_utc = datetime(cellstr(gull_range.day_datetime_utc),'TimeZone','UTC');
gull_range = table2timetable(gull_range);
gull_range.GroupCount = [];
gull_range.Properties.VariableNames = {'tidal'};

north_range = groupsummary(params.north,"datetime_utc","day","range","depth");
north_range(find(north_range.GroupCount < 144),:) = []; % Remove days that don't have all depth measurements (full day = 144 points)
north_range.day_datetime_utc = datetime(cellstr(north_range.day_datetime_utc),'TimeZone','UTC');
north_range = table2timetable(north_range);
north_range.GroupCount = [];
north_range.Properties.VariableNames = {'tidal'};

south_range = groupsummary(params.south,"datetime_utc","day","range","depth");
south_range(find(south_range.GroupCount < 144),:) = []; % Remove days that don't have all depth measurements (full day = 144 points)
south_range.day_datetime_utc = datetime(cellstr(south_range.day_datetime_utc),'TimeZone','UTC');
south_range = table2timetable(south_range); 
south_range.GroupCount = [];
south_range.Properties.VariableNames = {'tidal'};

%====Import windspeed data=========================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

%====Calculate daily means for each site===================================
gull_dailyAvg = retime(params.gull,'daily','mean');
north_dailyAvg = retime(params.north,'daily','mean');
south_dailyAvg = retime(params.south,'daily','mean');

wspd_dailyAvg = retime(era5Dat,'daily','mean');

%====Create time tables for each site======================================
%----Gull------------------------------------------------------------------
dt2 = dateshift(metab.gull.daystart_dt,'start','day');
metab.gull.daystart_dt = dt2;
daily.gull = synchronize(gull_dailyAvg,gull_range,wspd_dailyAvg,metab.gull);
daily.gull = removevars(daily.gull,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});
daily.gull(daily.gull.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily.gull(daily.gull.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Set anomalous GPP and ER values as NaN
anomER = find(daily.gull.ER_avg > 0);
anomGPP = find(daily.gull.GPP_avg < 0);
daily.gull.ER_avg(anomER) = NaN;
daily.gull.ER_sd(anomER) = NaN;
daily.gull.GPP_avg(anomGPP) = NaN;
daily.gull.GPP_sd(anomGPP) = NaN;

%----North-----------------------------------------------------------------
dt2 = dateshift(metab.north.daystart_dt,'start','day');
metab.north.daystart_dt = dt2;
daily.north = synchronize(north_dailyAvg,north_range,wspd_dailyAvg,metab.north);
daily.north = removevars(daily.north,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});
daily.north(daily.north.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily.north(daily.north.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Set anomalous GPP and ER values as NaN
anomER = find(daily.north.ER_avg > 0);
anomGPP = find(daily.north.GPP_avg < 0);
daily.north.ER_avg(anomER) = NaN;
daily.north.ER_sd(anomER) = NaN;
daily.north.GPP_avg(anomGPP) = NaN;
daily.north.GPP_sd(anomGPP) = NaN;

%----South-----------------------------------------------------------------
dt2 = dateshift(metab.south.daystart_dt,'start','day');
metab.south.daystart_dt = dt2;
daily.south = synchronize(south_dailyAvg,south_range,wspd_dailyAvg,metab.south);
daily.south = removevars(daily.south,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});
daily.south(daily.south.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily.south(daily.south.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Set anomalous GPP and ER values as NaN
anomER = find(daily.south.ER_avg > 0);
anomGPP = find(daily.south.GPP_avg < 0);
daily.south.ER_avg(anomER) = NaN;
daily.south.ER_sd(anomER) = NaN;
daily.south.GPP_avg(anomGPP) = NaN;
daily.south.GPP_sd(anomGPP) = NaN;

% See triad scheme in https://files-aje-com.s3.amazonaws.com/www/row/_assets/docs/Using_Color_In_Your_Manuscript_Figures.pdf
% and https://www.rapidtables.com/convert/color/rgb-to-hex.html?r=86&g=180&b=233
north_clr = '#5D3A9B';
gull_clr = '#019E73';
south_clr = '#D55E00';

t = datetime(2021,01,01):calyears(1):datetime(2023,01,01);
yr = year(t);

xaxisLbl = t;
xtk = 1:3;

cd([rootpath,'figures'])

%==========================================================================
% Gap fill each variable
%==========================================================================
varNames = {"temperature","salinity","DOsat","pH","depth","tidal","GPP_avg","ER_avg","NEM_avg"};
lblNames = {"Temperature","Salinity","DOsat","pH","Depth","Tidal","GPP","ER","NEM"};

% Create duplicate table of daily values to add dummy data to
for i = 1:length(site)
    daily_filled.(site{i}) = timetable(daily.(site{i}).datetime_utc,daily.(site{i}).temperature,...
        daily.(site{i}).salinity,daily.(site{i}).DOsat,daily.(site{i}).pH,daily.(site{i}).depth,...
        daily.(site{i}).tidal,daily.(site{i}).GPP_avg,daily.(site{i}).ER_avg,daily.(site{i}).NEM_avg,...
        'VariableNames',{'temperature','salinity','DOsat','pH','depth','tidal','GPP_avg','ER_avg','NEM_avg'});
end

for j = 1:length(varNames)
    for i = 1:length(site)
        % Calculate 10-day moving mean and use climatological mean of smoothed data as dummy data
        dummyVal = movmean(daily.(site{i}).(varNames{j}),days(10),"omitnan",'SamplePoints',daily.(site{i}).datetime_utc);
        dummyTbl = timetable(daily.(site{i}).datetime_utc,dummyVal,'VariableNames',(varNames{j}));
        doy_dummy = groupsummary(dummyTbl,"Time","dayofyear","mean");
        doy_dummy = renamevars(doy_dummy,"dayofyear_Time","doy");
        doy_dummy.doy = double(doy_dummy.doy);
        % Fill in gaps with dummy data
        ind_gap.(site{i}) = find(isnan(daily.(site{i}).(varNames{j})));
        doy_gap = day(daily.(site{i}).datetime_utc(ind_gap.(site{i})),'dayofyear');
        daily_filled.(site{i}).(varNames{j})(ind_gap.(site{i})) = doy_dummy.(strcat('mean_',varNames{j}))(doy_gap);
    end
    % Plot results all sites for each variable
    fig(j) = figure;clf
    fig(j).WindowState = 'maximized';
    t = tiledlayout(3,1,'Padding','compact','TileSpacing','tight');
    nexttile
    plot(daily.north.datetime_utc,daily.north.(varNames{j}),'.','Color',north_clr,'DisplayName','Data')
    hold on
    plot(daily_filled.north.Time(ind_gap.north),daily_filled.north.(varNames{j})(ind_gap.north),'x','Color',rgb('goldenrod'),'DisplayName','Climatological mean')
    ylabel(lblNames{j})
    legend('show','location','best')

    nexttile
    plot(daily.gull.datetime_utc,daily.gull.(varNames{j}),'.','Color',gull_clr,'Color',gull_clr,'DisplayName','Data')
    hold on
    plot(daily_filled.gull.Time(ind_gap.gull),daily_filled.gull.(varNames{j})(ind_gap.gull),'x','Color',rgb('goldenrod'),'DisplayName','Climatological mean')
    ylabel(lblNames{j})
    legend('show','location','best')

    nexttile
    plot(daily.south.datetime_utc,daily.south.(varNames{j}),'.','Color',south_clr,'DisplayName','Data')
    hold on
    plot(daily_filled.south.Time(ind_gap.south),daily_filled.south.(varNames{j})(ind_gap.south),'x','Color',rgb('goldenrod'),'DisplayName','Climatological mean')
    ylabel(lblNames{j})
    legend('show','location','best')

    title(t,'Gap-filled daily data','fontsize',20)
    
    pause
end
      
% Option to save results and figure
option = questdlg('Save the results?','Save results','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\all-sites'])
        save('gapFilled.mat','daily_filled')
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end

option = questdlg('Save the figures?','Save figures','Yes','No','Yes');
switch option
    case 'Yes'
        for j = 1:length(varNames)
            cd([rootpath,'figures\open-water-platform\all-sites\gap-filled'])
            saveas(fig(j),strcat(lblNames{j},'.fig'))
            saveas(fig(j),strcat(lblNames{j},'.png'))
        end
        disp('Figures saved!')
    case 'No'
        disp('Figures not saved.')
end