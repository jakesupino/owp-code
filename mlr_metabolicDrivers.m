%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mlr_metabolicDrivers.m
% This script creates multiple linear regression models to evaluate the influence
% of different environmental drivers on the estimated metabolic rates.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 7/22/2024
% Last updated: 11/6/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%==========================================================================
% Import data
%==========================================================================

%====Import final QC'd data & diel analysis results========================
fn = {'north','gull','south'};

for i = 1:length(fn)
    cd([rootpath,'open-water-platform-data\',fn{i},'\cleaned\final-qc'])
    load('finalQC.mat');
    params.(fn{i}) = finalQC;
    cd([rootpath,'diel-method\uncertainty-analysis\',fn{i}])
    load('MonteCarloResults')
    metab.(fn{i}) = diel_dtd_MC;
end
clearvars finalQC diel_dtd diel_obs

% Find daily tidal range for each site
gull_range = groupsummary(params.gull,"datetime_utc","day","range","depth");
gull_range.day_datetime_utc = datetime(cellstr(gull_range.day_datetime_utc),'TimeZone','UTC');
gull_range = table2timetable(gull_range);
gull_range.GroupCount = [];
gull_range.Properties.VariableNames = {'tidal'};

north_range = groupsummary(params.north,"datetime_utc","day","range","depth");
north_range.day_datetime_utc = datetime(cellstr(north_range.day_datetime_utc),'TimeZone','UTC');
north_range = table2timetable(north_range);
north_range.GroupCount = [];
north_range.Properties.VariableNames = {'tidal'};

south_range = groupsummary(params.south,"datetime_utc","day","range","depth");
south_range.day_datetime_utc = datetime(cellstr(south_range.day_datetime_utc),'TimeZone','UTC');
south_range = table2timetable(south_range);
south_range.GroupCount = [];
south_range.Properties.VariableNames = {'tidal'};

%====Import windspeed and PAR data=========================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

cd([rootpath,'physical-data\final-dataset'])
load('par.mat')

%====Calculate daily means for each site===================================
gull_dailyAvg = retime(params.gull,'daily','mean');
north_dailyAvg = retime(params.north,'daily','mean');
south_dailyAvg = retime(params.south,'daily','mean');

wspd_dailyAvg = retime(era5Dat,'daily','mean');
par_dailyAvg = retime(parDat,'daily','mean');

cd([rootpath,'figures\stats-analyses'])

%====Create time tables for each site======================================
%----Gull------------------------------------------------------------------
dt2 = dateshift(metab.gull.daystart_dt,'start','day');
metab.gull.daystart_dt = dt2;
gull_TT = synchronize(gull_dailyAvg,gull_range,wspd_dailyAvg,par_dailyAvg,metab.gull);
gull_TT = removevars(gull_TT,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg','datetime_local','Tair','light_lux'});
gull_TT(gull_TT.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
gull_TT(gull_TT.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Remove timepoints with anomalous GPP and ER values
anomER = find(gull_TT.ER_avg > 0);
anomGPP = find(gull_TT.GPP_avg < 0);
gull_TT([anomER;anomGPP],:) = [];

% Remove row if any explanatory or response variable is NaN
gull_TT = rmmissing(gull_TT,'DataVariables',{'tidal','salinity','temperature','DOsat','pH','wspd','GPP_avg','ER_avg','NEM_avg'});

% Get month numbers
mo = month(gull_TT.datetime_utc);

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);
isSpring = ismember(mo,[3 4 5]);
isSummer = ismember(mo,[6 7 8]);
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

% Add a column for season number
gull_TT.season(indWinter) = 1;
gull_TT.season(indSpring) = 2;
gull_TT.season(indSummer) = 3;
gull_TT.season(indFall) = 4;
gull_TT.season = categorical(gull_TT.season);

% Add a column for site number
gull_TT.site = repelem(2,height(gull_TT))';
gull_TT.site = categorical(gull_TT.site);

%----North-----------------------------------------------------------------
dt2 = dateshift(metab.north.daystart_dt,'start','day');
metab.north.daystart_dt = dt2;
north_TT = synchronize(north_dailyAvg,north_range,wspd_dailyAvg,par_dailyAvg,metab.north);
north_TT = removevars(north_TT,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg','datetime_local','Tair','light_lux'});
north_TT(north_TT.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
north_TT(north_TT.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Remove timepoints with anomalous GPP and ER values
anomER = find(north_TT.ER_avg > 0);
anomGPP = find(north_TT.GPP_avg < 0);
north_TT([anomER;anomGPP],:) = [];

% Remove row if any explanatory or response variable is NaN
north_TT = rmmissing(north_TT,'DataVariables',{'tidal','salinity','temperature','DOsat','pH','wspd','GPP_avg','ER_avg','NEM_avg'});

% Get month numbers
mo = month(north_TT.datetime_utc);

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);
isSpring = ismember(mo,[3 4 5]);
isSummer = ismember(mo,[6 7 8]);
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

% Add a column for season number
north_TT.season(indWinter) = 1;
north_TT.season(indSpring) = 2;
north_TT.season(indSummer) = 3;
north_TT.season(indFall) = 4;
north_TT.season = categorical(north_TT.season);

% Add a column for site number
north_TT.site = repelem(1,height(north_TT))';
north_TT.site = categorical(north_TT.site);

%----South-----------------------------------------------------------------
dt2 = dateshift(metab.south.daystart_dt,'start','day');
metab.south.daystart_dt = dt2;
south_TT = synchronize(south_dailyAvg,south_range,wspd_dailyAvg,par_dailyAvg,metab.south);
south_TT = removevars(south_TT,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg','datetime_local','Tair','light_lux'});
south_TT(south_TT.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
south_TT(south_TT.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Remove timepoints with anomalous GPP and ER values
anomER = find(south_TT.ER_avg > 0);
anomGPP = find(south_TT.GPP_avg < 0);
south_TT([anomER;anomGPP],:) = [];

% Remove row if any explanatory or response variable is NaN
south_TT = rmmissing(south_TT,'DataVariables',{'tidal','salinity','temperature','DOsat','pH','wspd','GPP_avg','ER_avg','NEM_avg'});

% Get month numbers
mo = month(south_TT.datetime_utc);

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);
isSpring = ismember(mo,[3 4 5]);
isSummer = ismember(mo,[6 7 8]);
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

% Add a column for season number
south_TT.season(indWinter) = 1;
south_TT.season(indSpring) = 2;
south_TT.season(indSummer) = 3;
south_TT.season(indFall) = 4;
south_TT.season = categorical(south_TT.season);

% Add a column for site number
south_TT.site = repelem(3,height(south_TT))';
south_TT.site = categorical(south_TT.site);

%====Create a table to store predictor variables for all sites=============
gull_dat.all = timetable2table(gull_TT);
gull_dat.all.datetime_utc = [];

north_dat.all = timetable2table(north_TT);
north_dat.all.datetime_utc = [];

south_dat.all = timetable2table(south_TT);
south_dat.all.datetime_utc = [];

allSites_dat.all = [gull_dat.all;north_dat.all;south_dat.all];

%==========================================================================
% Construct multiple linear regression models for ER, GPP, and NEM
%==========================================================================

%====MLRs across all sites: (1) for each season and (2) across all seasons

%----Original (not scaled) data--------------------------------------------
% Break data into seasons
isWinter = ismember(allSites_dat.all.season,'1');
isSpring = ismember(allSites_dat.all.season,'2');
isSummer = ismember(allSites_dat.all.season,'3');
isFall = ismember(allSites_dat.all.season,'4');

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

allSites_dat.winter = allSites_dat.all(indWinter,:);
allSites_dat.spring = allSites_dat.all(indSpring,:);
allSites_dat.summer = allSites_dat.all(indSummer,:);
allSites_dat.fall = allSites_dat.all(indFall,:);

% Compute models with original data
data = allSites_dat;

% (1) By season
% GPP
allSites_mdl.winter.GPP_avg = fitlm(data.winter,'GPP_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl.spring.GPP_avg = fitlm(data.spring,'GPP_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl.summer.GPP_avg = fitlm(data.summer,'GPP_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl.fall.GPP_avg = fitlm(data.fall,'GPP_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
% ER
allSites_mdl.winter.ER_avg = fitlm(data.winter,'ER_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl.spring.ER_avg = fitlm(data.spring,'ER_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl.summer.ER_avg = fitlm(data.summer,'ER_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl.fall.ER_avg = fitlm(data.fall,'ER_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
% NEM
allSites_mdl.winter.NEM_avg = fitlm(data.winter,'NEM_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl.spring.NEM_avg = fitlm(data.spring,'NEM_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl.summer.NEM_avg = fitlm(data.summer,'NEM_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl.fall.NEM_avg = fitlm(data.fall,'NEM_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');

% (2) Across all seasons
% GPP
allSites_mdl.all.GPP_avg = fitlm(data.all,'GPP_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + season*site'); % Include interactions
% ER
allSites_mdl.all.ER_avg = fitlm(data.all,'ER_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + season*site'); % Include interactions
% NEM
allSites_mdl.all.NEM_avg = fitlm(data.all,'NEM_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + season*site'); % Include interactions

%----Scaled data-----------------------------------------------------------
% Scale response and explanatory variables to center 0 and std 1 (Lowe et al., 2019)
VN = allSites_dat.all.Properties.VariableNames;
colNr = find(strcmp(VN,'NEM_avg'));   % Find column number for NEM
allSites_dat_norm.all = normalize(allSites_dat.all(:,1:colNr)); % Normalize all parameters except "season" and "site"
allSites_dat_norm.all.("season") = allSites_dat.all.season; % Add "season" column back in
allSites_dat_norm.all.("site") = allSites_dat.all.site; % Add "site" column back in

allSites_dat_norm.winter = allSites_dat_norm.all(indWinter,:);
allSites_dat_norm.spring = allSites_dat_norm.all(indSpring,:);
allSites_dat_norm.summer = allSites_dat_norm.all(indSummer,:);
allSites_dat_norm.fall = allSites_dat_norm.all(indFall,:);

% Compute models with normalized data
data_norm = allSites_dat_norm;

% (1) By season
% GPP
allSites_mdl_norm.winter.GPP_avg = fitlm(data_norm.winter,'GPP_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl_norm.spring.GPP_avg = fitlm(data_norm.spring,'GPP_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl_norm.summer.GPP_avg = fitlm(data_norm.summer,'GPP_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl_norm.fall.GPP_avg = fitlm(data_norm.fall,'GPP_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
% ER
allSites_mdl_norm.winter.ER_avg = fitlm(data_norm.winter,'ER_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl_norm.spring.ER_avg = fitlm(data_norm.spring,'ER_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl_norm.summer.ER_avg = fitlm(data_norm.summer,'ER_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl_norm.fall.ER_avg = fitlm(data_norm.fall,'ER_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
% NEM
allSites_mdl_norm.winter.NEM_avg = fitlm(data_norm.winter,'NEM_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl_norm.spring.NEM_avg = fitlm(data_norm.spring,'NEM_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl_norm.summer.NEM_avg = fitlm(data_norm.summer,'NEM_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');
allSites_mdl_norm.fall.NEM_avg = fitlm(data_norm.fall,'NEM_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + site');

% (2) Across all seasons
% GPP
allSites_mdl_norm.all.GPP_avg = fitlm(data_norm.all,'GPP_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + season*site'); % Include interactions
% ERf
allSites_mdl_norm.all.ER_avg = fitlm(data_norm.all,'ER_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + season*site'); % Include interactions
% NEM
allSites_mdl_norm.all.NEM_avg = fitlm(data_norm.all,'NEM_avg ~ tidal + wspd + temperature + salinity + DOsat + pH + season*site'); % Include interactions

%%
%==========================================================================
% Plot results for normalized data
%==========================================================================
data = allSites_dat_norm;
mdl = allSites_mdl_norm;

respNames = {'GPP_avg','ER_avg','NEM_avg'};
predNames = {'salinity','temperature','DOsat','pH','tidal','wspd'};
seasonNames = {'all','winter','spring','summer','fall'};

% Find means of each variable for plotting individual regression lines
for k = 1:length(predNames)
    meanVal.all.(predNames{k}) = mean(data.all.(predNames{k}),'omitmissing');
    meanVal.winter.(predNames{k}) = mean(data.winter.(predNames{k}),'omitmissing');
    meanVal.spring.(predNames{k}) = mean(data.spring.(predNames{k}),'omitmissing');
    meanVal.summer.(predNames{k}) = mean(data.summer.(predNames{k}),'omitmissing');
    meanVal.fall.(predNames{k}) = mean(data.fall.(predNames{k}),'omitmissing');
end
xlblNames = {'Salinity','Temperature','DO % saturation','pH','Tidal range','Wind speed'}; % Normalized, so no units!
ylblNames = {'GPP','ER','NEM'};

%====GPP, ER, NEM vs. Explanatory Variables================================
for i = 1:length(respNames)   % Loop through plots for GPP, ER, and NEM
    fig(i) = figure(i);clf
    fig(i).Position = [1.0000    1.0000  439.2000  700.8000];

    t = tiledlayout(3,2,'Padding','none','TileSpacing','compact');

    % Coefficient estimates for each "season"
    for m = 1:length(seasonNames)
        % Model with all data
        b0.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(1); % Intercept
        b1.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(2); % Salinity
        b2.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(3); % Temperature
        b3.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(4); % DOsat
        b4.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(5); % pH
        b5.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(6); % Tidal
        b6.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(7); % Wind speed
    end

    for j = 1:length(predNames) % Loop through plotting each predictor
        t1 = tiledlayout(t,2,2,'Padding','none','Tilespacing','compact');
        t1.Layout.Tile = j;

        x = data.all.(predNames{j});
        w = data.winter.(predNames{j});
        sp = data.spring.(predNames{j});
        su = data.summer.(predNames{j});
        f = data.fall.(predNames{j});

        % Winter
        ax1 = nexttile(t1);
        % Plot All data in grey
        plot(data.all.(predNames{j}),data.all.(respNames{i}),'.','color',rgb('lightgrey'),'MarkerSize',8)
        hold on
        % Plot regression line for All data in grey
        switch predNames{j}
            case 'salinity'
                plot(x, b0.all + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.all*meanVal.all.(predNames{2}) + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'),'LineWidth',1.5)
            case 'temperature'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'),'LineWidth',1.5)
            case 'DOsat'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'),'LineWidth',1.5)
            case 'pH'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'),'LineWidth',1.5)
            case 'tidal'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'),'LineWidth',1.5)
            case 'wspd'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x,'color',rgb('grey'),'LineWidth',1.5)
        end
        % Plot Winter data in black
        plot(data.winter.(predNames{j}),data.winter.(respNames{i}),'.k','MarkerSize',8)
        % Plot regression line for Winter data in red
        switch predNames{j}
            case 'salinity'
                plot(w, b0.winter + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*w...
                    + b2.winter*meanVal.winter.(predNames{2}) + b3.winter*meanVal.winter.(predNames{3}) + b4.winter*meanVal.winter.(predNames{4})...
                    + b5.winter*meanVal.winter.(predNames{5}),'r','LineWidth',1.5)
            case 'temperature'
                plot(w, b0.winter + b1.winter*meanVal.winter.(predNames{1}) + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*w...
                    + b3.winter*meanVal.winter.(predNames{3}) + b4.winter*meanVal.winter.(predNames{4})...
                    + b5.winter*meanVal.winter.(predNames{5}),'r','LineWidth',1.5)
            case 'DOsat'
                plot(w, b0.winter + b1.winter*meanVal.winter.(predNames{1}) + b2.winter*meanVal.winter.(predNames{2})...
                    + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*w + b4.winter*meanVal.winter.(predNames{4})...
                    + b5.winter*meanVal.winter.(predNames{5}) + b6.winter*meanVal.winter.(predNames{6}),'r','LineWidth',1.5)
            case 'pH'
                plot(w, b0.winter + b1.winter*meanVal.winter.(predNames{1}) + b2.winter*meanVal.winter.(predNames{2})...
                    + b3.winter*meanVal.winter.(predNames{3}) + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*w...
                    + b5.winter*meanVal.winter.(predNames{5}) + b6.winter*meanVal.winter.(predNames{6}),'r','LineWidth',1.5)
            case 'tidal'
                plot(w, b0.winter + b1.winter*meanVal.winter.(predNames{1}) + b2.winter*meanVal.winter.(predNames{2})...
                    + b3.winter*meanVal.winter.(predNames{3}) + b4.winter*meanVal.winter.(predNames{4})...
                    + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*w + b6.winter*meanVal.winter.(predNames{6}),'r','LineWidth',1.5)
            case 'wspd'
                plot(w, b0.winter + b1.winter*meanVal.winter.(predNames{1}) + b2.winter*meanVal.winter.(predNames{2})...
                    + b3.winter*meanVal.winter.(predNames{3}) + b4.winter*meanVal.winter.(predNames{4})...
                    + b5.winter*meanVal.winter.(predNames{5}) + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*w,'r','LineWidth',1.5)
        end
        pbaspect(ax1,[1 1 1]); % Make relative lengths of axes equal
        % make a text object for the title
        xl = get(gca(),'Xlim');
        yl = get(gca(),'Ylim');
        switch i
            case 1
                text(xl(1)+.2,yl(2)+1.3,'Winter','VerticalAlignment','top','HorizontalAlignment','left','FontSize',9)
            case 2
                text(xl(1)+.2,yl(1)-.9,'Winter','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9)
            case 3
                text(xl(1)+.2,yl(1)-1.7,'Winter','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9)
        end
        hAx = gca;
        hAx.LineWidth = 1;

        % Spring
        ax2 = nexttile(t1);
        plot(data.all.(predNames{j}),data.all.(respNames{i}),'.','color',rgb('lightgrey'),'MarkerSize',8)
        hold on
        % Plot regression line for All data in grey
        switch predNames{j}
            case 'salinity'
                plot(x, b0.all + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.all*meanVal.all.(predNames{2}) + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'temperature'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'DOsat'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'pH'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'tidal'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'wspd'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x,'color',rgb('grey'))
        end
        % Plot Spring data in black
        plot(data.spring.(predNames{j}),data.spring.(respNames{i}),'.k','MarkerSize',8)
        % Plot regression line for Spring data in red
        switch predNames{j}
            case 'salinity'
                plot(sp, b0.spring + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*sp...
                    + b2.spring*meanVal.spring.(predNames{2}) + b3.spring*meanVal.spring.(predNames{3}) + b4.spring*meanVal.spring.(predNames{4})...
                    + b5.spring*meanVal.spring.(predNames{5}),'r')
            case 'temperature'
                plot(sp, b0.spring + b1.spring*meanVal.spring.(predNames{1}) + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*sp...
                    + b3.spring*meanVal.spring.(predNames{3}) + b4.spring*meanVal.spring.(predNames{4})...
                    + b5.spring*meanVal.spring.(predNames{5}),'r')
            case 'DOsat'
                plot(sp, b0.spring + b1.spring*meanVal.spring.(predNames{1}) + b2.spring*meanVal.spring.(predNames{2})...
                    + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*sp + b4.spring*meanVal.spring.(predNames{4})...
                    + b5.spring*meanVal.spring.(predNames{5}) + b6.spring*meanVal.spring.(predNames{6}),'r')
            case 'pH'
                plot(sp, b0.spring + b1.spring*meanVal.spring.(predNames{1}) + b2.spring*meanVal.spring.(predNames{2})...
                    + b3.spring*meanVal.spring.(predNames{3}) + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*sp...
                    + b5.spring*meanVal.spring.(predNames{5}) + b6.spring*meanVal.spring.(predNames{6}),'r')
            case 'tidal'
                plot(sp, b0.spring + b1.spring*meanVal.spring.(predNames{1}) + b2.spring*meanVal.spring.(predNames{2})...
                    + b3.spring*meanVal.spring.(predNames{3}) + b4.spring*meanVal.spring.(predNames{4})...
                    + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*sp + b6.spring*meanVal.spring.(predNames{6}),'r')
            case 'wspd'
                plot(sp, b0.spring + b1.spring*meanVal.spring.(predNames{1}) + b2.spring*meanVal.spring.(predNames{2})...
                    + b3.spring*meanVal.spring.(predNames{3}) + b4.spring*meanVal.spring.(predNames{4})...
                    + b5.spring*meanVal.spring.(predNames{5}) + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*sp,'r')
        end
        pbaspect(ax2,[1 1 1]); % Make relative lengths of axes equal
        % make a text object for the title
        xl = get(gca(),'Xlim');
        yl = get(gca(),'Ylim');
        switch i
            case 1
                text(xl(1)+.2,yl(2)+1.3,'Spring','VerticalAlignment','top','HorizontalAlignment','left','FontSize',9)
            case 2
                text(xl(1)+.2,yl(1)-.9,'Spring','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9)
            case 3
                text(xl(1)+.2,yl(1)-1.7,'Spring','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9)
        end
        hAx = gca;
        hAx.LineWidth = 1;

        % Summer
        ax3 = nexttile(t1);
        plot(data.all.(predNames{j}),data.all.(respNames{i}),'.','color',rgb('lightgrey'),'MarkerSize',8)
        hold on
        % Plot regression line for All data in grey
        switch predNames{j}
            case 'salinity'
                plot(x, b0.all + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.all*meanVal.all.(predNames{2}) + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'temperature'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'DOsat'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'pH'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'tidal'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'wspd'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x,'color',rgb('grey'))
        end
        % Plot Summer data in black
        plot(data.summer.(predNames{j}),data.summer.(respNames{i}),'.k','MarkerSize',8)
        % Plot regression line for Summer data in red
        switch predNames{j}
            case 'salinity'
                plot(su, b0.summer + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*su...
                    + b2.summer*meanVal.summer.(predNames{2}) + b3.summer*meanVal.summer.(predNames{3}) + b4.summer*meanVal.summer.(predNames{4})...
                    + b5.summer*meanVal.summer.(predNames{5}),'r')
            case 'temperature'
                plot(su, b0.summer + b1.summer*meanVal.summer.(predNames{1}) + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*su...
                    + b3.summer*meanVal.summer.(predNames{3}) + b4.summer*meanVal.summer.(predNames{4})...
                    + b5.summer*meanVal.summer.(predNames{5}),'r')
            case 'DOsat'
                plot(su, b0.summer + b1.summer*meanVal.summer.(predNames{1}) + b2.summer*meanVal.summer.(predNames{2})...
                    + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*su + b4.summer*meanVal.summer.(predNames{4})...
                    + b5.summer*meanVal.summer.(predNames{5}) + b6.summer*meanVal.summer.(predNames{6}),'r')
            case 'pH'
                plot(su, b0.summer + b1.summer*meanVal.summer.(predNames{1}) + b2.summer*meanVal.summer.(predNames{2})...
                    + b3.summer*meanVal.summer.(predNames{3}) + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*su...
                    + b5.summer*meanVal.summer.(predNames{5}) + b6.summer*meanVal.summer.(predNames{6}),'r')
            case 'tidal'
                plot(su, b0.summer + b1.summer*meanVal.summer.(predNames{1}) + b2.summer*meanVal.summer.(predNames{2})...
                    + b3.summer*meanVal.summer.(predNames{3}) + b4.summer*meanVal.summer.(predNames{4})...
                    + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*su + b6.summer*meanVal.summer.(predNames{6}),'r')
            case 'wspd'
                plot(su, b0.summer + b1.summer*meanVal.summer.(predNames{1}) + b2.summer*meanVal.summer.(predNames{2})...
                    + b3.summer*meanVal.summer.(predNames{3}) + b4.summer*meanVal.summer.(predNames{4})...
                    + b5.summer*meanVal.summer.(predNames{5}) + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*su,'r')
        end
        pbaspect(ax3,[1 1 1]); % Make relative lengths of axes equal
        % make a text object for the title
        xl = get(gca(),'Xlim');
        yl = get(gca(),'Ylim');
        switch i
            case 1
                text(xl(1)+.2,yl(2)+1.3,'Summer','VerticalAlignment','top','HorizontalAlignment','left','FontSize',9)
            case 2
                text(xl(1)+.2,yl(1)-.9,'Summer','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9)
            case 3
                text(xl(1)+.2,yl(1)-1.7,'Summer','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9)
        end
        hAx = gca;
        hAx.LineWidth = 1;
     
        % Fall
        ax4 = nexttile(t1);
        plot(data.all.(predNames{j}),data.all.(respNames{i}),'.','color',rgb('lightgrey'),'MarkerSize',8)
        hold on
        % Plot regression line for All data in grey
        switch predNames{j}
            case 'salinity'
                plot(x, b0.all + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.all*meanVal.all.(predNames{2}) + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'temperature'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'DOsat'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'pH'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'tidal'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'wspd'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x,'color',rgb('grey'))
        end
        % Plot Fall data in black
        plot(data.fall.(predNames{j}),data.fall.(respNames{i}),'.k','MarkerSize',8)
        % Plot regression line for Fall data in red
        switch predNames{j}
            case 'salinity'
                plot(f, b0.fall + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*f...
                    + b2.fall*meanVal.fall.(predNames{2}) + b3.fall*meanVal.fall.(predNames{3}) + b4.fall*meanVal.fall.(predNames{4})...
                    + b5.fall*meanVal.fall.(predNames{5}),'r')
            case 'temperature'
                plot(f, b0.fall + b1.fall*meanVal.fall.(predNames{1}) + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*f...
                    + b3.fall*meanVal.fall.(predNames{3}) + b4.fall*meanVal.fall.(predNames{4})...
                    + b5.fall*meanVal.fall.(predNames{5}),'r')
            case 'DOsat'
                plot(f, b0.fall + b1.fall*meanVal.fall.(predNames{1}) + b2.fall*meanVal.fall.(predNames{2})...
                    + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*f + b4.fall*meanVal.fall.(predNames{4})...
                    + b5.fall*meanVal.fall.(predNames{5}) + b6.fall*meanVal.fall.(predNames{6}),'r')
            case 'pH'
                plot(f, b0.fall + b1.fall*meanVal.fall.(predNames{1}) + b2.fall*meanVal.fall.(predNames{2})...
                    + b3.fall*meanVal.fall.(predNames{3}) + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*f...
                    + b5.fall*meanVal.fall.(predNames{5}) + b6.fall*meanVal.fall.(predNames{6}),'r')
            case 'tidal'
                plot(f, b0.fall + b1.fall*meanVal.fall.(predNames{1}) + b2.fall*meanVal.fall.(predNames{2})...
                    + b3.fall*meanVal.fall.(predNames{3}) + b4.fall*meanVal.fall.(predNames{4})...
                    + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*f + b6.fall*meanVal.fall.(predNames{6}),'r')
            case 'wspd'
                plot(f, b0.fall + b1.fall*meanVal.fall.(predNames{1}) + b2.fall*meanVal.fall.(predNames{2})...
                    + b3.fall*meanVal.fall.(predNames{3}) + b4.fall*meanVal.fall.(predNames{4})...
                    + b5.fall*meanVal.fall.(predNames{5}) + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*f,'r')
        end
        pbaspect(ax4,[1 1 1]); % Make relative lengths of axes equal
        % make a text object for the title
        xl = get(gca(),'Xlim');
        yl = get(gca(),'Ylim');
        switch i
            case 1
                text(xl(1)+.2,yl(2)+1.3,'Fall','VerticalAlignment','top','HorizontalAlignment','left','FontSize',9)
            case 2
                text(xl(1)+.2,yl(1)-.9,'Fall','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9)
            case 3
                text(xl(1)+.2,yl(1)-1.7,'Fall','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9)
        end
        hAx = gca;
        hAx.LineWidth = 1;

        % Control x- and y-tick labels
        ax = t1.Children;
        xlbl = get(ax,'XLabel');
        ylbl = get(ax,'YLabel');
        set([xlbl{:}],'String','');
        set([ylbl{:}],'String','');
        set(ax,'XTickLabel',{})
        set(ax,'YTickLabel',{})
        set([ax(1),ax(2)],'XTickLabelMode','auto','FontSize',8)
        set([ax(2),ax(4)],'YTickLabelMode','auto','FontSize',8)

        xlabel(t1,xlblNames(j),'FontSize',11);
        t1.XLabel.VerticalAlignment = 'middle'; % Control distance of xlabel to xaxis

        switch i
            case 1
                set([ax1 ax2 ax3 ax4],'YLim',[-2.5 6.5])
            case 2
                set([ax1 ax2 ax3 ax4],'YLim',[-6.5 2.5])
            case 3
                set([ax1 ax2 ax3 ax4],'YLim',[-9 3])
        end

        switch j
            case {1 3 5}
                ylabel(t1,ylblNames{i},'FontSize',11)
            otherwise
                ylabel(t1,'')
        end
    end
end

%====Save the plots========================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\stats-analyses\Metabolic Drivers'])
        saveas(fig(1),'GPP_norm.fig')
        exportgraphics(fig(1),'GPP_norm.jpg','Resolution',600)
        saveas(fig(2),'ER_norm.fig')
        exportgraphics(fig(2),'ER_norm.jpg','Resolution',600)
        saveas(fig(3),'NEM_norm.fig')
        exportgraphics(fig(3),'NEM_norm.jpg','Resolution',600)
        disp('Plots saved!')
    case 'No'
        disp('Plots not saved.')
end