%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dielAnalysis_final.m
% This script uses the diel oxygen method to calculate GPP, ER, and NEM
% following the methods descriptions of Caffrey et al. (2014) and Beck et
% al. (2015), but using the k parameterization of Wanninkhof (2014).
% 
% Computes summary statistics for both the observed and detided data.
%
% Output diel analysis results from the detided data (diel_dtd) include
% anomalous values, but have endpoints around big gaps of data removed.
%
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 1/3/2024
% Last updated: 9/25/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

cd([rootpath,'\diel-method\owp-data\final-qc'])
load([site,'_obs.mat'])   % Load the table

% Load R results
cd([rootpath,'diel-method\R-results\final-qc\',site])
wtreg_res = readtable('wtreg_res.csv');
wtreg_res.Properties.DimensionNames{1} = 'datetime_utc';
wtreg_res.DateTimeStamp.TimeZone = "UTC";
wtreg_res = table2timetable(wtreg_res);

% Retime weighted regression data to same datetimes as sonde/physical data
wtreg_res.solar_period = [];
newTimes = dat.datetime_utc(1):minutes(10):dat.datetime_utc(end);
wtreg_res_rt = retime(wtreg_res,newTimes,'mean');

metab_obs = table2timetable(readtable('metab_obs.csv'));
metab_dtd = table2timetable(readtable('metab_dtd.csv'));

%====Define input variables================================================
dt_utc = wtreg_res_rt.DateTimeStamp;

% Water
S = wtreg_res_rt.Sal;
T = wtreg_res_rt.Temp;
p = wtreg_res_rt.Tide;  % Pressure in dbar and depth in meters are approx. equal
d = wtreg_res_rt.Tide;
% Mean water column depth, H = d + D [m] - see Collab Lab Notebook, Table 2 for manual measurements of D for each site
switch site
    case 'Gull'
        H_time = dat.depth + 0.42;
    case 'North'
        H_time = dat.depth + 0.47;
    case 'South'
        H_time = dat.depth + 0.80;
end
% Air
Tair = wtreg_res_rt.ATemp;  % [deg C]
patm = wtreg_res_rt.BP;     % [hPa]
U10 = wtreg_res_rt.WSpd;    % [m/s]

% DO concentration conversions
DO_obs = wtreg_res_rt.DO_obs*1000/32;   % Observed DO concentration [mmol m-3]
DO_nrm = wtreg_res_rt.DO_nrm*1000/32;   % Detided DO concentration [mmol m-3]

%==========================================================================
%   Conduct diel analysis with observed DO data
%==========================================================================
DO_conc = DO_obs;

%====STEP 1: Determine gas exchange at air-water interface=================
DO_sat = O2sol(S,T);                % [umol kg-1]

% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(S,T,p);            % [kg m-3]c

% Convert DO_sat units
DO_sat = DO_sat.*rho_sw/1000;       % [mmol m-3]

% Calculate the percent oxygen saturation
% DO_per_sat = DO_conc./DO_sat*100;   % [%]

% Calculate gas exchange coefficient, k
C = DO_conc/1000;   % [mol m-3]
slp = patm/1013.25; % [atm]
param = 'W14';      % Wanninkhof (2014) parameterization

[Fd, k] = fas_Fd(C,U10,S,T,slp,'O2',param); % Fd: [mol m-2 s-1]; k: [m s-1]
k = 1./H_time * 3600 .* k;  % Convert to [h-1]

% Calculate air-water diffusive flux
D = k.*(DO_sat - DO_conc); % [mmol m-3 h-1]

% % Sanity check plot of calculated D
% figure,clf
% plot(dat.datetime_utc,Fd./H_time*1000*3600,'.')
% hold on
% plot(dat.datetime_utc,D,'.')

% Make a table to save the air-water exchange terms
fas = table(wtreg_res_rt.DateTimeStamp,k,D,'VariableNames',{'datetime_utc','k','D'});
fas = table2timetable(fas);
fas.Properties.VariableUnits = {'h-1','mmol m-3 h-1'};

%====STEP 2: Calculate HOURLY rates of nighttime respiration (R) and apparent primary production (P)
% Determine day/night times - Get day/night indices based on lat and lon and Hilary's indexDayNight function
lat = 39.08;
lon = -74.78;
time_in = dt_utc;
tol = 0;
UTCoffset = 0;  % Input datetime vector is in UTC

[dayind,nightind] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);

% Find the indices for when each day starts and stops
daystart = dayind(find(diff(dayind) > 1) + 1);
dayend = dayind(find(diff(dayind) > 1));

% Length of each day
daylength = dt_utc(dayend(2:end)) - dt_utc(daystart(1:end-1) - 1);
daylength = hours(daylength);

% Find the actual datetimes for when each day starts and stops
daystart_dt = dt_utc(daystart(1:end-1));
dayend_dt = dt_utc(dayend(2:end));

% Calculate change in DO over sampling intervals (10 min)
dCdt = nan(length(DO_conc),1);
dCdt(2:end,1) = diff(DO_conc) ./ hours(diff(dt_utc));  % [mmol m-3 h-1]

% Mean hourly rates
for i = 1:length(daylength)
    % During night hours, P = 0
    R_hourly(i,1) = mean(dCdt(dayend(i+1):daystart(i+1)) - D(dayend(i+1):daystart(i+1)),'omitnan'); % Mean hourly rate of nighttime respiration; [mmol m-3 h-1]
    % During day hours, P != 0
    P_hourly(i,1) = mean(dCdt(daystart(i):dayend(i+1)) - D(daystart(i):dayend(i+1)),'omitnan');     % Mean hourly rate of apparent/net production; [mmol m-3 h-1]
end

%====STEP 3: Calculate DAILY rates of respiration and gross production=====
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

%====STEP 4: Convert volumetric rates to depth-integrated (areal) estimates
daily_depth = groupsummary(dat,"datetime_utc","day","mean","depth");
switch site
    case 'Gull'
        H_daily = daily_depth.mean_depth + 0.42;
    case 'North'
        H_daily = daily_depth.mean_depth + 0.47;
    case 'South'
        H_daily = daily_depth.mean_depth + 0.80;
end
GPP = P_daily .* H_daily(2:end-1);  % [mmol O2 m-2 d-1]
ER = R_daily .* H_daily(2:end-1);   % [mmol O2 m-2 d-1]
NEM = GPP + ER;                     % [mmol O2 m-2 d-1]

date = dateshift(daystart_dt,'start','day');
date = datetime(date,'TimeZone','UTC');
diel_obs = table(date,daystart_dt,dayend_dt,daylength,R_hourly,P_hourly,R_daily,P_daily,GPP,ER,NEM);
diel_obs = table2timetable(diel_obs);

%   Make summary statistics table
% Remove rows with missing data
diel_obs = rmmissing(diel_obs);

% Detided data
anomGPP = length(find(diel_obs.GPP < 0)) / length(diel_obs.GPP) * 100;
anomER = length(find(diel_obs.ER > 0)) / length(diel_obs.ER) * 100;
meanGPP = mean(diel_obs.GPP);
sdGPP = std(diel_obs.GPP);
meanER = mean(diel_obs.ER);
sdER = std(diel_obs.ER);
summStats.obs = table(meanGPP,sdGPP,anomGPP,meanER,sdER,anomER);

%==========================================================================
%   Conduct diel analysis with detided DO data
%==========================================================================
DO_conc = DO_nrm;

%====STEP 1: Determine gas exchange at air-water interface=================
DO_sat = O2sol(S,T);              % [umol kg-1]

% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(S,T,p);          % [kg m-3]

% Convert DO_sat units
DO_sat = DO_sat.*rho_sw/1000;     % [mmol m-3]

% Calculate the percent oxygen saturation
% DO_per_sat = DO_conc./DO_sat*100; % [%]

% Calculate gas exchange coefficient, k
C = DO_conc/1000;    % [mol m-3]
slp = patm/1013.25; % [atm]
param = 'W14';      % Wanninkhof (2014) parameterization
[Fd, k] = fas_Fd(C,U10,S,T,slp,'O2',param); % Fd: [mol m-2 s-1]; k: [m s-1]
k = 1./H_time * 3600 .* k;  % Convert to [h-1]

% Calculate air-water diffusive flux
D = k.*(DO_sat - DO_conc); % [mmol m-3 h-1]

% Make a table to save the air-water exchange terms
fas = table(wtreg_res_rt.DateTimeStamp,k,D,'VariableNames',{'datetime_utc','k','D'});
fas = table2timetable(fas);
fas.Properties.VariableUnits = {'h-1','mmol m-3 h-1'};

%====STEP 2: % Calculate HOURLY rates of nighttime respiration (R) and apparent primary production (P)
% Determine day/night times - Get day/night indices based on lat and lon and Hilary's indexDayNight function
lat = 39.08;
lon = -74.78;
time_in = dt_utc;
tol = 0;
UTCoffset = 0;  % Input datetime vector is in UTC

[dayind,nightind] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);

% Find the indices for when each day starts and stops
daystart = dayind(find(diff(dayind) > 1) + 1);
dayend = dayind(find(diff(dayind) > 1));

% Find corresponding datetimes
daystart_dt = dt_utc(daystart(1:end-1));
dayend_dt = dt_utc(dayend(2:end));

% Length of each day
daylength = dt_utc(dayend(2:end)) - dt_utc(daystart(1:end-1) - 1);
daylength = hours(daylength);

% Calculate change in DO over sampling intervals (10 min)
dCdt = nan(length(DO_conc),1);
dCdt(2:end,1) = diff(DO_conc) ./ hours(diff(dt_utc));  % [mmol m-3 h-1]

% Mean hourly rates
for i = 1:length(daylength)
    % During night hours, P = 0
    R_hourly(i,1) = mean(dCdt(dayend(i+1):daystart(i+1)) - D(dayend(i+1):daystart(i+1)),'omitnan'); % Hourly rate of respiration; [mmol m-3 h-1]
    % During day hours, P != 0
    P_hourly(i,1) = mean(dCdt(daystart(i):dayend(i+1)) - D(daystart(i):dayend(i+1)),'omitnan');     % Hourly rate of net production; [mmol m-3 h-1]
end

%====STEP 3: Calculate DAILY rates of respiration and gross production=====
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

%====STEP 4: Convert volumetric rates to depth-integrated (areal) estimates
daily_depth = groupsummary(dat,"datetime_utc","day","mean","depth");
switch site
    case 'Gull'
        H_daily = daily_depth.mean_depth + 0.42;
    case 'North'
        H_daily = daily_depth.mean_depth + 0.47;
    case 'South'
        H_daily = daily_depth.mean_depth + 0.80;
end
GPP = P_daily .* H_daily(2:end-1);  % [mmol O2 m-2 d-1]
ER = R_daily .* H_daily(2:end-1);   % [mmol O2 m-2 d-1]
NEM = GPP + ER;                     % [mmol O2 m-2 d-1]

date = dateshift(daystart_dt,'start','day');
date = datetime(date,'TimeZone','UTC');
diel_dtd = table(date,daystart_dt,dayend_dt,daylength,R_hourly,P_hourly,R_daily,P_daily,GPP,ER,NEM);
diel_dtd = table2timetable(diel_dtd);

%   Make summary statistics table
% Remove rows with missing data
diel_dtd = rmmissing(diel_dtd);

% Detided data
anomGPP = length(find(diel_dtd.GPP < 0)) / length(diel_dtd.GPP) * 100;
anomER = length(find(diel_dtd.ER > 0)) / length(diel_dtd.ER) * 100;
meanGPP = mean(diel_dtd.GPP);
sdGPP = std(diel_dtd.GPP);
meanER = mean(diel_dtd.ER);
sdER = std(diel_dtd.ER);
summStats.dtd = table(meanGPP,sdGPP,anomGPP,meanER,sdER,anomER);

%==========================================================================
%   Delete endpoint values of GPP, ER, and NEM around big gaps in data 
%==========================================================================
% (happens especially for North and South -- gaps accompanied by large NEM absolute values)
threshold = duration(days(7));
gap = diff(diel_dtd.date);
idx = find(gap > threshold);
idx_delete = [idx;idx+1];

diel_dtd(idx_delete,:) = [];

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(diel_dtd.date,diel_dtd.GPP,'.-','MarkerSize',12,'LineWidth',1,'DisplayName','GPP (Detided)')
hold on
plot(diel_dtd.date,diel_dtd.ER,'k.-','MarkerSize',12,'LineWidth',1,'DisplayName','ER (Detided)')
plot(diel_dtd.date,diel_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1,'DisplayName','NEM (Detided)')
plot(diel_obs.date,diel_obs.GPP,':','Color',rgb('blue'),'MarkerSize',12,'LineWidth',1,'DisplayName','GPP (Observed)')
plot(diel_obs.date,diel_obs.ER,':','Color',rgb('red'),'MarkerSize',12,'LineWidth',1,'DisplayName','ER (Observed)')
plot(diel_obs.date,diel_obs.NEM,':','Color',rgb('black'),'MarkerSize',12,'LineWidth',1,'DisplayName','NEM (Observed)')
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend('show','FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,': MATLAB Results'])
ylim([-500 500])

%==========================================================================
%   Option to save the figures & results
%==========================================================================
%====Save the figure=======================================================
option = questdlg('Save results figure?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'figures\diel-analysis\matlab-results\',site])
        saveas(fig1,'dielAnalysis_daily.fig')
        saveas(fig1,'dielAnalysis_daily.png')
        disp('Figure saved!')
    case 'No'
        disp('Figure not saved.')
end

%====Save the data=========================================================
option = questdlg('Save diel analysis results?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'diel-method\matlab-results\final-qc\',site])
        save('diel_res.mat','diel_obs','diel_dtd','summStats','fas')
        disp('Files saved!')
    case 'No'
        disp('Files not saved.')
end
