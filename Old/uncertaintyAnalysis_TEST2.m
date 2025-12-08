%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uncertaintyAnalysis.m
% This script...
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 9/23/2024
% Last updated: 10/8/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

%==========================================================================
% Import data
%==========================================================================
% Load final QC'd timeseries data
cd([rootpath,'diel-method\owp-data\final-qc'])
load([site,'_obs.mat'])

% Load diel analysis results
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')

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

% Load results (which include uncertainties) from duplicate sensor comparison
cd([rootpath,'open-water-platform-data\',site,'\cleaned\dupcheck'])
load([site,'-bestGuess.mat'])   % Output from dupCheck scripts
% Define errors for each term
S_err = two_sigma.S; % [PSU]
T_err = two_sigma.T; % [degC]
p_err = two_sigma.d; % [dbar] (pressure in [dbar] and depth in [m] are approx. equal)
d_err = two_sigma.d; % [m]
DOconc_2sig = two_sigma.DOconc;  % [umol/L]

cd([rootpath,'physical-data\final-dataset'])
load('BP&AT.mat')
slp_err = two_sigma.BP; % [hPa]

load('windspeed.mat')
U10_err = two_sigma.wspd; % [m/s]

% Load model output from DO linear regression with Winkler data
cd([rootpath,'open-water-platform-data\',site,'\cleaned\validation'])
load([site,'-winkler_mdl.mat'])
DOconc_mdl = mdl;
clear mdl

% Load diel analysis results from different k parameterizations
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\ro_hunt_wann_vary'])
R2b = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\wann_vary'])
W2b = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\emer_vary'])
Eb = load('diel_res.mat');

%==========================================================================
% Monte Carlo simulations
%==========================================================================
counter = 1000;         % Number of simulations

%====STEP 1: Determine gas exchange at air-water interface=================
% DO conc - take larger of errors
DOconc_SEM = DOconc_mdl.RMSE / sqrt(height(DOconc_mdl.Variables)); % Standard error of the mean from Winkler regression [umol/L]
DOconc_err = max([DOconc_2sig DOconc_SEM]);   % [umol/L]

% Choice of k parameterization - find mean absolute differences between three k params
k_err(1) = mean(abs(diff([R2b.fas.k, W2b.fas.k],1,2)),'omitmissing');  % R&H-2008 vs W-2014
k_err(2) = mean(abs(diff([R2b.fas.k, Eb.fas.k],1,2)),'omitmissing');   % R&H-2008 vs E-2019
k_err(3) = mean(abs(diff([W2b.fas.k, Eb.fas.k],1,2)),'omitmissing');   % W-2014 vs E-2019

% Take largest difference between parameterizations
k_err = max(k_err);   % [h-1]

%====Define input variables================================================
dt_utc = wtreg_res_rt.DateTimeStamp;

% Water
S = wtreg_res_rt.Sal;   % [PSU]
T = wtreg_res_rt.Temp;  % [degC]
p = wtreg_res_rt.Tide;  % [m] (pressure in [dbar] and depth in [m] are approx. equal)
d = wtreg_res_rt.Tide;  % [m]
% Mean water column depth, H = d + D [m] - see Collab Lab Notebook, Table 2 for manual measurements of D for each site
switch site
    case 'Gull'
        H = dat.depth + 0.42;
    case 'North'
        H = dat.depth + 0.47;
    case 'South'
        H = dat.depth + 0.80;
end

% Air
% AT = wtreg_res_rt.ATemp;        % [deg C]
slp = wtreg_res_rt.BP/1013.25; % Convert from [hPa] to [atm]
U10 = wtreg_res_rt.WSpd;       % [m/s]

% DO concentration conversions
DOconc = wtreg_res_rt.DO_nrm*1000/32;   % Detided DO concentration [mmol m-3]

% Define parameterization to use in fas_Fd      
param = 'W14';          % Wanninkhof (2014) parameterization

% Initialize vectors to hold outputs
DOsat = NaN(length(dt_utc),counter);
k = NaN(length(dt_utc),counter);
D = NaN(length(dt_utc),counter);

% Calculate the solutions for D
% 1000 simulations = 45 seconds
tic
for j = 1:counter
    % Randomly draw an input value for each variable every time the calculation is performed
    DOconc_randErr = randn*DOconc_err;
    S_randErr = randn*S_err;
    T_randErr = randn*T_err;
    p_randErr = randn*p_err;
    d_randErr = randn*d_err;
    k_randErr = randn*k_err;
    U10_randErr = randn*U10_err;
    slp_randErr = randn*slp_err;

    % Calculate DO saturation concentration
    DOsat = O2sol(S+S_randErr, T+T_randErr);                % [umol kg-1]
    % Use Gibbs Seawater toolbox to calculate seawater density
    rho_sw = gsw_rho(S+S_randErr, T+T_randErr, p+p_randErr);  % [kg m-3]
    % Convert DOsat units
    DOsat = DOsat.*rho_sw/1000;       % [mmol m-3]

    % Calculate gas exchange coefficient, k
    [Fd, k_ms] = fas_Fd(DOconc+DOconc_randErr, U10+U10_randErr, S+S_randErr, T+T_randErr, slp+slp_randErr, 'O2', param); % Fd: [mol m-2 s-1]; k: [m s-1]
    k = 1./(H+d_randErr) * 3600 .* k_ms;  % Convert to [h-1]

    % Calculate air-water diffusive flux, D
    D(:,j) = (k+k_randErr).*(DOsat - (DOconc+DOconc_randErr)); % [mmol m-3 h-1]
end
toc

% Calculate statistics of solution population
D_med = median(D,2);
D_avg = mean(D,2);
D_sd = std(D,0,2);

% Sanity check plots with original diel analysis results for k and D
figure,clf
plot(dt_utc,k,'.','DisplayName','MC average')
hold on
plot(fas.datetime_utc,fas.k,'.','DisplayName','Original diel analysis')
legend('show')
ylabel('k (h^{-1})')

figure,clf
plot(dt_utc,D_avg,'.','DisplayName','MC average')
hold on
plot(fas.datetime_utc,fas.D,'.','DisplayName','Original diel analysis')
legend('show')
ylabel('D (mmol m^{-3} h^{-1})')

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
dCdt = NaN(length(DOconc),1);
dCdt(2:end,1) = diff(DOconc) ./ hours(diff(dt_utc));  % [mmol m-3 h-1]

% Initialize vectors to hold outputs
R_hourly = NaN(length(daylength),counter);
P_hourly = NaN(length(daylength),counter);

% Calculate the solutions for R_hourly and P_hourly (mean hourly rates)
% 1000 simulations = 4 seconds
tic
for i = 1:length(daylength)
    for j = 1:counter        
        % During night hours, P = 0
        D_randErr = rand*D_sd(dayend(i+1):daystart(i+1));
        R_hourly(i,j) = mean(dCdt(dayend(i+1):daystart(i+1)) - (D_avg(dayend(i+1):daystart(i+1))+D_randErr),'omitnan'); % Mean hourly rate of nighttime respiration; [mmol m-3 h-1]
        
        % During day hours, P != 0
        D_randErr = rand*D_sd(daystart(i):dayend(i+1));
        P_hourly(i,j) = mean(dCdt(daystart(i):dayend(i+1)) - (D_avg(daystart(i):dayend(i+1))+D_randErr),'omitnan');     % Mean hourly rate of apparent/net production; [mmol m-3 h-1]
    end
end
toc

R_hourly_med = median(R_hourly,2);  % [mmol m-3 h-1]
R_hourly_avg = mean(R_hourly,2);    % [mmol m-3 h-1]
R_hourly_sd = std(R_hourly,0,2);    % [mmol m-3 h-1]

P_hourly_med = median(P_hourly,2);  % [mmol m-3 h-1]
P_hourly_avg = mean(P_hourly,2);    % [mmol m-3 h-1]
P_hourly_sd = std(P_hourly,0,2);    % [mmol m-3 h-1]

%====STEP 3: Calculate DAILY rates of respiration and gross production=====
R_daily_avg = R_hourly_avg .* 24;                          % Daily rate of respiration; [mmol m-3 d-1]
P_daily_avg = (P_hourly_avg - R_hourly_avg) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

R_daily_sd = R_hourly_sd .* 24;  % [mmol m-3 h-1]
P_daily_sd = P_hourly_sd .* 24;  % [mmol m-3 h-1]

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
GPP_avg = P_daily_avg .* H_daily(2:end-1);  % [mmol O2 m-2 d-1]
ER_avg = R_daily_avg .* H_daily(2:end-1);   % [mmol O2 m-2 d-1]

% Convert hourly volumetric stds to daily areal units
GPP_sd = P_hourly_sd .* 24 .* H_daily(2:end-1); % [mmol O2 m-2 d-1]
ER_sd = R_hourly_sd .* 24 .* H_daily(2:end-1);  % [mmol O2 m-2 d-1]

% Final Monte Carlo to calculate NEM
NEM = NaN(length(daystart_dt),counter);
% 1000 simulations << 1 second
tic
for j = 1:counter
    GPP_randErr = rand*GPP_sd;
    ER_randErr = rand*ER_sd;

    NEM(:,j) = (GPP_avg+GPP_randErr) + (ER_avg+ER_randErr);  % [mmol O2 m-2 d-1]
end
toc

NEM_med = median(NEM,2);
NEM_avg = mean(NEM,2);
NEM_sd = std(NEM,0,2);

% Make a table of Monte Carlo results
date = dateshift(daystart_dt,'start','day');
date = datetime(date,'TimeZone','UTC');
diel_dtd_MC = table(date,daystart_dt,dayend_dt,daylength,...
    R_hourly_avg,P_hourly_avg,R_daily_avg,P_daily_avg,GPP_avg,ER_avg,NEM_avg,...
    R_hourly_sd,P_hourly_sd,R_daily_sd,P_daily_sd,GPP_sd,ER_sd,NEM_sd);
diel_dtd_MC = table2timetable(diel_dtd_MC);
% Remove rows with missing data
diel_dtd_MC = rmmissing(diel_dtd_MC);
% Delete endpoint values of GPP, ER, and NEM around big gaps in data 
threshold = duration(days(7));
gap = diff(diel_dtd.date);
idx = find(gap > threshold);
idx_delete = [idx;idx+1];
diel_dtd_MC(idx_delete,:) = [];

% Comparison of original diel analysis results and MC averages
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(diel_dtd.date,diel_dtd.GPP,'.-','MarkerSize',12,'LineWidth',1,'DisplayName','GPP (Original)')
hold on
plot(diel_dtd.date,diel_dtd.ER,'k.-','MarkerSize',12,'LineWidth',1,'DisplayName','ER (Original)')
plot(diel_dtd.date,diel_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1,'DisplayName','NEM (Original)')
errorbar(diel_dtd_MC.date,diel_dtd_MC.GPP_avg,diel_dtd_MC.GPP_sd,':','Color',rgb('blue'),'MarkerSize',12,'LineWidth',1,'DisplayName','GPP (Monte Carlo)')
errorbar(diel_dtd_MC.date,diel_dtd_MC.ER_avg,diel_dtd_MC.ER_sd,':','Color',rgb('black'),'MarkerSize',12,'LineWidth',1,'DisplayName','ER (Monte Carlo)')
errorbar(diel_dtd_MC.date,diel_dtd_MC.NEM_avg,diel_dtd_MC.NEM_sd,':','Color',rgb('red'),'MarkerSize',12,'LineWidth',1,'DisplayName','NEM (Monte Carlo)')
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend('show','FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,': MATLAB Results (Detided Data)'])
% ylim([-500 500])

% Original diel analysis results with errorbars for uncertainty from MC analysis
fig2 = figure(2);clf
fig2.WindowState = 'maximized';
errorbar(diel_dtd.date,diel_dtd.GPP,diel_dtd_MC.GPP_sd,'.-b','MarkerSize',12,'LineWidth',1,'DisplayName','GPP (Original)')
hold on
errorbar(diel_dtd.date,diel_dtd.ER,diel_dtd_MC.ER_sd,'.-k','MarkerSize',12,'LineWidth',1,'DisplayName','ER (Original)')
errorbar(diel_dtd.date,diel_dtd.NEM,diel_dtd_MC.NEM_sd,'.-r','MarkerSize',12,'LineWidth',1,'DisplayName','NEM (Original)')
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend('show','FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,': MATLAB Results (Detided Data)'])
ylim([-500 500])
