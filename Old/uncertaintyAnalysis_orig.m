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
% Load diel analysis results
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')

% Load final QC'd timeseries data
cd([rootpath,'diel-method\owp-data\final-qc'])
load([site,'_obs.mat'])

% Load results (which include standard deviations) from duplicate sensor comparison
cd([rootpath,'open-water-platform-data\',site,'\cleaned\dupcheck'])
load([site,'-bestGuess.mat'])   % Output from dupCheck scripts

% Load model output from DO linear regression with Winkler data
cd([rootpath,'open-water-platform-data\',site,'\cleaned\validation'])
load([site,'-winkler_mdl.mat'])
DOconc_mdl = mdl;
clear mdl

% Load model output from S linear regression with DIC/TA discrete sample data


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

%====STEP 1: Determine gas exchange at air-water interface=================
% Define errors for each term
% For AquaTroll parameters, use mean std
S_sd = stdev.S; % [PSU]
T_sd = stdev.T; % [degC]
p_sd = stdev.d; % [m] (pressure in [dbar] and depth in [m] are approx. equal)
% Not sure what to set as ERA5 wind speed error
% U10_sd = ; 
% Not sure what to set as barometric pressure error
% slp_sd = ;

% DO conc - take larger of errors
DOconc_SEM = DOconc_mdl.RMSE / sqrt(height(DOconc_mdl.Variables)); % Standard error of the mean from Winkler regression [umol/L]
DOconc_sd = max([stdev.DOconc DOconc_SEM]);   % [umol/L]

% Choice of k parameterization - find average of daily stds between three k params
k_sd(1) = mean(std([R2b.fas.k, W2b.fas.k],0,2),'omitmissing');  % R&H-2008 vs W-2014
k_sd(2) = mean(std([R2b.fas.k, Eb.fas.k],0,2),'omitmissing');   % R&H-2008 vs E-2019
k_sd(3) = mean(std([W2b.fas.k, Eb.fas.k],0,2),'omitmissing');   % W-2014 vs E-2019
% Take largest difference between parameterizations
k_sd = max(k_sd);   % [h-1]

% Define inputs for calculations
counter = 1000;         % Number of simulations
param = 'W14';          % Wanninkhof (2014) parameterization
% DOconc = dat.DOconc/1000;   % Convert from [umol/L] to [mol m-3]
DOconc = dat.DOconc;    % [mmol m-3]
U10 = dat.wspd;         % [m/s]
S = dat.salinity;       % [PSU]
T = dat.temperature;    % [degC]
slp = dat.patm/1013.25; % Convert from [hPa] to [atm]
p = dat.depth;          % [m] (pressure in [dbar] and depth in [m] are approx. equal)
% Mean water column depth, H = d + D [m] (see Collab Lab Notebook, Table 2 for manual measurements for D for each site)
switch site
    case 'Gull'
        H = dat.depth + 0.42;
    case 'North'
        H = dat.depth + 0.47;
    case 'South'
        H = dat.depth + 0.80;
end

% Initialize vectors to hold outputs
DOsat = NaN(length(dat.datetime_utc),counter);
k = NaN(length(dat.datetime_utc),counter);
D = NaN(length(dat.datetime_utc),counter);

% Calculate the solutions for D
% 1000 simulations = 53 seconds
tic
for j = 1:counter
    % Randomly draw an input value for each variable every time the calculation is performed
    DOconc_randErr = randn*DOconc_sd;
    S_randErr = randn*S_sd;
    T_randErr = randn*T_sd;
    p_randErr = randn*p_sd;
    k_randErr = randn*k_sd;

    % Calculate DO saturation concentration
    DOsat = O2sol(S+S_randErr, T+T_randErr);                % [umol kg-1]
    % Use Gibbs Seawater toolbox to calculate seawater density
    rho_sw = gsw_rho(S+S_randErr, T+T_randErr, p+p_randErr);  % [kg m-3]
    % Convert DOsat units
    DOsat = DOsat.*rho_sw/1000;       % [mmol m-3]

    % Calculate gas exchange coefficient, k
    [Fd, k_ms] = fas_Fd(DOconc+DOconc_randErr, U10, S+S_randErr, T+T_randErr, slp, 'O2', param); % Fd: [mol m-2 s-1]; k: [m s-1]
    k = 1./H * 3600 .* k_ms;  % Convert to [h-1]

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
plot(dat.datetime_utc,k,'.','DisplayName','MC average')
hold on
plot(fas.datetime_utc,fas.k,'.','DisplayName','Original diel analysis')
legend('show')
ylabel('k (h^{-1})')

figure,clf
plot(dat.datetime_utc,D_avg,'.','DisplayName','MC average')
hold on
plot(fas.datetime_utc,fas.D,'.','DisplayName','Original diel analysis')
legend('show')
ylabel('D (mmol m^{-3} h^{-1})')

%====STEP 2: Calculate HOURLY rates of nighttime respiration (R) and apparent primary production (P)
% Determine day/night times - Get day/night indices based on lat and lon and Hilary's indexDayNight function
lat = 39.08;
lon = -74.78;
time_in = dat.datetime_utc;
tol = 0;
UTCoffset = 0;  % Input datetime vector is in UTC

[dayind,nightind] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);

% Find the indices for when each day starts and stops
daystart = dayind(find(diff(dayind) > 1) + 1);
dayend = dayind(find(diff(dayind) > 1));

% Length of each day
daylength = dat.datetime_utc(dayend(2:end)) - dat.datetime_utc(daystart(1:end-1) - 1);
daylength = hours(daylength);

% Find the actual datetimes for when each day starts and stops
daystart_dt = dat.datetime_utc(daystart(1:end-1));
dayend_dt = dat.datetime_utc(dayend(2:end));

% LOOKS LIKE THIS PART ISN'T NECESSARY
% % Calculate change in DO over sampling intervals (10 min)
% % dCdt_orig = nan(length(DOconc),1);
% % dCdt_orig(2:end,1) = diff(DOconc) ./ hours(diff(dat.datetime_utc));  % [mmol m-3 h-1]

% % Initialize vector to hold output
% dCdt = NaN(length(DOconc),counter);
% % Calculate the solutions for dC/dt
% % 1000 simulations < 3 seconds
% tic
% for j = 1:counter
%     % Randomly draw an input value for each variable every time the calculation is performed
%     DOconc_randErr = randn*DOconc_sd;
% 
%     dCdt(2:end,j) = (diff(DOconc)+DOconc_randErr) ./ hours(diff(dat.datetime_utc));  % [mmol m-3 h-1]
% end
% toc
% 
% dCdt_med = median(dCdt,2);
% dCdt_avg = mean(dCdt,2);
% dCdt_sd = std(dCdt,0,2);
% 
% % Sanity check plots with original diel analysis results for k and D
% figure,clf
% plot(dat.datetime_utc,dCdt_med,'.','DisplayName','MC average')
% hold on
% plot(fas.datetime_utc,dCdt_orig,'.','DisplayName','Original diel analysis')
% legend('show')
% ylabel('dC/dt ()')

% Calculate change in DO over sampling intervals (10 min)
dCdt = NaN(length(DOconc),1);
dCdt(2:end,1) = diff(DOconc) ./ hours(diff(dat.datetime_utc));  % [mmol m-3 h-1]

% Initialize vectors to hold outputs
R_hourly = NaN(length(daylength),counter);
P_hourly = NaN(length(daylength),counter);

% Calculate the solutions for R_hourly and P_hourly (mean hourly rates)
% 1000 simulations = 103 seconds (4 seconds on 10/10?)
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

R_hourly_med = median(R_hourly,2);
R_hourly_avg = mean(R_hourly,2);
R_hourly_sd = std(R_hourly,0,2);

P_hourly_med = median(P_hourly,2);
P_hourly_avg = mean(P_hourly,2);
P_hourly_sd = std(P_hourly,0,2);

%====STEP 3: Calculate DAILY rates of respiration and gross production=====
R_daily_avg = R_hourly_avg .* 24;                          % Daily rate of respiration; [mmol m-3 d-1]
P_daily_avg = (P_hourly_avg - R_hourly_avg) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

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

fig=figure;clf
fig.WindowState = 'maximized';
errorbar(daystart_dt,GPP_avg,GPP_sd,'.-','Color',rgb('blue'),'MarkerSize',12,'LineWidth',1)
hold on
errorbar(daystart_dt,ER_avg,ER_sd,'.-','Color',rgb('black'),'MarkerSize',12,'LineWidth',1)
errorbar(daystart_dt,NEM_avg,NEM_sd,'.-','Color',rgb('red'),'MarkerSize',12,'LineWidth',1)
