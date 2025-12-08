%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dielAnalysis_compare_k.m
% This script uses the diel oxygen method to calculate GPP, ER, and NEM
% using different k parameterizations.
%
% Plot results in "plotDielResults_compare_k.m"
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 9/5/2024
% Last updated: 9/11/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

% Load the sonde data
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

metab_obs.Date.TimeZone = 'UTC';
metab_dtd.Date.TimeZone = 'UTC';

%====Define input variables================================================
dt_utc = wtreg_res_rt.DateTimeStamp;
dt_local = dt_utc;
dt_local.TimeZone = 'America/New_York';

% Water
S = wtreg_res_rt.Sal;
T = wtreg_res_rt.Temp;
% p = dat.p;
p = wtreg_res_rt.Tide;  % Pressure in dbar and depth in meters are approx. equal
d = wtreg_res_rt.Tide;
% Mean water column depth, H = d + D [m]
% See Collab Lab Notebook - Table 2 for manual measurements for D for each site
switch site
    case 'Gull'
        H = mean(dat.depth,'omitnan') + 0.42;
        H_time = dat.depth + 0.42;
    case 'North'
        H = mean(dat.depth,'omitnan') + 0.47;
        H_time = dat.depth + 0.47;
    case 'South'
        H = mean(dat.depth,'omitnan') + 0.80;
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
%   Conduct diel analysis with detided DO data
%==========================================================================
DO_conc = DO_nrm;

%====STEP 1: Calculate DO concentration at equilibrium (DO_sat)============
DO_sat = O2sol(S,T);              % [umol kg-1]

% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(S,T,p);          % [kg m-3]

% Convert DO_sat units
DO_sat = DO_sat.*rho_sw/1000;     % [mmol m-3]

% Calculate the percent oxygen saturation
DO_per_sat = DO_conc./DO_sat*100; % [%]

%====STEP 2: Calculate gas exchange at air-water interface=================

% Try different k parameterizations - see the cited sources

%----Parameterization 1: From Ro & Hunt (2006)-----------------------------
Tkair = Tair + 273.15;                  % Absolute air temperature [K]
R_specific = 287.0500676;               % Specific gas constant for dry air [J kg-1 K-1]
rho_a = patm*100 ./ (R_specific*Tkair); % Density of dry air [kg m-3]
% Calculate Vw (kinematic viscosity of sw)
for i = 1:length(T)
    Uw(i,1) = swp('m',S(i),T(i));  % Dynamic viscosity of sw [kg m-1 s-1]
end
Vw = Uw./rho_sw;    % Kinematic viscosity [m2 s-1]
kB = 1.3806503E-23; % Boltzmann's constant [m2 kg s-2 K-1]
R0 = 1.72E-10;      % Radius of O2 molecule [m]
Tk = T + 273.15;    % Absolute temperature [K]

% 0) Thebault et al. (2008) implementation using Caffrey et al. (2014) formula with typo
k_R0 = 1/H * 0.1706 .* (kB*Tk./(4*rho_sw.*Vw.^2.*R0)).^0.5 .* (rho_a./rho_sw).^0.5 .* U10.^1.81; % [h-1]

% 1a) Correct Thebault et al. (2008) implementation and const H
k_R1a = 1/H * 1.706 .* (kB*Tk./(4*rho_sw.*Vw.^2.*R0)).^0.5 .* (rho_a./rho_sw).^0.5 .* U10.^1.81; % [h-1]

% 1b) Correct Thebault et al. (2008) implementation and time-varying H
k_R1b = 1./H_time .* 1.706 .* (kB*Tk./(4*rho_sw.*Vw.^2.*R0)).^0.5 .* (rho_a./rho_sw).^0.5 .* U10.^1.81; % [h-1]

% 2a) With Sc formulation from Wanninkhof (2014) and const H
Sc = 1920.4 - 135.6*T + 5.2122*T.^2 - 0.10939*T.^3 + 0.00093777*T.^4; % Schmidt number for O2 for seawater at 35‰ (see Table 1 in W2014)
k_R2a = 1/H * 1.706 .* Sc.^-0.5 .* (rho_a./rho_sw).^0.5 .* U10.^1.81;

% 2b) With Sc formulation from Wanninkhof (2014) and time-varying H
Sc = 1920.4 - 135.6*T + 5.2122*T.^2 - 0.10939*T.^3 + 0.00093777*T.^4; % Schmidt number for O2 for seawater at 35‰ (see Table 1 in W2014)
k_R2b = 1./H_time * 1.706 .* Sc.^-0.5 .* (rho_a./rho_sw).^0.5 .* U10.^1.81;

%----Parameterization 2: From Wanninkhof (2014)----------------------------
% 1) My code
k = 0.251 * U10.^2 .* (Sc/660).^-0.5;    % [cm h-1]
k_W1 = 1/H * k/100;   % Convert to [h-1]

% 2) Via gas_toolbox
C = DO_nrm/1000;    % [mol m-3]
slp = patm/1013.25; % [atm]
param = 'W14';      % Wanninkhof (2014) parameterization
[Fd, k] = fas_Fd(C,U10,S,T,slp,'O2',param); % Fd: [mol m-2 s-1]; k: [m s-1]

% 2a) With const H
k_W2a = 1/H * 3600 * k;  % Convert to [h-1]

% 2b) With time-varying H
k_W2b = 1./H_time * 3600 .* k;  % Convert to [h-1]

%----Parameterization 3: Emerson et al. (2019) via gas_toolbox-------------
[Fd,Fc,Fp,Deq,k] = fas_E19(C,U10,S,T,slp,'O2'); % Fluxes: [mol m-2 s-1]; k: [m s-1]

% a) With const H
k_Ea = 1/H * 3600 * k;    % Convert to [h-1]

% a) With time-varying H
k_Eb = 1./H_time * 3600 .* k;    % Convert to [h-1]

%====Plot results to compare k parameterizations===========================
cd([rootpath,'\figures\diel-analysis\sensitivity-analysis'])

dl = datetime('12-Jun-2023','TimeZone','UTC');
dr = datetime('01-Jul-2023','TimeZone','UTC');

%----Initial sanity check figure: Compare wrong and right Caffrey formulas, and my code vs. gas_toolbox for Wanninkhof (2019)
fig = figure;clf
fig.WindowState = 'maximized';
yyaxis left
plot(wtreg_res_rt.DateTimeStamp,wtreg_res_rt.WSpd,':','color',rgb('sandybrown'),'HandleVisibility','off')
ylabel('Wind speed (m/s)')
ax = gca;
ax.YColor = rgb('sandybrown');
xlim([dl dr])

yyaxis right
plot(wtreg_res_rt.DateTimeStamp,k_R0,'-','color',rgb('mediumblue'),'DisplayName','Ro & Hunt 2006 (Thebault 2008 implementation w/ wrong formula)')
hold on
plot(wtreg_res_rt.DateTimeStamp,k_R1a,':','color',rgb('royalblue'),'DisplayName','Ro & Hunt 2006 (Thebault 2008 implementation w/ right formula)')
plot(wtreg_res_rt.DateTimeStamp,k_R2a,':','color',rgb('skyblue'),'DisplayName','Ro & Hunt 2006 (with Sc from Wanninkhof 2014)')
plot(wtreg_res_rt.DateTimeStamp,k_W1,'-','color',rgb('indigo'),'DisplayName','Wanninkhof 2014 (my code)')
plot(wtreg_res_rt.DateTimeStamp,k_W2a,':','color',rgb('mediumorchid'),'DisplayName','Wanninkhof 2014 (gas toolbox)')
ylabel('k (h^{-1})')
legend('show','location','best')
ax = gca;
ax.YColor = rgb('royalblue');
title('Sanity check of k formulas')
xlim([dl dr])

%----Compare select k parameterizations------------------------------------
fig = figure;clf
fig.WindowState = 'maximized';
yyaxis left
plot(wtreg_res_rt.DateTimeStamp,wtreg_res_rt.WSpd,':','color',rgb('sandybrown'),'HandleVisibility','off')
ylabel('Wind speed (m/s)')
ax = gca;
ax.YColor = rgb('sandybrown');
xlim([dl dr])

yyaxis right
plot(wtreg_res_rt.DateTimeStamp,k_R1a,'-','color',rgb('mediumblue'),'DisplayName','Ro & Hunt 2006 (Thebault 2008 implementation w/ const H)')
hold on
plot(wtreg_res_rt.DateTimeStamp,k_R1b,':','color',rgb('darkcyan'),'DisplayName','Ro & Hunt 2006 (Thebault 2008 implementation w/ time-varying H)')
plot(wtreg_res_rt.DateTimeStamp,k_R2a,'-','color',rgb('skyblue'),'DisplayName','Ro & Hunt 2006 (with Sc from Wanninkhof 2014) w/ const H')
plot(wtreg_res_rt.DateTimeStamp,k_W2a,'-','color',rgb('mediumorchid'),'DisplayName','Wanninkhof 2014 (gas toolbox) w/ const H')
plot(wtreg_res_rt.DateTimeStamp,k_Ea,'-','color',rgb('dodgerblue'),'DisplayName','Emerson et al. 2019 (gas toolbox) w/ const H')
ylabel('k (h^{-1})')
legend('show','location','best')
ax = gca;
ax.YColor = rgb('royalblue');
title('Compare k parameterizations')
xlim([dl dr])

%====Pick a k parameterization and continue on with diel analysis==========
msg = 'Pick a k parameterization to conduct diel analysis';
opts = ["R&H-06 (T-08 w/ const H)", "R&H-06 (T-08 w/ time-varying H)",...
    "R&H-06 (Sc from W-14 w/ const H)", "R&H-06 (Sc from W-14 w/ time-varying H)",...
    "W-14 (w/ const H)", "W-14 (w/ time-varying H)"...
    "E-19 (w/ const H)", "E-19 (w/ time-varying H)"];
choice = menu(msg,opts);

switch choice
    case 1
        k = k_R1a; 
    case 2
        k = k_R1b;
    case 3
        k = k_R2a;
    case 4
        k = k_R2b;
    case 5
        k = k_W2a;
    case 6
        k = k_W2b;
    case 7
        k = k_Ea;
    case 8
        k = k_Eb;
end

% Calculate air-water diffusive flux
D = k.*(DO_sat - DO_conc); % [mmol m-3 h-1]

% Make a table to save the air-water exchange terms
fas = table(wtreg_res_rt.DateTimeStamp,k,D,'VariableNames',{'datetime_utc','k','D'});
fas = table2timetable(fas);
fas.Properties.VariableUnits = {'h-1','mmol m-3 h-1'};

%====STEP 3: Determine day/night times=====================================
% Get day/night indices based on lat and lon and Hilary's indexDayNight function
lat = 39.08;
lon = -74.78;
time_in = dt_utc;
tol = 0;
UTCoffset = 0;  % Input datetime vector is in UTC

[dayind,nightind] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);

% Manually define DST start and end times
start1 = datenum('03/14/2021 02:00','mm/dd/yyyy HH:MM');
end1 = datenum('11/07/2021 02:00','mm/dd/yyyy HH:MM');
start2 = datenum('03/13/2022 02:00','mm/dd/yyyy HH:MM');
end2 = datenum('11/6/2022 02:00','mm/dd/yyyy HH:MM');
start3 = datenum('03/12/2023 02:00','mm/dd/yyyy HH:MM');
end3 = datenum('11/5/2023 02:00','mm/dd/yyyy HH:MM');

% For plotting, create datetimes for when EDT begins/ends each year during data record
end1_dt = datetime(end1,'ConvertFrom','datenum','TimeZone','UTC');
start2_dt = datetime(start2,'ConvertFrom','datenum','TimeZone','UTC');
end2_dt = datetime(end2,'ConvertFrom','datenum','TimeZone','UTC');
start3_dt = datetime(start3,'ConvertFrom','datenum','TimeZone','UTC');
end3_dt = datetime(end3,'ConvertFrom','datenum','TimeZone','UTC');

% Find the indices for when each day starts and stops
daystart = dayind(find(diff(dayind) > 1) + 1);
dayend = dayind(find(diff(dayind) > 1));

% Find corresponding datetimes
daystart_dt = dt_utc(daystart(1:end-1));
dayend_dt = dt_utc(dayend(2:end));

% Length of each day
daylength = dt_utc(dayend(2:end)) - dt_utc(daystart(1:end-1) - 1);
daylength = hours(daylength);

%====STEP 4: Calculate rates===============================================
dCdt = nan(length(DO_conc),1);
dCdt(2:end,1) = diff(DO_conc) ./ hours(diff(dt_utc));  % [mmol m-3 h-1]

% Hourly rates of respiration and net production
for i = 1:length(daylength)
    % During night hours, P = 0
    R_hourly(i,1) = mean(dCdt(dayend(i+1):daystart(i+1)) - D(dayend(i+1):daystart(i+1)),'omitnan'); % Hourly rate of respiration; [mmol m-3 h-1]
    % During day hours, P != 0
    P_hourly(i,1) = mean(dCdt(daystart(i):dayend(i+1)) - D(daystart(i):dayend(i+1)),'omitnan');     % Hourly rate of net production; [mmol m-3 h-1]
end

% Daily rates
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

% Convert volumetric rates to depth-integrated (areal) estimates
% If statement to use constant mean H or time-varying H depending on parameterization
if choice == 2 || 4 || 6 || 8
    daily_depth = groupsummary(dat,"datetime_utc","day","mean","depth");
    switch site
        case 'Gull'
            H_daily = daily_depth.mean_depth + 0.42;
        case 'North'
            H_daily = daily_depth.mean_depth + 0.47;
        case 'South'
            H_daily = daily_depth.mean_depth + 0.80;
    end
    GPP = P_daily .* H_daily(2:end-1);
    ER = R_daily .* H_daily(2:end-1);
    NEM = GPP + ER;
else
    GPP = P_daily * H;      % [mmol O2 m-2 d-1]
    ER = R_daily * H;       % [mmol O2 m-2 d-1]
    NEM = GPP + ER;         % [mmol O2 m-2 d-1]
end
date = dateshift(daystart_dt,'start','day');
date = datetime(date,'TimeZone','UTC');
diel_dtd = table(date,daystart_dt,dayend_dt,daylength,R_hourly,P_hourly,R_daily,P_daily,GPP,ER,NEM);
diel_dtd = table2timetable(diel_dtd);

%==========================================================================
%   Make summary statistics tables
%==========================================================================
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
%   Option to save the sensitivity analysis results
%==========================================================================
%====Save the data=========================================================
option = questdlg('Save diel analysis results?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        switch choice
            case 1
                cd([rootpath,'diel-method\sensitivity-analysis\',site,'\ro_hunt_theb_const'])
            case 2
                cd([rootpath,'diel-method\sensitivity-analysis\',site,'\ro_hunt_theb_vary'])
            case 3
                cd([rootpath,'diel-method\sensitivity-analysis\',site,'\ro_hunt_wann_const'])
            case 4
                cd([rootpath,'diel-method\sensitivity-analysis\',site,'\ro_hunt_wann_vary'])
            case 5
                cd([rootpath,'diel-method\sensitivity-analysis\',site,'\wann_const'])
            case 6
                cd([rootpath,'diel-method\sensitivity-analysis\',site,'\wann_vary'])
            case 7
                cd([rootpath,'diel-method\sensitivity-analysis\',site,'\emer_const'])
            case 8
                cd([rootpath,'diel-method\sensitivity-analysis\',site,'\emer_vary'])
        end
        save('diel_res.mat','diel_dtd','fas','summStats')
        disp('Files saved!')
    case 'No'
        disp('Files not saved.')
end
