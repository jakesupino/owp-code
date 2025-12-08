%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dielAnalysis_salinityTest.m
% This script uses the diel oxygen method to calculate GPP, ER, and NEM
% following Beck et al. (2015).
% The detided BC dissolved oxygen data are run through the diel analysis
% using either the BC- or ERDC-measured salinity.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 1/3/2024
% Last updated: 4/15/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
% Import sonde/physical data (output from prepDielData.m)
site = 'Gull';
sondename = 'BC';

cd([rootpath,'\diel-method\owp-data\initial-qc'])

bc = load([site,'-bc_obs.mat']);       % Load BC data
erdc = load([site,'-erdc_obs.mat']);   % Load ERDC data

% Find indices of deployment changes
dep = fillmissing(bc.dat.deployment,'previous'); % Replace NaNs in BC deployment column with previous deployment #
bc.dat.deployment = dep;
ind_dep = find(diff(bc.dat.deployment) > 0);

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = bc.dat.datetime_utc(end);
red = [0.8500 0.3250 0.0980];   % BC sonde
darkred = rgb('DarkRed');
blue = [0 0.4470 0.7410];       % ERDC sonde
midnightblue = rgb('MidnightBlue');
fontsize = 14;
linewidth = 1;
dotsize = 12;
circlesize = 6;

label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
    'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
    'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
    'Deployment 15','Deployment 16','Deployment 17'};

% Load R weighted regression results
cd([rootpath,'diel-method\R-results\',site,'-',sondename])
wtreg_res = readtable('wtreg_res.csv');
wtreg_res.Properties.DimensionNames{1} = 'datetime_utc';
wtreg_res.DateTimeStamp.TimeZone = "UTC";
wtreg_res = table2timetable(wtreg_res);

% Retime weighted regression data to same datetimes as sonde/physical data
wtreg_res.solar_period = [];
newTimes = bc.dat.datetime_utc(1):minutes(10):bc.dat.datetime_utc(end);
wtreg_res_rt = retime(wtreg_res,newTimes,'mean'); 

%====Define input variables================================================
dt_utc = bc.dat.datetime_utc;
dt_local = bc.dat.datetime_local;

% Water
S_bc = bc.dat.salinity;
S_erdc = erdc.dat.salinity;
T = bc.dat.temperature;

p = bc.dat.depth;  % Pressure in dbar and depth in meters are approx. equal
d = bc.dat.depth; 
H = mean(d,'omitnan') + 0.1524;   % Mean water depth [m] plus 6 inches to account for height of sonde above bottom

% Air
Tair = bc.dat.Tair; % [deg C]
patm = bc.dat.patm; % [hPa]
u = bc.dat.wspd;    % [m/s]

% DO concentration conversions
DO_obs = wtreg_res_rt.DO_obs*1000/32;   % Observed DO concentration [mmol m-3]
DO_nrm = wtreg_res_rt.DO_nrm*1000/32;   % Detided DO concentration [mmol m-3]

cd([rootpath,'\figures\diel-analysis'])

% Plot observed & detided DO concentration and tidal level (c.f. Beck Fig. 6)
fig = figure(1);clf
fig.WindowState = 'maximized';
t = tiledlayout(2,1,'TileSpacing','compact');
ax1 = nexttile;
plot(dt_local,DO_obs,'.-','MarkerSize',6,'Linewidth',1)
hold on
plot(dt_local,DO_nrm,'.-','MarkerSize',6,'Linewidth',1)
legend('Observed','Detided')
ylabel('DO conc. (mmol m^{-3})','FontSize',14)
% title([site,' ',sondename,' Sonde'])
ax2 = nexttile;
plot(dt_local,d,'.-','MarkerSize',6,'Linewidth',1)
xlabel('Local Time')
ylabel('Water depth (m)')
linkaxes([ax1 ax2],'x')

%% Conduct diel analysis with BC detided DO data and BC salinity
DO_conc = DO_nrm;
S = S_bc;

%====STEP 1: Calculate DO concentration at equilibrium (DO_sat)============
DO_sat = O2sol(S,T);              % [umol kg-1]

% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(S,T,p);          % [kg m-3]

% Convert DO_sat units
DO_sat = DO_sat.*rho_sw/1000;     % [mmol m-3]

% Calculate the percent oxygen saturation
DO_per_sat = DO_conc./DO_sat*100; % [%]

%====STEP 2: Calculate gas exchange at air-water interface=================
% Caffrey/Beck (uses "reaeration coefficient", ka)
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
ka = 1/H * 0.1706 .* (kB*Tk./(4*rho_sw.*Vw.^2.*R0)).^0.5 .* (rho_a./rho_sw).^0.5 .* u.^1.81;

% Air-water flux
D = ka.*(DO_sat - DO_conc); % [mmol m-3 h-1]

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

D_bc = D;

% Daily rates
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

% Convert volumetric rates to depth-integrated (areal) estimates
GPP = P_daily * H;      % [mmol O2 m-2 d-1]
ER = R_daily * H;       % [mmol O2 m-2 d-1]
NEM = GPP + ER;         % [mmol O2 m-2 d-1]

metab = table(daystart_dt,dayend_dt,daylength,R_hourly,P_hourly,R_daily,P_daily,GPP,ER,NEM);
metab = table2timetable(metab);

bc.('D') = timetable(dt_utc,D);
bc.('metab') = metab;

%% Conduct diel analysis with BC detided DO data and ERDC salinity
DO_conc = DO_nrm;
S = S_erdc;

S = [S;NaN(length(S_bc)-length(S_erdc),1)];

%====STEP 1: Calculate DO concentration at equilibrium (DO_sat)============
DO_sat = O2sol(S,T);              % [umol kg-1]

% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(S,T,p);          % [kg m-3]

% Convert DO_sat units
DO_sat = DO_sat.*rho_sw/1000;     % [mmol m-3]

% Calculate the percent oxygen saturation
DO_per_sat = DO_conc./DO_sat*100; % [%]

%====STEP 2: Calculate gas exchange at air-water interface=================
% Caffrey/Beck (uses "reaeration coefficient", ka)
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
ka = 1/H * 0.1706 .* (kB*Tk./(4*rho_sw.*Vw.^2.*R0)).^0.5 .* (rho_a./rho_sw).^0.5 .* u.^1.81;

% Air-water flux
D = ka.*(DO_sat - DO_conc); % [mmol m-3 h-1]

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

% Mean hourly rates of respiration and net production
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
GPP = P_daily * H;      % [mmol O2 m-2 d-1]
ER = R_daily * H;       % [mmol O2 m-2 d-1]
NEM = GPP + ER;         % [mmol O2 m-2 d-1]

metab = table(daystart_dt,dayend_dt,daylength,R_hourly,P_hourly,R_daily,P_daily,GPP,ER,NEM);
metab = table2timetable(metab);

erdc.('D') = timetable(dt_utc,D);
erdc.('metab') = metab;

%%
% Stats for Deployment 10
dep10_start = bc.dat.datetime_utc(ind_dep(7));
dep10_end = bc.dat.datetime_utc(ind_dep(8));

% Find nearest matching indices for these datetimes
diffs = abs(bc_diel_dtd.daystart_dt - dep10_start);
[~, ind_start] = min(diffs);
diffs = abs(bc_diel_dtd.daystart_dt - dep10_end);
[~, ind_end] = min(diffs);

bc_meanGPP = mean(bc_diel_dtd.GPP(ind_start:ind_end));
erdc_meanGPP = mean(erdc_diel_dtd.GPP(ind_start:ind_end));

bc_meanER = mean(bc_diel_dtd.ER(ind_start:ind_end));
erdc_meanER = mean(erdc_diel_dtd.ER(ind_start:ind_end));

bc_meanNEM = mean(bc_diel_dtd.NEM(ind_start:ind_end));
erdc_meanNEM = mean(erdc_diel_dtd.NEM(ind_start:ind_end));

table(bc_meanGPP,bc_meanER,bc_meanNEM,erdc_meanGPP,erdc_meanER,erdc_meanNEM)

%%
fig4 = figure(4);clf
fig4.WindowState = 'maximized';
% plot(dt_utc,dCdt,'.k','markersize',circlesize,'DisplayName','dDO/dt')
ylabel('mmol m^{-3} h^{-1}')
plot(dt_utc,bc.D,'.','color',red,'markersize',12,'DisplayName','D (BC Salinity)')
hold on
plot(dt_utc,erdc.D,'.','color',blue,'markersize',12,'DisplayName','D (ERDC Salinity)')
ylabel('D (mmol m^{-3} h^{-1})')
xline([bc.dat.datetime_utc(1); bc.dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend('show')
xlabel('UTC')
title([site,' ',sondename,': Detided DO Concentration'])

%%
fig5 = figure(5);clf
fig5.WindowState = 'maximized';
plot(dt_utc(daystart(2:end)),bc.metab.R_hourly,'.-','color',red,'markersize',12,'DisplayName','R_{hourly} (BC)')
hold on
plot(dt_utc(daystart(2:end)),erdc.metab.R_hourly,'.-','color',blue,'markersize',12,'DisplayName','R_{hourly} (ERDC)')
ylabel('Mean hourly R (mmol m^{-3} h^{-1})')
xline([bc.dat.datetime_utc(1); bc.dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend('show')

fig6 = figure(6);clf
fig6.WindowState = 'maximized';
plot(dt_utc(daystart(2:end)),bc.metab.P_hourly,'.-','color',red,'markersize',12,'DisplayName','P_{hourly} (BC)')
hold on
plot(dt_utc(daystart(2:end)),erdc.metab.P_hourly,'.-','color',blue,'markersize',12,'DisplayName','P_{hourly} (ERDC)')
ylabel('Mean hourly P (mmol m^{-3} h^{-1})')
xline([bc.dat.datetime_utc(1); bc.dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend('show')

%%
fig7 = figure(7);clf
fig7.WindowState = 'maximized';
plot(dt_utc(daystart(2:end)),bc.metab.R_daily,'.-','color',red,'markersize',12,'DisplayName','R_{daily} (BC)')
hold on
plot(dt_utc(daystart(2:end)),erdc.metab.R_daily,'.-','color',blue,'markersize',12,'DisplayName','R_{daily} (ERDC)')
ylabel('Mean hourly R (mmol m^{-3} d^{-1})')
xline([bc.dat.datetime_utc(1); bc.dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend('show')

fig8 = figure(8);clf
fig8.WindowState = 'maximized';
plot(dt_utc(daystart(2:end)),bc.metab.P_daily,'.-','color',red,'markersize',12,'DisplayName','P_{daily} (BC)')
hold on
plot(dt_utc(daystart(2:end)),erdc.metab.P_daily,'.-','color',blue,'markersize',12,'DisplayName','P_{daily} (ERDC)')
ylabel('Daily P (mmol m^{-3} d^{-1})')
xline([bc.dat.datetime_utc(1); bc.dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend('show')
%%
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
tiledlayout(2,1)

t1 = nexttile;
plot(bc.metab.daystart_dt,bc.metab.GPP,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(bc.metab.daystart_dt,bc.metab.ER,'k.-','MarkerSize',12,'LineWidth',1)
plot(bc.metab.daystart_dt,bc.metab.NEM,'r.-','MarkerSize',12,'LineWidth',1)
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
xline([bc.dat.datetime_utc(1); bc.dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,' ',sondename,': MATLAB Results Using Detided Data and BC Salinities'])
grid on

t2 = nexttile;
plot(erdc.metab.daystart_dt,erdc.metab.GPP,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(erdc.metab.daystart_dt,metab.ER,'k.-','MarkerSize',12,'LineWidth',1)
plot(erdc.metab.daystart_dt,metab.NEM,'r.-','MarkerSize',12,'LineWidth',1)
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
xline([erdc.dat.datetime_utc(1); erdc.dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,' ',sondename,': MATLAB Results Using Detided Data and ERDC Salinities'])
grid on

ylim([-400 400])
linkaxes([t1 t2])
