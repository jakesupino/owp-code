%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dielAnalysis_beck.m
% This script uses the diel oxygen method to calculate GPP, ER, and NEM
% following Beck et al. (2015) and using their dataset, available at https://github.com/fawda123/WtRegDO
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% January 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

cd([rootpath,'diel-method\example-data'])

%====Import example data================================================
% Data from Beck et al. (2015) on GitHub
% varNames = ["date","T","S","DO_conc","ATemp","BP","windspeed","tide"];
% varUnits = ["","degC","psu","g m-3","degC","mbar","m s-1","m"];
sapelo = readtable('beck-data.csv');
metab_obs = readtable('beck-metab_obs.csv');
metab_dtd = readtable('beck-metab_dtd.csv');
wtreg_res = readtable('beck-wtreg_res.csv');

metab_obs = table2timetable(metab_obs);
metab_dtd = table2timetable(metab_dtd);

%% STEP 1
%====Calculate DO concentration at equilibrium (DO_sat)=================
S = wtreg_res.Sal;
T = wtreg_res.Temp;
d = wtreg_res.Tide;
p = d;      % Pressure in [dbar] and depth in [m] are approximately equal

DO_sat = O2sol(S,T);     % [umol kg-1]

% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(S,T,p); % [kg m-3]

% Convert DO_sat units
DO_sat = DO_sat.*rho_sw./1000;   % [mmol m-3]

%% STEP 2
%====Calculate gas exchange at air-water interface=======================
% Caffrey/Beck (uses "reaeration coefficient", ka)

% H = mean(d) + 0.5;    % Mean water depth plus 0.5 m to account for height of sonde above bottom [m]
H = 1.593;
rho_a = 1.293;  % Density of air [kg m-3]
u = wtreg_res.WSpd;

% Calculate Vw
for i = 1:length(T)
    Uw(i,1) = swp('m',S(i),T(i));  % Dynamic viscosity [kg m-1 s-1]
end
Vw = Uw./rho_sw; % Kinematic viscosity [m2 s-1]
kB = 1.3806503E-23; % Boltzmann's constant [m2 kg s-2 K-1]
R0 = 1.72E-10;      % Radius of O2 molecule [m]
Tk = T + 273.15;    % [K]
ka = 1/H * 0.1706 .* (kB*Tk./(4*rho_sw.*Vw.^2.*R0)).^0.5 .* (rho_a./rho_sw).^0.5 .* u.^1.81;

% Convert DO concentration from [g m-3] to [mmol m-3]
DO_obs = wtreg_res.DO_obs*1000/32;  % [mmol m-3]
DO_nrm = wtreg_res.DO_nrm*1000/32;  % [mmol m-3]

% INPUT
DO_conc = DO_obs;   % If using observed DO data
% DO_conc = DO_nrm; % If using detided DO data

% Air-water flux
D = ka.*(DO_sat - DO_conc); % [mmol m-3 h-1]

%% STEP 3
%====Calculate rates=====================================================
dCdt = nan(length(DO_conc),1);
dCdt(2:end,1) = diff(DO_conc) ./ hours(diff(wtreg_res.DateTimeStamp));  % [mmol m-3 h-1]

% Get day/night indices based on lat and lon and Hilary's indexDayNight function
% Or use PAR data for SMIIL analysis
lat = 31.39;
lon = -81.28;
UTCoffset = -5;
time_in = wtreg_res.DateTimeStamp;
tol = 0;
[dayind,nightind] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);
% figure(1),clf;
% plot(wtreg_res.DateTimeStamp(dayind),DO_conc(dayind),'b.')
% hold on;
% plot(wtreg_res.DateTimeStamp(nightind),DO_conc(nightind),'k.');

% daylength = length(dayind);        % Length of daytime [h]; measurements are at 1-h intervals
daystart = nightind(find(diff(nightind)>1)) + 1;
dayend = dayind(find(diff(dayind)>1));
daylength = wtreg_res.DateTimeStamp(dayend) - wtreg_res.DateTimeStamp(daystart(1:end-1));
daylength = hours(daylength);

% Hourly rates of respiration and net production
for i = 1:length(daylength)
    % During night hours, P = 0
    R_hourly(i,1) = mean(dCdt(dayend(i):daystart(i+1)) - D(dayend(i):daystart(i+1)),'omitnan'); % Hourly rate of respiration; [mmol m-3 h-1]
    % During day hours, P != 0
    P_hourly(i,1) = mean(dCdt(daystart(i):dayend(i)) - D(daystart(i):dayend(i)),'omitnan');     % Hourly rate of net production; [mmol m-3 h-1]
end

% Daily rates
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

% Convert volumetric rates to depth-integrated (aereal) estimates
GPP = P_daily * H;      % [mmol O2 m-2 d-1]
ER = R_daily * H;       % [mmol O2 m-2 d-1]
NEM = GPP + ER;         % [mmol O2 m-2 d-1]

% Plot results
cd([rootpath,'figures\diel-analysis-figures\Beck Example Data'])
figure(1),clf
plot(wtreg_res.DateTimeStamp(daystart(2:end)),GPP,'.-','MarkerSize',12,'LineWidth',2)
hold on
plot(wtreg_res.DateTimeStamp(daystart(2:end)),ER,'k.-','MarkerSize',12,'LineWidth',2)
plot(wtreg_res.DateTimeStamp(daystart(2:end)),NEM,'r.-','MarkerSize',12,'LineWidth',2)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title('Observed')
% title('Detided')
ylim([-600 600])