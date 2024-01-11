%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dielAnalysis_beck.m
% This script uses the diel oxygen method to calculate GPP, ER, and NEM
% following Beck et al. (2015) and using their dataset, available at https://github.com/fawda123/WtRegDO
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% January 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import example data================================================
% Data from Beck et al. (2015) on GitHub
varNames = ["date","T","S","DO_conc","ATemp","BP","windspeed","tide"];
varUnits = ["","degC","psu","g m-3","degC","mbar","m s-1","m"];
% dat = readtable('G:\My Drive\Postdoc\Work\SMIIL\diel-method\example-data\beck-data.csv');
dat = readtable('G:\My Drive\Postdoc\Work\SMIIL\diel-method\example-data\sonde1_all-gull.csv');
dat = rmmissing(dat);   % Remove any rows that contain missing data
dat.Properties.VariableNames = varNames;
dat.Properties.VariableUnits = varUnits;

%% STEP 1
%====Calculate DO concentration at equilibrium (DO_sat)=================
DO_sat = O2sol(dat.S,dat.T);     % [umol kg-1]
p = 0;
% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(dat.S,dat.T,p); % [kg m-3]
% Convert DO_sat units
DO_sat = DO_sat.*rho_sw./1000;   % [mmol m-3]
% Calculate the percent oxygen saturation
DO_per_sat = dat.DO_conc./DO_sat*100;   % [%]

%% STEP 2
%====Calculate gas exchange at air-water interface=======================
% Caffrey/Beck (uses "reaeration coefficient", ka)

H = 1.593;   % Mean water depth [m]
rho_a = 1.293;  % Density of air [kg m-3]
T = dat.T;
S = dat.S;
u = dat.windspeed;

% Using constant Dw and Vw as placeholders; see
% https://unisense.com/wp-content/uploads/2021/10/Seawater-Gases-table.pdf
% Dw = 2E-9;   % Oxygen diffusivity [m2 s-1]
% Vw = 1.5E-6; % Kinematic viscosity [m2 s-1]
% ka = 1/H * 0.1706 * (Dw/Vw)^0.5 * (rho_a./rho_sw).^0.5 .* u.^1.81; % Reaeration coefficient [h-1]

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
DO_conc = dat.DO_conc./31.998;  % [mmol m-3]

% Air-water flux
D = ka.*(DO_sat - DO_conc); % [mmol m-3 h-1]

%% STEP 3
%====Apply filter========================================================
% Needoba uses HBF
% Caffrey et al. (2014) doesn't bother with filtering
% Beck et al. (2015) uses their R package to remove tidal advection

%% STEP 4
%====Calculate rates=====================================================
dCdt = nan(length(dat.DO_conc),1);
dCdt(2:end,1) = diff(dat.DO_conc)/0.5;  % [g m-3 h-1]; measurements are at 0.5-h intervals
dCdt = dCdt*1000/31.998;    % [mmol m-3 h-1]

% Get day/night indices based on lat and lon and Hilary's indexDayNight function
% Or use PAR data for SMIIL analysis
lat = 31.39;
lon = -81.28;
UTCoffset = -5;
time_in = dat.date;
tol = 0;
[dayind,nightind] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);
figure(1),clf;
plot(dat.date(dayind),dat.DO_conc(dayind),'b.')
hold on;
plot(dat.date(nightind),dat.DO_conc(nightind),'k.');

% daylength = length(dayind);        % Length of daytime [h]; measurements are at 1-h intervals
daystart = nightind(find(diff(nightind)>1)) + 1;
dayend = dayind(find(diff(dayind)>1));
daylength = dat.date(dayend) - dat.date(daystart(1:end-1));
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
%%
% Plot results
figure(1),clf
plot(dat.date(daystart(2:end)),GPP,'.-','MarkerSize',12,'LineWidth',2)
hold on
plot(dat.date(daystart(2:end)),ER,'k.-','MarkerSize',12,'LineWidth',2)
plot(dat.date(daystart(2:end)),NEM,'r.-','MarkerSize',12,'LineWidth',2)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)