%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dielAnalysis.m
% This script uses the diel oxygen method to calculate GPP, ER, and NEM
% following Beck et al. (2015)
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% January 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
cd([rootpath,'open-water-platform-data\gull\cleaned'])
load gull-cleaned.mat

dat_all = sonde1_cleaned;

cd([rootpath,'physical-data\wind-speed'])
load metData.mat

newTimes = sonde1_cleaned.datetime_utc(1):minutes(10):sonde1_cleaned.datetime_utc(end);
metDat2 = retime(metDat,newTimes,'mean'); % Calculate mean of values in each time bin

dat = synchronize(dat_all,metDat2);

dat = rmmissing(dat,'DataVariables',{'DO_conc','wspd'});

figure(1),clf
plot(dat.datetime_utc,dat.DO_conc,'.')
ylabel('DO conc')

%% STEP 1
%====Calculate DO concentration at equilibrium (DO_sat)=================
DO_sat = O2sol(dat.salinity,dat.temperature);     % [umol kg-1]
% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(dat.salinity,dat.temperature,dat.p); % [kg m-3]
% Convert DO_sat units
DO_sat = DO_sat.*rho_sw/1000;   % [mmol m-3]
% Convert DO_conc units
DO_conc = dat.DO_conc*1000/1000;    % [mmol m-3]
% Calculate the percent oxygen saturation
DO_per_sat = dat.DO_conc./DO_sat*100;   % [%]

%% STEP 2
%====Calculate gas exchange at air-water interface=======================
% Caffrey/Beck (uses "reaeration coefficient", ka)

H = 1.593;      % Mean water depth [m] -- CHANGE THIS!!
rho_a = 1.293;  % Density of air [kg m-3]
T = dat.temperature;
S = dat.salinity;
u = dat.wspd;

% Calculate Vw
for i = 1:length(T)
    Uw(i,1) = swp('m',S(i),T(i));  % Dynamic viscosity [kg m-1 s-1]
end
Vw = Uw./rho_sw;    % Kinematic viscosity [m2 s-1]
kB = 1.3806503E-23; % Boltzmann's constant [m2 kg s-2 K-1]
R0 = 1.72E-10;      % Radius of O2 molecule [m]
Tk = T + 273.15;    % [K]
ka = 1/H * 0.1706 .* (kB*Tk./(4*rho_sw.*Vw.^2.*R0)).^0.5 .* (rho_a./rho_sw).^0.5 .* u.^1.81;

% Air-water flux
D = ka.*(DO_sat - dat.DO_conc); % [mmol m-3 h-1]

%% STEP 3
%====Apply filter========================================================
% Needoba uses HBF
% Caffrey et al. (2014) doesn't bother with filtering
% Beck et al. (2015) uses their R package to remove tidal advection

%% STEP 4
%====Calculate rates=====================================================
dCdt = nan(length(dat.DO_conc),1);
dCdt(2:end,1) = diff(dat.DO_conc) ./ hours(diff(dat.datetime_utc));  % [mmol m-3 h-1]

% Get day/night indices based on lat and lon and Hilary's indexDayNight function
% Or use PAR data for SMIIL analysis
lat = 39.04;
lon = -74.76;
UTCoffset = -5; % Need to change to -4/-5 depending if it's DST or not
time_in = dat.datetime_utc;
tol = 0;
[dayind,nightind] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);

figure,clf;
plot(dat.datetime_utc(dayind),dat.DO_conc(dayind),'b.')
hold on;
plot(dat.datetime_utc(nightind),dat.DO_conc(nightind),'k.');
ylabel('DO Concentration (mmol m^{-3})')
legend('Day','Night')

daystart = nightind(find(diff(nightind)>1)) + 1;
dayend = dayind(find(diff(dayind)>1));
daylength = dat.datetime_utc(dayend(2:end)) - dat.datetime_utc(daystart(1:end));
daylength = hours(daylength);

% Hourly rates of respiration and net production
for i = 1:length(daylength)
    % During night hours, P = 0
    R_hourly(i,1) = mean(dCdt(dayend(i):daystart(i)) - D(dayend(i):daystart(i)),'omitnan'); % Hourly rate of respiration; [mmol m-3 h-1]
    % During day hours, P != 0
    P_hourly(i,1) = mean(dCdt(daystart(i):dayend(i+1)) - D(daystart(i):dayend(i+1)),'omitnan');     % Hourly rate of net production; [mmol m-3 h-1]
end

% Daily rates
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

% Convert volumetric rates to depth-integrated (aereal) estimates
GPP = P_daily * H;      % [mmol O2 m-2 d-1]
ER = R_daily * H;       % [mmol O2 m-2 d-1]
NEM = GPP + ER;         % [mmol O2 m-2 d-1]

% Plot results
figure(1),clf
plot(dat.datetime_utc(daystart),GPP,'.-','MarkerSize',12,'LineWidth',2)
hold on
plot(dat.datetime_utc(daystart),ER,'k.-','MarkerSize',12,'LineWidth',2)
plot(dat.datetime_utc(daystart),NEM,'r.-','MarkerSize',12,'LineWidth',2)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)