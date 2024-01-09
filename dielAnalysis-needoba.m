%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dielAnalysis.m
% This script uses the diel oxygen method to calculate GPP, ER, and NEM
% using example dataset for Columbia River Estuary from Needoba et al. (2012).
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% 12/21/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import example data================================================
% Table 2 from Needoba
varNames = ["date","T","S","DO_conc","DO_sat","windspeed","watervel","dDO_hbf","F_O2","BDO_d","BDO_n"];
varUnits = ["","degC","psu","mmol m-3","mmol m-3","m s-1","m s-1","mmol m-2","mmol m-2","mmol m-2","mmol m-2"];
dat = readtable('G:\My Drive\Postdoc\Work\SMIIL\diel-method\example-data\needoba-table2.csv');
dat.Properties.VariableNames = varNames;
dat.Properties.VariableUnits = varUnits;
dat.S = zeros(length(dat.S),1); % S < 0.1 in this example

%% STEP 1
%====Calculate DO concentration at equilibrium (DO_sat)====================
DO_sat = O2sol(dat.S,dat.T);     % [umol kg-1]
p = 0;
% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(dat.S,dat.T,p); % [kg m-3]
% Convert DO_sat units
DO_sat = DO_sat.*rho_sw./1000;   % [mmol m-3]
% Calculate the percent oxygen saturation
DO_per_sat = dat.DO_conc./DO_sat*100;   % [%]

%% STEP 2
%====Calculate gas exchange at air-water interface======================
% Method 1: Needoba Eq. 4.8 ("air-water diffusion flux", F_O2)
% According to the box in Sect 3.2.3, the CRE experiences tidal flows and is
% sufficiently wide, so both k_wind and k_flow must be used to estimate piston velocity

% Calculate k_wind
T = dat.T;
S = dat.S;
u = dat.windspeed;
Sc0 = 1800.6 - 120.1.*T + 3.7818.*T.^2 - 0.047608.*T.^3; % Schmidt number for O2 in freshwater
Sc = Sc0.*(1 + 3.14E-3.*S);   % Schmidt number corrected for salinity 
k_wind = 0.31*u.^2.*(Sc./660).^-0.5;    % [cm h-1]
k_wind = k_wind/100;          % Convert to [m h-1]

% Calculate k_flow
U = dat.watervel;
Tk = T + 273.15;    % Temperature in [K]
h = 7;  % [m] Average depth of CRE; see box in Sect 3.2.3 
x = 2.26;   % Association factor of water
M = 18; % Molar weight of water
nu = (2.414E-2)*10.^(247.8./(Tk - 140)); % Dynamic viscosity of water [centipoise]
V = 25.6;   % Molar volume at the boiling point
D = 7.4E-8*((x*M)^0.5*Tk)./(nu*V^0.6);
k_flow = (U.*D./h).^0.5;  % [cm s-1]
k_flow = k_flow*36;   % [m h-1]

% Piston velocity
v_O2 = k_wind + k_flow; % [m h-1]

% Air-water flux
F_O2 = -v_O2.*(dat.DO_conc - DO_sat);    % [mmol m-2 h-1]  

%% STEP 3
%====Apply filter========================================================
% Needoba uses HBF
% Caffrey et al. (2014) doesn't bother with filtering
% Beck et al. (2015) uses their R package to remove tidal advection

%% STEP 4
%====Calculate rates=====================================================

% Needoba method ("biological oxygen change")
% BDO_t(1,1) = NaN;
% for i = 2:length(dat.DO_conc)
%     BDO_t(i,1) = (dat.DO_conc(i) - dat.DO_conc(i-1)) * h - F_O2(i);
% end

for i = 1:length(dat.DO_conc)-1
    BDO_t(i,1) = (dat.DO_conc(i+1) - dat.DO_conc(i)) * h - F_O2(i);
end
BDO_t(end+1,1) = NaN;

% Check results
figure(2),clf;bar(dat.BDO_d)
hold on;bar(dat.BDO_n)
plot(BDO_t,'k','linewidth',1)

ind_day = find(~isnan(dat.BDO_d));  % Get day/night indices from Needoba's table (where BDO is paritioned into day/night values)
ind_night = find(~isnan(dat.BDO_n));
daylength = length(ind_day)*1;      % Length of daytime [h]
nightlength = length(ind_night)*1;  % Length of daytime [h]

% Total daily net ecosystem metabolism
NEM = sum(BDO_t,'omitnan');   % [mmol O2 m-2 d-1]

% Net ecoystem production = Total BDO during photoperiod
NEP = sum(BDO_t(ind_day),'omitnan');   % [mmol O2 m-2 d-1]; aka daily production rate

% Night respiration = Total BDO during dark period
NR = sum(BDO_t(ind_night),'omitnan'); % [mmol O2 m-2]

% Hourly respiration rate
R_hourly = sum(BDO_t(ind_night),'omitnan')/nightlength;   % [mmol O2 m-2 h-1]
ER = R_hourly*24;  % [mmol O2 m-2 d-1]

% Gross primary production
GPP = NEP + abs(R_hourly)*daylength;    % [mmol O2 m-2 d-1]
