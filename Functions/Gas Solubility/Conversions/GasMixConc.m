% Script to calculate gas concentrations [umol/kg sw] for a gas standard with known
% composition at a given S and T

clear all
close all
clc

%% Parameters
% Atmospheric mol % (const.)
x0_CO2 = 0.033;
x0_CH4 = 0.000179;
x0_Ar = 0.934;
x0_O2 = 20.946;
x0_N2 = 78.084;

% Sample water
S = 17.42;      % ppt
T = 25.1;       % deg C

% Gas standard mol %
x_CO2 = 0.75;
x_CH4 = 0.25;
x_Ar = 1.3;
x_O2 = 21.02;
x_N2 = 76.68;

%% Calculations
CO2conc = x_CO2/x0_CO2 * CO2sol(S,T);
CH4conc = x_CH4/x0_CH4 * CH4sol(S,T);
Arconc = x_Ar/x0_Ar * Arsol(S,T);
O2conc = x_O2/x0_O2 * O2sol(S,T);
N2conc = x_N2/x0_N2 * N2sol(S,T);

concs = [CO2conc; CH4conc; Arconc; O2conc; N2conc];
