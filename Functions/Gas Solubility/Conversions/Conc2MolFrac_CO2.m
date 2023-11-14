%% Script to convert CO2 concentration in seawater (umol/L SW) to mol% and wt%

% Uses function CO2sol.m

clear all
close all
clc
set(0,'defaulttextinterpreter','latex')

cd('C:\Users\Emily\Dropbox\Work_EJC\MATLAB Scripts\Gas Solubility\Conversions')
%% INPUTS
% Physical parameters
S = 30;              % Sample salinity [ppt]
T = 25;              % Sample temperature [degree C]
C_CO2 = 1000;          % Target concentration [umol/L SW]

%% CONSTANTS
MW_CO2 = 44.01;      % [g/mol]
MW_N2 = 28.0134;     % [g/mol]

%% Atmospheric Composition
% Source: https://www.engineeringtoolbox.com/air-composition-d_212.html
xstar_CO2 = 0.00033;
xstar_N2 = 0.78084;

%% CALCULATIONS 

Cstar_CO2 = CO2sol(S,T);             % Conc. of CO2 in eqbm w/ atmos. [umol/kg SW] 

x_CO2 = xstar_CO2 * C_CO2/Cstar_CO2; % Mole fraction of CO2 [mol CO2/mol mixture]
x_N2 = 1 - x_CO2;                    % Mole fraction of N2 [mol CO2/mol mixture]

MW_mix = x_CO2*MW_CO2 + x_N2*MW_N2;  % Molecular weight of mixture [g/mol]

display(['C_CO2 = ',num2str(C_CO2,3),' uM corresponds to x_CO2 = ',num2str(x_CO2*100,2), ' mol%']);
