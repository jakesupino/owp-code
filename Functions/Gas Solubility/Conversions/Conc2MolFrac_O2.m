%% Script to convert O2 concentration in seawater (umol/L SW) to mol% and wt%

% Uses function O2sol

clear all
close all
clc
set(0,'defaulttextinterpreter','latex')

cd('C:\Users\ejchu\Dropbox\Work_EJC\MATLAB Scripts\Gas Solubility\Conversions')

% INPUTS
S = 0;              % Sample salinity [ppt]
T = 23;              % Sample temperature [degree C]
C_O2 = 350;          % Target concentration [umol/L SW]

% CONSTANTS
MW_O2 = 31.999;      % [g/mol]
MW_N2 = 28.0134;     % [g/mol]

% Atmospheric Composition
% Source: https://www.engineeringtoolbox.com/air-composition-d_212.html
xstar_O2 = 0.20946;
xstar_N2 = 0.78084;

% CALCULATIONS 

Cstar_O2 = O2sol(S,T);                  % Conc. of O2 in eqbm w/ atmos. [umol/kg SW] 

x_O2 = xstar_O2 * C_O2/Cstar_O2;        % Mole fraction of O2 [mol O2/mol mixture]
x_N2 = 1 - x_O2;                        % Mole fraction of N2 [mol O2/mol mixture]

MW_mix = x_O2*MW_O2 + x_N2*MW_N2;       % Molecular weight of mixture [g/mol]

display(['At S = ',num2str(S,3),' ppt and T = 25oC, C_O2 = ',num2str(C_O2,3),' uM corresponds to x_O2 = ',num2str(x_O2*100,2), ' mol%']);
