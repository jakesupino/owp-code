%% Script to convert gas concentrations in seawater (umol/L SW) to mol% and wt%

% Uses functions CH4sol

clear all
close all
clc
set(0,'defaulttextinterpreter','latex')

cd('C:\Users\ejchu\Dropbox\Work_EJC\MATLAB Scripts\Gas Solubility\Conversions')

% INPUTS
% Physical parameters
S = 35;              % Sample salinity [ppt]
T = 20;              % Sample temperature [degree C]
C_CH4 = 2500;         % Target concentration [nmol/L SW]

% CONSTANTS
MW_CH4 = 16.04;      % [g/mol]
MW_N2 = 28.0134;     % [g/mol]

% Atmospheric Composition
% Source: https://www.engineeringtoolbox.com/air-composition-d_212.html
xstar_CH4 = 0.00000179;
xstar_N2 = 0.78084;

% CALCULATIONS 
C_CH4 = C_CH4*10^-3;                 % Use [umol/L SW] for calculations
Cstar_CH4 = CH4sol(S,T);             % Conc. of CH4 in eqbm w/ atmos. [umol/kg SW] 

x_CH4 = xstar_CH4 * C_CH4/Cstar_CH4; % Mole fraction of CH4 [mol CH4/mol mixture]
x_N2 = 1 - x_CH4;                    % Mole fraction of N2 [mol CH4/mol mixture]

MW_mix = x_CH4*MW_CH4 + x_N2*MW_N2;  % Molecular weight of mixture [g/mol]

display(['At S = ',num2str(S,3),' ppt and T = ',num2str(T,3),'^oC C_CH4 = ',num2str(C_CH4*1000,3),' nM corresponds to x_CH4 = ',num2str(x_CH4*100,2), ' mol%']);
