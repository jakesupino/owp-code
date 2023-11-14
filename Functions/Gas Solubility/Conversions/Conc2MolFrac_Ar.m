%% Script to convert Ar concentration in seawater (umol/L SW) to mol% and wt%

% Uses function Arsol

clear all;close all;clc
set(0,'defaulttextinterpreter','latex')
cd('C:\Users\ejchu\Dropbox\Work_EJC\MATLAB Scripts\Gas Solubility\Conversions')

% INPUTS
S = 35;              % Sample salinity [ppt]
T = 22;              % Sample temperature [degree C]
C_Ar = 11;          % Target concentration [umol/L SW]

% CONSTANTS
MW_Ar = 39.948;      % [g/mol]
MW_N2 = 28.0134;     % [g/mol]

% Atmospheric Composition
% Source: https://www.engineeringtoolbox.com/air-composition-d_212.html
xstar_Ar = 0.00934;
xstar_N2 = 0.78084;

% CALCULATIONS 
Cstar_Ar = Arsol(S,T);                  % Conc. of Ar in eqbm w/ atmos. [umol/kg SW] 

x_Ar = xstar_Ar * C_Ar/Cstar_Ar;        % Mole fraction of Ar [mol Ar/mol mixture]
x_N2 = 1 - x_Ar;                        % Mole fraction of N2 [mol N2/mol mixture]

MW_mix = x_Ar*MW_Ar + x_N2*MW_N2;       % Molecular weight of mixture [g/mol]

display(['At S = ',num2str(S,3),' ppt and T = ',num2str(T),'oC, C_Ar = ',num2str(C_Ar,3),' uM corresponds to x_Ar = ',num2str(x_Ar*100,2), ' mol%']);
