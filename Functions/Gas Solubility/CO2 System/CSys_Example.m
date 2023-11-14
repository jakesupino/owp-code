clear all
close all
clc

%% INPUTS
xCO2 = 0.000404;    % Mole fraction of CO2 in the dry mixture
Alk = 2388;         % Total alkalinity (umol/kgSW)

%% Set up function and calculate

pCO2_in = xCO2*1E6;    % pCO2 = xCO2*patm (uatm)

par1type = 1;          % Choose C system parameter #1
par1     = Alk;        % Value of parameter #1
par2type = 4;          % Choose C system parameter #2
par2     = pCO2_in;    % Value of parameter #2
sal      = 36.1;       % Salinity of the sample
tempin   = 25.75;      % Temperature at input conditions
presin   = 10.1325;    % Pressure at input conditions (in dbar)
tempout  = 45;         % Temperature at output conditions 
presout  = 10.1325;    % Pressure at output conditions (in dbar)
sil      = 0;          % Concentration of silicate in the sample (in umol/kg) - 0 if unknown
po4      = 0;          % Concentration of phosphate in the sample (in umol/kg) - 0 if unknown
pHscale  = 1;          % pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c    = 4;          % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    = 1;          % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

%% Do the calculation. See CO2SYS's help for syntax and output format
A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);

%% OUTPUTS
DIC = A(2);         % (umol/kgSW)
pCO2 = A(19);       % (uatm)
fCO2 = A(20);       % (uatm)
conc_CO2 = A(23);   % (umol/kgSW)
