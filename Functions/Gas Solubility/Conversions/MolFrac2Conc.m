% Script to convert partial pressures (mole fraction) of each gas in a 
% dry calibration mixture to concentrations in seawater (umol/kgSW) 

% Uses functions Arsol, O2sol, N2sol, N2Osol, CH4sol, and CO2SYS

clear all;close all;clc
cd('C:\Users\ejchu\Dropbox\Work_EJC\MATLAB Scripts\Gas Solubility\Conversions')
% INPUTS

% Physical parameters
%S = 36.1;        % Sample salinity [ppt]
%T = 25.75;       % Sample temperature [degree C]
S=35;
T=20;

% cd('C:\Users\Emily\Dropbox\Work_EJC\Lab Experiments\New Setup\GasRatios\MolFrac')
% [xCO2,xCH4,xAr,xO2,xN2] = importfile('Mix4.txt');

% Mole fractions in dry gas
xO2 = 0; 
xAr = 0; 
xCO2 = 0; 
xCH4 = 0.0025/100; 
xN2 = 1 - (xO2 + xAr + xCO2 + xCH4);   % N2 is usually the balance gas

% CONSTANTS
% Atmospheric composition
% Source: http://ossfoundation.us/projects/environment/global-warming/atmospheric-composition
amb_O2 = 0.2095;
amb_N2 = 0.78084;
amb_Ar = 0.0093;
amb_CO2 = 0.000350;
amb_CH4 = 0.00000179;
amb_N2O = 0.000000325;

% CALCULATIONS
% Original formulas assume the concs. are atmospheric, so need to SCALE them linearly
% IGNORE OUTPUT IN COMMAND WINDOW -- THESE ARE THE UNSCALED CONCENTRATIONS!!!

conc_O2 = xO2/amb_O2*O2sol(S,T);           % [umol/kgSW]
conc_N2 = xN2/amb_N2*N2sol(S,T);           % [umol/kgSW]
conc_Ar = xAr/amb_Ar*Arsol(S,T);           % [umol/kgSW]
conc_CO2 = xCO2/amb_CO2*CO2sol(S,T);       % [umol/kgSW]    
conc_CH4 = xCH4/amb_CH4*CH4sol(S,T);       % [umol/kgSW]
%conc_N2O = xN2O/amb_N2O*N2Osol(S,T);

N2Ar = conc_N2/conc_Ar;
O2Ar = conc_O2/conc_Ar;

%% WRITE OUTPUT TO .TXT FILE
cd('C:\Users\Emily\Dropbox\Work_EJC\Lab Experiments\New Setup\GasRatios\Concentrations')

output = [conc_CO2 conc_CH4 conc_Ar conc_O2 conc_N2 N2Ar O2Ar];
fileID = fopen('Atmos.txt','w');
fprintf(fileID,'C_CO2 \t C_CH4 \t C_Ar \t C_O2 \t C_N2 \t N2/Ar \t O2/Ar\r\n');
fprintf(fileID,'%.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f\n',output);
fclose(fileID);
%type Atmos_conc.txt             % Display file contents

fileID = fopen('Atmos.txt','w');
fprintf(fileID,'%.2f \r\n%.2f \r\n%.2f \r\n%.2f \r\n%.2f \r\n%.2f \r\n%.2f\n',output');
fclose(fileID);