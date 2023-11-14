function [conc_CO2] = CO2sol(S,T)

%CO2sol   Solubility of CO2 in sea water

%=========================================================================
% Author: Emily Chua
%
% USAGE:  concCO2 = CO2sol(S,T)
%
% DESCRIPTION:
%    Solubility (saturation) of carbon dioxide (CO2) in sea water
%    at 1-atm pressure of air including saturated water vapor
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%
% OUTPUT:
%   concCO2 = solubility of CO2  [umol/kg] 
% 
% AUTHOR:  Emily Chua
%
% REFERENCES:
%    Weiss, 1974
%    "CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A NON-IDEAL GAS"

%=========================================================================

% CALLER: general purpose
% CALLEE: sw_dens_0.m

%----------------------
% Check input parameters
%----------------------
if nargin ~=2
   error('CO2sol.m: Must pass 2 parameters')
end %if

% Check S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) | (ns~=nt)) & (ms+ns>2) & (mt+nt>2)
   error('CO2sol: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

Tabs = T + 273.15;     % Absolute temperature

% FUGACITY CALCULATION

xCO2 = 0.000397;        % Mole fraction of CO2 in dry atmosphere; from http://ossfoundation.us/projects/environment/global-warming/atmospheric-composition

p = 101325;             % Pressure [Pa]

R = 8.31447;            % Gas constant 

B = -1636.75+12.0408*Tabs-3.27957E-2*Tabs^2+3.16528E-5*Tabs^3;   % Eq. 6: Virial cofficient of pure CO2 gas [cm^3 mol^-1]
delta = 57.7 - 0.118*Tabs;                                       % Eq. 11: Virial cofficient of CO2 in air [cm^3 mol^-1]

fg = xCO2*p/101325*exp((B*10^-6+2*(1-xCO2)^2*delta*10^-6)*p/(R*Tabs));   % Eq. 9: Fugacity [atm]

% Constants from Table 1 of Weiss (1974)
A1 = -60.2409;
A2 = 93.4517;
A3 = 23.3585;
B1 = 0.023517;
B2 = -0.023656;
B3 = 0.0047036;

% Eqn (12) of Weiss (1974)
conc_CO2 = fg*exp(A1 + A2*(100/Tabs) + A3*log(Tabs/100) + S*[B1 + B2*(Tabs/100) + B3*(Tabs/100)^2])*10^6;

display(['CO2 concentration = ',num2str(conc_CO2,3),' umol/kg']);

return
