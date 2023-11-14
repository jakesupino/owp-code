function [conc_NO] = NOsol(S,T)

% NOsol   Solubility of NO in freshwater
%=========================================================================
% Author: Emily Chua
%
% USAGE:  concNO = NOsol(S,T)
%
% DESCRIPTION:
%    Solubility (saturation) of nitric oxide (NO) in freshwater
%    at 1-atm pressure of air including saturated water vapor
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%
% OUTPUT:
%   concCNO = solubility of NO  [umol/kg] 
%
% REFERENCE:
%    Tian et al., 2021
%    "Continuous Chemiluminescence Measurements of Dissolved Nitric Oxide (NO) and Nitrogen Dioxide (NO2) in the Ocean Surface Layer of the East China Sea"
%
%=========================================================================

% CALLER: general purpose
% CALLEE: sw_dens_0.m

%----------------------
% Check input parameters
%----------------------
if nargin ~=2
   error('NOsol.m: Must pass 2 parameters')
end %if

% Check S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) | (ns~=nt)) & (ms+ns>2) & (mt+nt>2)
   error('NOsol: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

% Constants from Eq. 5 of Tian et al. (2021)

% [nmol/L]
B0 = 23.9;
B1 = 212.3;
B2 = 3732.8;
B3 = 18976.8;
B4 = 52512.6;
B5 = 73102.2;
B6 = 40212.7;

Tabs = T + 273.15;     % Absolute temperature
theta = Tabs/100;
R = 8.31446;           % Gas constant [m3 Pa K-1 mol-1]

patm = 101325;         % [Pa]
pNO = patm;            % Pure NO gas

% Eqn (4) of Tian et al. (2021)
f_NO = pNO.*exp(((B0 + B1./theta - B2./theta.^2 + B3./theta.^3 - B4./theta.^4 + B5./theta.^5 + B6./theta.^6).*patm*10^-6)./(R*Tabs));

% Eqn (8) of Tian et al. (2021)
HCP = 1.9E-5*exp(1600*(1./Tabs - 1/298.15)); % mol m-3 Pa-1

conc_NO = f_NO.*HCP.*10^3; % [uM]   

return