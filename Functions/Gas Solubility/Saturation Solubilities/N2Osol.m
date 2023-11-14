function [conc_N2O] = N2Osol(S,T)

% N2Osol   Solubility of N2O in sea water
%=========================================================================
% AUTHOR: Emily Chua
%
% USAGE:  concN2O = N2Osol(S,T)
%
% DESCRIPTION:
%    Solubility (saturation) of nitrous oxide (N2O) in sea water
%    at 1-atm pressure of air including saturated water vapor
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%
% OUTPUT:
%   concN2O = solubility of N2O  [umol/kg] 
% 
% AUTHOR:  Emily Chua
%
% REFERENCE:
%    Weiss and Price, 1980
%    "NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER"

%=========================================================================

% CALLER: general purpose
% CALLEE: sw_dens_0.m

%----------------------
% Check input parameters
%----------------------
if nargin ~=2
   error('N2Osol.m: Must pass 2 parameters')
end %if

% Check S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) | (ns~=nt)) & (ms+ns>2) & (mt+nt>2)
   error('N2Osol: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

Tabs = T + 273.15;     % Absolute temperature

% Constants from Table 2 of Weiss and Price (1980)
% [mol/L]
A1 = -165.8806;
A2 = 222.8743;
A3 = 92.0792;
A4 = -1.48425;
B1 = -0.056235;
B2 = 0.031619;
B3 = -0.0048472;

% [mol/kg]
% A1 = -168.2459;
% A2 = 226.0894;
% A3 = 93.2817;
% A4 = -1.48693;
% B1 = -0.060361;
% B2 = 0.033765;
% B3 = -0.0051862;

fg = 0.0000325/100;  % Atmospheric N2O concentration; from http://ossfoundation.us/projects/environment/global-warming/atmospheric-composition

% Eqn (13) of Weiss and Price (1980)
conc_N2O = fg*exp(A1 + A2*(100/Tabs) + A3*log(Tabs/100) + A4*(Tabs/100)^2 + S*[B1 + B2*(Tabs/100) + B3*(Tabs/100)^2])*10^6;

display(['N2O concentration = ',num2str(conc_N2O*1000,3),' nmol/kg']);

return