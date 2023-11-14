function [conc_CH4] = CH4sol(S,T)

% CH4sol   Solubility of CH4 in sea water
%=========================================================================
% Author: Emily Chua
%
% USAGE:  concCH4 = CH4sol(S,T)
%
% DESCRIPTION:
%    Solubility (saturation) of methane (CH4) in sea water
%    at 1-atm pressure of air including saturated water vapor
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%
% OUTPUT:
%   concCH4 = solubility of CH4  [umol/kg] 
%
% REFERENCE:
%    Wiesenburg and Garcia, 1979
%    "Equilibrium Solubilities of Methane, Carbon Monoxide, and Hydrogen in Water and Sea Water"
%
%=========================================================================

% CALLER: general purpose
% CALLEE: sw_dens_0.m

%----------------------
% Check input parameters
%----------------------
if nargin ~=2
   error('CH4sol.m: Must pass 2 parameters')
end %if

% Check S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) | (ns~=nt)) & (ms+ns>2) & (mt+nt>2)
   error('CH4sol: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

% Constants from Table 6 of Wiesenburg and Guinasso (1979)

% [nmol/L]
A1 = -415.2807;
A2 = 596.8104;
A3 = 379.2599;
A4 = -62.0757;
B1 = -0.059160;
B2 = 0.032174;
B3 = -0.0048198;

% [nmol/kg]
% A1 = -417.5053;
% A2 = 599.8626;
% A3 = 380.3636;
% A4 = -62.0764;
% B1 = -0.064236;
% B2 = 0.034980;
% B3 = -0.0052732;

Tabs = T + 273.15;     % Absolute temperature

fg = 0.00000179;    % Atmospheric CH4 concentration; from http://ossfoundation.us/projects/environment/global-warming/atmospheric-composition

% Eqn (7) of Wiesenburg and Guinasso (1979)
conc_CH4 = fg*exp(A1 + A2*(100./Tabs) + A3*log(Tabs./100) + A4*(Tabs./100) + S.*[B1 + B2*(Tabs./100) + B3*(Tabs./100).^2])./10^3;

% display(['CH4 concentration = ',num2str(conc_CH4*1000,3),' nmol/kg']);

return