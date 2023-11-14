%% INPUTS
t = 25.75;             % Temperature of solution [deg C]
S = 36.1;             % Salinity of solution [ppt]
x = 0.000404;        % Mole fraction of CO2 in dry atmosphere

T = t + 273.15;     % [K]
p = 1;              % [atm]
p = p*101325;       % [Pa]

R = 8.31447;
%% CALCULATIONS
% SOP 24 from "Guide to Best Practices for Ocean CO2 Measurements"
% (Dickson, Sabine, & Christian, 2007)

B = -1636.75+12.0408*T-3.27957E-2*T^2+3.16528E-5*T^3;   % Eq. 19: Virial coefficeint of pure CO2 gas [cm^3 mol^-1]
delta = 57.7 - 0.118*T;                                 % Eq. 20: Virial cofficient of CO2 in air [cm^3 mol^-1]

f = x*p/101325*exp((B*10^-6+2*(1-x)^2*delta*10^-6)*p/(R*T));   % Fugacity [atm]
