%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getSfromAC.m
% This script converts the discrete sample actual conductivity (AC) measurements
% to salinity values using the AquaTroll equations*, and adds these values
% to the "Lab_Salinity" column in the OWP salinity inventory sheet.
%
% *See technical note: https://in-situ.com/pub/media/support/documents/Aqua-TROLL-200-Measurement-Methodology-Tech-Note.pdf?srsltid=AfmBOorUxzlwtKlwS-Q1qCIQ19E_crNfiskl3MdtCMqksoul7m_96G-q
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 9/26/2024
% Last updated: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

cd('G:\Shared drives\SMIIL\Shared Data')

sample_sal = readtable('OWP_bottle_no_JS_withAC.xlsx');

% Convert AC to S using AquaTroll equations
% Constants
a0 = 0.0080;
a1 = -0.1692;
a2 = 25.3851;
a3 = 14.0941;
a4 = -7.0261;
a5 = 2.7081;
b0 = 0.0005;
b1 = -0.0056;
b2 = -0.0066;
b3 = 0.0375;
b4 = 0.0636;
b5 = -0.0144;
r0 = 29752.63;
r1 = 830.5102;
r2 = 3.429338;
r3 = -0.02193934;

% Inputs
AC = sample_sal.MeasuredAC * 10^3; % [mS/cm] --> [uS/cm]
T = sample_sal.Temp;        % [degC]

R = AC ./ (r0 + r1*T + r2*T.^2 + r3*T.^3);
X = 400*R;
Y = 100 *R;
f = (T-15) ./ (1 + 0.0162*(T-15));

S_calc = a0 + a1*R.^(1/2) + a2*R + a3*R.^(3/2) + a4*R.^2 + a5*R.^(5/2)...
    + f.*(b0 + b1*R.^(1/2) + b2*R + b3*R.^(3/2) + b4*R.^2 + b5*R.^(5/2))...
    - a0 ./ (1 + 1.5*X + X.^2)...
    - b0 * f ./ (1 + Y.^(1/2) + Y.^(3/2)); % [PSU]

sample_sal.S_calc = S_calc;

% Create a new inventory sheet with only OWP samples and add the calculated lab S values
inventory = readtable('SMIIL DIC_TA Sample Salinities.xlsx');

ind_owp = ismember(inventory.Location_ID,{'Gull','North','South'});

owp_inventory = inventory(ind_owp,:);

ind_match = find(ismember(owp_inventory.Bottle_No,sample_sal.Bottle_No));

owp_inventory.Lab_Salinity(ind_match) = S_calc;

cd('G:\Shared drives\SMIIL\Shared Data')
writetable(owp_inventory,'SMIIL DIC_TA Sample Salinities - OWP.xlsx')
