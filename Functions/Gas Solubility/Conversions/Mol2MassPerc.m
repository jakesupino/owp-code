% Script to convert mole percent to mass percent

clear all
close all
clc

%% INPUTS
% Mole percent in dry gas
mpCO2 = 0; 
mpCH4 = 0.25; 
mpAr = 0; 
mpO2 = 0; 
mpN2 = 100 - (mpO2 + mpAr + mpCO2 + mpCH4);   % N2 is the balance gas

%% CONSTANTS
% Molar mass [g/mol]
mmCO2 = 44.01;
mmCH4 = 16.04;
mmAr = 39.948;
mmO2 = 31.998;
mmN2 = 28.0134;

%% CONVERSION

% Individual masses [g]
mCO2 = mpCO2*mmCO2;
mCH4 = mpCH4*mmCH4;
mAr = mpAr*mmAr;
mO2 = mpO2*mmO2;
mN2 = mpN2*mmN2;

% Total mass [g]
%mtot = mCO2 + mCH4 + mAr + mO2 + mN2;
mtot = mCO2 + (100-mpCO2)*mmN2;

% Mass percent [%]
CO2 = mCO2/mtot*100;
CH4 = mCH4/mtot*100;
Ar = mAr/mtot*100;
O2 = mO2/mtot*100;
N2 = mN2/mtot*100;

% Display results
names = {'CO2 %','CH4 %','Ar %','O2 %','N2 %'};
val = {CO2, CH4, Ar, O2, N2};
format bank;                % Display 3 decimal places (MFC resolution is 0.01%)
output = [names; val]
