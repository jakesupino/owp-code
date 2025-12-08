%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcMeansSeasonality.m
% This script calculates the overall mean T, S, DO%sat, pH, depth, and
% tidal range across the entire time series for each site, and also the
% seasonality for the first four parameters.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 12/16/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%==========================================================================
% Import gap-filled data
%==========================================================================
site = {'north','gull','south'};

for i = 1:length(site)
    % Load parameter values
    cd([rootpath,'open-water-platform-data\all-sites'])
    load('gapFilled.mat');
    % Load two sigma table and define errors (σ = 1 standard deviation) for each term
    cd([rootpath,'open-water-platform-data\',site{i},'\cleaned\dupcheck'])
    load([site{i},'-bestGuess.mat'])   % Output from dupCheck scripts
    S_err = two_sigma.S/2; % [PSU]
    T_err = two_sigma.T/2; % [degC]
    DOconc_err = two_sigma.DOconc/2;  % [umol/L]
    pH_err = two_sigma.pH/2;
    % Load metab results - only need for the uncertainty estimates
    cd([rootpath,'diel-method\uncertainty-analysis\',site{i}])
    load('MonteCarloResults')
    % metab.(site{i}) = diel_dtd_MC;
    GPP_err = mean(diel_dtd_MC.GPP_sd);
    ER_err = mean(diel_dtd_MC.ER_sd);
    NEM_err = mean(diel_dtd_MC.NEM_sd);
    sigma.(site{i}) = table(T_err,S_err,DOconc_err,pH_err,GPP_err,ER_err,NEM_err,'VariableNames',{'T_err','S_err','DO_conc_err','pH_err','GPP_err','ER_err','NEM_err'});
end
clearvars finalQC bestguess two_sigma sigma_tbl diel_dtd diel_obs diel_dtd_MC counter S_err T_err DOconc_err pH_err GPP_err ER_err NEM_err

%====================================================================================
% Calculate parameter means and seasonal ranges across entire time series for Table 1
%====================================================================================
for i = 1:length(site)
    % Overall means
    T_mean = mean(daily_filled.(site{i}).temperature,'omitmissing');
    S_mean = mean(daily_filled.(site{i}).salinity,'omitmissing');
    DOsat_mean = mean(daily_filled.(site{i}).DOsat,'omitmissing');
    pH_mean = mean(daily_filled.(site{i}).pH,'omitmissing');
    depth_mean = mean(daily_filled.(site{i}).depth,'omitmissing');
    tidal_mean = mean(daily_filled.(site{i}).tidal,'omitmissing');
    GPP_mean = mean(daily_filled.(site{i}).GPP_avg,'omitmissing');
    ER_mean = mean(daily_filled.(site{i}).ER_avg,'omitmissing');
    NEM_mean = mean(daily_filled.(site{i}).NEM_avg,'omitmissing');

    % Find seasonality (only 2022 and 2023 are complete years)
    yr = year(daily_filled.(site{i}).Time);
    yr_21 = find(yr==2021);
    yr_22 = find(yr==2022);
    yr_23 = find(yr==2023);
    yr_24 = find(yr==2024);

    % Temperature
    hi_21 = prctile(daily_filled.(site{i}).temperature(yr_21(1):yr_21(end)),95);
    lo_22 = prctile(daily_filled.(site{i}).temperature(yr_22(1):yr_22(end)),5);
    hi_22 = prctile(daily_filled.(site{i}).temperature(yr_22(1):yr_22(end)),95);
    lo_23 = prctile(daily_filled.(site{i}).temperature(yr_23(1):yr_23(end)),5);
    hi_23 = prctile(daily_filled.(site{i}).temperature(yr_23(1):yr_23(end)),95);
    lo_24 = prctile(daily_filled.(site{i}).temperature(yr_24(1):yr_24(end)),5);
    T_lo = mean([lo_22, lo_23, lo_24]);
    T_hi = mean([hi_21, hi_22, hi_23]);

    % Salinity
    hi_21 = prctile(daily_filled.(site{i}).salinity(yr_21(1):yr_21(end)),95);
    lo_22 = prctile(daily_filled.(site{i}).salinity(yr_22(1):yr_22(end)),5);
    hi_22 = prctile(daily_filled.(site{i}).salinity(yr_22(1):yr_22(end)),95);
    lo_23 = prctile(daily_filled.(site{i}).salinity(yr_23(1):yr_23(end)),5);
    hi_23 = prctile(daily_filled.(site{i}).salinity(yr_23(1):yr_23(end)),95);
    lo_24 = prctile(daily_filled.(site{i}).salinity(yr_24(1):yr_24(end)),5);
    S_lo = mean([lo_22, lo_23, lo_24]);
    S_hi = mean([hi_21, hi_22, hi_23]);

    % DOsat
    hi_21 = prctile(daily_filled.(site{i}).DOsat(yr_21(1):yr_21(end)),95);
    lo_22 = prctile(daily_filled.(site{i}).DOsat(yr_22(1):yr_22(end)),5);
    hi_22 = prctile(daily_filled.(site{i}).DOsat(yr_22(1):yr_22(end)),95);
    lo_23 = prctile(daily_filled.(site{i}).DOsat(yr_23(1):yr_23(end)),5);
    hi_23 = prctile(daily_filled.(site{i}).DOsat(yr_23(1):yr_23(end)),95);
    lo_24 = prctile(daily_filled.(site{i}).DOsat(yr_24(1):yr_24(end)),5);
    DOsat_lo = mean([lo_22, lo_23, lo_24]);
    DOsat_hi = mean([hi_21, hi_22, hi_23]);

    % pH
    hi_21 = prctile(daily_filled.(site{i}).pH(yr_21(1):yr_21(end)),95);
    lo_22 = prctile(daily_filled.(site{i}).pH(yr_22(1):yr_22(end)),5);
    hi_22 = prctile(daily_filled.(site{i}).pH(yr_22(1):yr_22(end)),95);
    lo_23 = prctile(daily_filled.(site{i}).pH(yr_23(1):yr_23(end)),5);
    hi_23 = prctile(daily_filled.(site{i}).pH(yr_23(1):yr_23(end)),95);
    lo_24 = prctile(daily_filled.(site{i}).pH(yr_24(1):yr_24(end)),5);
    pH_lo = mean([lo_22, lo_23, lo_24]);
    pH_hi = mean([hi_21, hi_22, hi_23]);

    % GPP
    hi_21 = prctile(daily_filled.(site{i}).GPP_avg(yr_21(1):yr_21(end)),95);
    lo_22 = prctile(daily_filled.(site{i}).GPP_avg(yr_22(1):yr_22(end)),5);
    hi_22 = prctile(daily_filled.(site{i}).GPP_avg(yr_22(1):yr_22(end)),95);
    lo_23 = prctile(daily_filled.(site{i}).GPP_avg(yr_23(1):yr_23(end)),5);
    hi_23 = prctile(daily_filled.(site{i}).GPP_avg(yr_23(1):yr_23(end)),95);
    lo_24 = prctile(daily_filled.(site{i}).GPP_avg(yr_24(1):yr_24(end)),5);
    GPP_lo = mean([lo_22, lo_23, lo_24]);
    GPP_hi = mean([hi_21, hi_22, hi_23]);

    % ER
    hi_21 = prctile(daily_filled.(site{i}).ER_avg(yr_21(1):yr_21(end)),95);
    lo_22 = prctile(daily_filled.(site{i}).ER_avg(yr_22(1):yr_22(end)),5);
    hi_22 = prctile(daily_filled.(site{i}).ER_avg(yr_22(1):yr_22(end)),95);
    lo_23 = prctile(daily_filled.(site{i}).ER_avg(yr_23(1):yr_23(end)),5);
    hi_23 = prctile(daily_filled.(site{i}).ER_avg(yr_23(1):yr_23(end)),95);
    lo_24 = prctile(daily_filled.(site{i}).ER_avg(yr_24(1):yr_24(end)),5);
    ER_lo = mean([lo_22, lo_23, lo_24]);
    ER_hi = mean([hi_21, hi_22, hi_23]);

    % NEM
    hi_21 = prctile(daily_filled.(site{i}).NEM_avg(yr_21(1):yr_21(end)),95);
    lo_22 = prctile(daily_filled.(site{i}).NEM_avg(yr_22(1):yr_22(end)),5);
    hi_22 = prctile(daily_filled.(site{i}).NEM_avg(yr_22(1):yr_22(end)),95);
    lo_23 = prctile(daily_filled.(site{i}).NEM_avg(yr_23(1):yr_23(end)),5);
    hi_23 = prctile(daily_filled.(site{i}).NEM_avg(yr_23(1):yr_23(end)),95);
    lo_24 = prctile(daily_filled.(site{i}).NEM_avg(yr_24(1):yr_24(end)),5);
    NEM_lo = mean([lo_22, lo_23, lo_24]);
    NEM_hi = mean([hi_21, hi_22, hi_23]);

    means.(site{i}) = table(T_mean,S_mean,DOsat_mean,pH_mean,depth_mean,tidal_mean,GPP_mean,ER_mean,NEM_mean,'VariableNames',{'T','S','DOsat','pH','depth','tidal','GPP','ER','NEM'});
    seas_ranges.(site{i}) = table(T_lo,T_hi,S_lo,S_hi,DOsat_lo,DOsat_hi,pH_lo,pH_hi,GPP_lo,GPP_hi,ER_lo,ER_hi,NEM_lo,NEM_hi);
end
%%
% Find overall site-wide mean
param = "GPP_avg";
mean([daily_filled.north.(param),daily_filled.gull.(param),daily_filled.south.(param)],"all")

% Find mean site-wide seasonal range
lower = "GPP_lo";
upper = "GPP_hi";
mean([(seas_ranges.north.(upper) - seas_ranges.north.(lower)),(seas_ranges.gull.(upper) - seas_ranges.gull.(lower)),(seas_ranges.south.(upper) - seas_ranges.south.(lower))],"all")

%% OLD - Calculate metrics for Years 2 & 3 only
% June 28, 2022 starts on 15663 (North), 31227 (Gull), 33062 (South) for daily_filled
% June 28, 2022 starts on 63 (North), 157 (Gull), 176 (South) for metab

% %====Find parameter means and seasonality for Table 1======================
% for i = 1:length(site)
%     if strcmp(site(i),'north')
%         start = 15663;
%     elseif strcmp(site(i),'gull')
%         start = 31227;
%     elseif strcmp(site(i),'south')
%         start = 33062;
%     end
% 
%     % Overall means
%     T_mean = mean(daily_filled.(site{i}).temperature(start:end),'omitmissing');
%     S_mean = mean(daily_filled.(site{i}).salinity(start:end),'omitmissing');
%     DOsat_mean = mean(daily_filled.(site{i}).DOsat(start:end),'omitmissing');
%     pH_mean = mean(daily_filled.(site{i}).pH(start:end),'omitmissing');
%     depth_mean = mean(daily_filled.(site{i}).depth(start:end),'omitmissing');
%     tidal_mean = mean(tidal.(site{i}).daily_range,'omitmissing');
% 
%     % Find seasonality (only 2022 and 2023 are complete years)
%     yr = year(daily_filled.(site{i}).datetime_utc(start:end));
%     % yr_21 = find(yr==2021);
%     yr_22 = find(yr==2022);
%     yr_23 = find(yr==2023);
%     yr_24 = find(yr==2024);
% 
%     % Temperature
%     % hi_21 = prctile(daily_filled.(site{i}).temperature(yr_21(1):yr_21(end)),95);
%     lo_22 = prctile(daily_filled.(site{i}).temperature(yr_22(1):yr_22(end)),5);
%     hi_22 = prctile(daily_filled.(site{i}).temperature(yr_22(1):yr_22(end)),95);
%     lo_23 = prctile(daily_filled.(site{i}).temperature(yr_23(1):yr_23(end)),5);
%     hi_23 = prctile(daily_filled.(site{i}).temperature(yr_23(1):yr_23(end)),95);
%     lo_24 = prctile(daily_filled.(site{i}).temperature(yr_24(1):yr_24(end)),5);
%     T_lo = mean([lo_22, lo_23]);
%     % T_hi = mean([hi_21, hi_22, hi_23]);
%     T_hi = mean([hi_22, hi_23]);
% 
%     % Salinity
%     hi_21 = prctile(daily_filled.(site{i}).salinity(yr_21(1):yr_21(end)),95);
%     lo_22 = prctile(daily_filled.(site{i}).salinity(yr_22(1):yr_22(end)),5);
%     hi_22 = prctile(daily_filled.(site{i}).salinity(yr_22(1):yr_22(end)),95);
%     lo_23 = prctile(daily_filled.(site{i}).salinity(yr_23(1):yr_23(end)),5);
%     hi_23 = prctile(daily_filled.(site{i}).salinity(yr_23(1):yr_23(end)),95);
%     lo_24 = prctile(daily_filled.(site{i}).salinity(yr_24(1):yr_24(end)),5);
%     S_lo = mean([lo_22, lo_23, lo_24]);
%     S_hi = mean([hi_21, hi_22, hi_23]);
% 
%     % DOsat
%     hi_21 = prctile(daily_filled.(site{i}).DOsat(yr_21(1):yr_21(end)),95);
%     lo_22 = prctile(daily_filled.(site{i}).DOsat(yr_22(1):yr_22(end)),5);
%     hi_22 = prctile(daily_filled.(site{i}).DOsat(yr_22(1):yr_22(end)),95);
%     lo_23 = prctile(daily_filled.(site{i}).DOsat(yr_23(1):yr_23(end)),5);
%     hi_23 = prctile(daily_filled.(site{i}).DOsat(yr_23(1):yr_23(end)),95);
%     lo_24 = prctile(daily_filled.(site{i}).DOsat(yr_24(1):yr_24(end)),5);
%     DOsat_lo = mean([lo_22, lo_23, lo_24]);
%     DOsat_hi = mean([hi_21, hi_22, hi_23]);
% 
%     % pH
%     hi_21 = prctile(daily_filled.(site{i}).pH(yr_21(1):yr_21(end)),95);
%     lo_22 = prctile(daily_filled.(site{i}).pH(yr_22(1):yr_22(end)),5);
%     hi_22 = prctile(daily_filled.(site{i}).pH(yr_22(1):yr_22(end)),95);
%     lo_23 = prctile(daily_filled.(site{i}).pH(yr_23(1):yr_23(end)),5);
%     hi_23 = prctile(daily_filled.(site{i}).pH(yr_23(1):yr_23(end)),95);
%     lo_24 = prctile(daily_filled.(site{i}).pH(yr_24(1):yr_24(end)),5);
%     pH_lo = mean([lo_22, lo_23, lo_24]);
%     pH_hi = mean([hi_21, hi_22, hi_23]);
% 
%     param_means.(site{i}) = table(T_mean,S_mean,DOsat_mean,pH_mean,depth_mean,tidal_mean,'VariableNames',{'T','S','DOsat','pH','depth','tidal'});
%     param_ranges.(site{i}) = table(T_lo,T_hi,S_lo,S_hi,DOsat_lo,DOsat_hi,pH_lo,pH_hi);
% end