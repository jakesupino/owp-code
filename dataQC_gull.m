%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataQC_gull.m
% This script performs data quality control on the adjusted merged sonde
% data for Gull.
%
% Code that requires manual input is commented with "INPUTS".
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 11/9/2023
% Last updated: 1/29/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

site = 'gull';

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

cd([rootpath,'open-water-platform-data\',site,'\adjusted\merged'])

load(['alldeps-',site,'-adj.mat'])

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = max(sonde1_all.datetime_utc(end), sonde2_all.datetime_utc(end));
XLabel = 'UTC';
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
FontSize = 14;
LineWidth = 1;

switch site
    case 'gull'
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
            'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
            'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
            'Deployment 15','Deployment 16'};
    % case 'north'
    %     label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8',...
    %         'Deployment 9','Deployment 10','Deployment 11','Deployment 12',...
    %         'Deployment 13','Deployment 14','Deployment 15'};
    % case 'south'
    %     % BC
    %     label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
    %         'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
    %         'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
    %         'Deployment 14'};
    %     % ERDC
    %     label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 7',...
    %         'Deployment 8','Deployment 9','Deployment 10','Deployment 11',...
    %         'Deployment 12','Deployment 13','Deployment 14','Deployment 15'};
end

% Find indices of deployment changes
ind_dep = find(diff(sonde1_all.deployment) > 0);

%====Create a table of flags for each sonde================================
% Start with everything as passing (flag = 1)
flags1 = ones(height(sonde1_all),width(sonde1_all));
flags1 = array2table(flags1);
flags1.Properties.VariableNames = {'deployment' 'datetime_utc' 'datetime_local' ...
    'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
    'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
    'pH' 'pH_raw' 'ORP' 'chla' 'nitrate' 'external_voltage' 'battery_capacity'};    % Or create flag columns for variables of interest only??

%====(0) Manual removals (see notes in Collab Lab Notebook)================
% Define "out-of-water" times
% These times were identified by running the QC tests (Gross Range and Spike) on depth, and then
% visually picking out the "obviously" out-of-water points.
oow1 = (34474:34476)'; % Dep 5 --> 6
oow2 = (51466:51486)'; % Dep 6 --> 7
oow3 = (60821:60829)'; % Dep 7 --> 8
oow4 = (67605:67643)'; % Dep 8 --> 9
oow5 = (71503:71528)'; % Dep 9 --> 10
oow6 = (78899); % Dep 10 --> 11 
oow7 = (87962:87965)'; % Dep 11 --> 12
oow8 = (104299:104326)'; % Dep 12 --> 13
oow9 = (113499:113501)'; % Dep 13 --> 14
oow10 = (122403:122405)'; % Dep 14 --> 15
oow11 = (130652:130661)'; % Dep 15 --> 16
ind_oow = [oow1;oow2;oow3;oow4;oow5;oow6;oow7;oow8;oow9;oow10;oow11];

% INPUTS
% Isolated points during large gap from 7/30/21 (during Dep 1) to start of Dep 2
ind_manual1 = (7420:12590)';
ind_manual2 = find(sonde1_all.depth(1:ind_dep(2)) < 0);
ind_manual_global = unique([ind_manual1;ind_manual2]);

clearvars oow1 oow2 oow3 oow4 oow5 oow6 oow7 oow8 oow9 oow10 oow11

cd([rootpath,'\figures\open-water-platform-figures\gull\data-qc\bc'])

% Plotting options
dotsize = 2;
circlesize = 6;
%% DEPTH
depth1_orig = sonde1_all.depth;  % Preserve original depth data for plotting

%====(1) Gross Range Test==================================================
% INPUTS
d_low = 0;      % Lower limit (m) - flag all negative values as FAIL
d_high = 5;     % Upper limit (m) - not sure what should be

ind_low = find(sonde1_all.depth < d_low);
ind_high = find(sonde1_all.depth > d_high);
% Don't duplicate points already defined to be manually removed
ind_low = setdiff(ind_low,ind_manual_global);   
ind_high = setdiff(ind_high,ind_manual_global);
ind_grossRange = [ind_low;ind_high];
ind_grossRange = setdiff(ind_grossRange,ind_oow);

ind_manual = [ind_manual_global;ind_oow];

%====(2) Spike Test========================================================
% INPUTS
spike_threshold = 0.1;    % FAIL threshold (m)

d_depth = zeros(length(sonde1_all.depth),1);

for i = 2:length(sonde1_all.depth)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (sonde1_all.depth(i+1)+sonde1_all.depth(i-1))/2; % Calculate the average of i+1 and i-1
        d_depth(i) = sonde1_all.depth(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_depth) >= spike_threshold);
[C,ig,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of depth differences
fig = figure;clf
fig.WindowState = 'maximized';
histogram(d_depth)
title('Spike Test Histogram - Depth')
xlabel('d_{depth} (m)')
set(gca,'YScale','log','FontSize',FontSize)

% Plot original data with flagged points
fig = figure;clf
fig.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.depth,'.k','markersize',dotsize);
hold on
h1 = plot(sonde1_all.datetime_utc(ind_manual_global),sonde1_all.depth(ind_manual_global),'og','markersize',circlesize);
h2 = plot(sonde1_all.datetime_utc(ind_oow),sonde1_all.depth(ind_oow),'oc','markersize',circlesize);
h3 = plot(sonde1_all.datetime_utc(ind_grossRange),sonde1_all.depth(ind_grossRange),'or','markersize',circlesize);
h4 = plot(sonde1_all.datetime_utc(ind_spike),sonde1_all.depth(ind_spike),'om','markersize',circlesize);
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
hold off
ylabel('Depth (m)')
legend([h0 h1 h2 h3 h4],'Original Data','Manual Removal','Out of Water','Gross Range Test','Spike Test','location','best')
title('QC Tests: Gull BC Sonde - Depth')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)

%====Clean data============================================================
% Discard points that failed Gross Range and Spike Tests and do Manual Removals
sonde1_all.depth(ind_manual_global) = NaN;
sonde1_all.depth(ind_grossRange) = NaN;
sonde1_all.depth(ind_spike) = NaN;
sonde1_all.depth(ind_oow) = NaN;

% Flag all NaNs as FAIL
ind_fail = find(isnan(sonde1_all.depth));
flags1.depth(ind_fail) = 4;

% Create structure to save flagged depth indices
depths_flagged = struct('ind_manual',ind_manual_global,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'ind_oow',ind_oow,'d_high',d_high,'d_low',d_low,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike C ig is

% Check that there aren't redundant points flagged in the Depth QC tests
% depths_flagged_all = unique([depths_flagged.ind_grossRange;depths_flagged.ind_manual;depths_flagged.ind_oow;depths_flagged.ind_spike]);
% depths_flagged_all0 = [depths_flagged.ind_grossRange;depths_flagged.ind_manual;depths_flagged.ind_oow;depths_flagged.ind_spike];

%===Clean pressure data based on cleaned depth data========================
% p1_adj = sonde1.density.*1000*9.81.*depth1_adj/6894.76;
sonde1_all.p = mean(sonde1_all.density,'omitmissing').*1000*9.81.*sonde1_all.depth/6894.76;
flags1.p(ind_fail) = 4;

%% TEMPERATURE
T1_orig = sonde1_all.temperature;  % Preserve original T data for plotting

% Plot the temperature data with the depth flags
% fig=figure;clf
% fig.WindowState = 'maximized';
% h0 = plot(sonde1_all.datetime_utc,T1_orig,'.k');
% hold on
% h1 = plot(sonde1_all.datetime_utc(depths_flagged_all),T1_orig(depths_flagged_all),'or');
% xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
% hold off
% legend([h0 h1],'Original Data','Depth Flags','location','best')
% ylabel('Temperature (^oC)')
% title('QC Tests: Gull BC Sonde - Temperature')
% xlim([dt1 dt2])                 % Use same x limits for comparing sites 
% set(gca,'FontSize',FontSize)

%====(1) Gross Range Test==================================================
% INPUTS -- could make fancier by changing limits based on season??
T_low = -2;      % Lower limit (oC) - Freezing pt of SW at 35 ppt (http://www.csgnetwork.com/h2ofreezecalc.html)
T_high = 35;     % Upper limit (oC)

ind_low = find(sonde1_all.temperature < T_low);
ind_high = find(sonde1_all.temperature > T_high);
ind_low = setdiff(ind_low,ind_manual_global);
ind_high = setdiff(ind_high,ind_manual_global);
ind_grossRange = [ind_low;ind_high];

%====(2) Spike Test========================================================
% INPUTS
spike_threshold = 1;    % FAIL threshold (deg C)

d_T = zeros(length(sonde1_all.temperature),1);

for i = 2:length(sonde1_all.temperature)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (sonde1_all.temperature(i+1)+sonde1_all.temperature(i-1))/2; % Calculate the average of i+1 and i-1
        d_T(i) = sonde1_all.temperature(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_T) >= spike_threshold);
[C,ig,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of temperature differences
fig = figure;clf
fig.WindowState = 'maximized';
histogram(d_T)
title('Spike Test Histogram - Temperature')
xlabel('d_{temperature} (^oC)')
set(gca,'YScale','log','FontSize',FontSize)

% Plot original data with flagged points
fig = figure;clf
fig.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.temperature,'.k','MarkerSize',dotsize);
hold on
h1 = plot(sonde1_all.datetime_utc(ind_manual_global),sonde1_all.temperature(ind_manual_global),'og','MarkerSize',circlesize);
h2 = plot(sonde1_all.datetime_utc(ind_grossRange),sonde1_all.temperature(ind_grossRange),'or','MarkerSize',circlesize);
h3 = plot(sonde1_all.datetime_utc(ind_spike),sonde1_all.temperature(ind_spike),'om','MarkerSize',circlesize);
h4 = plot(sonde1_all.datetime_utc(ind_oow),sonde1_all.temperature(ind_oow),'oc','MarkerSize',circlesize);
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
hold off
ylabel('Temperature (^oC)')
legend([h0 h1 h2 h3 h4],'Original Data','Manual Removal','Gross Range Test','Spike Test','Out of Water','location','best')
title('QC Tests: Gull BC Sonde - Temperature')
xlim([dt1 dt2])                 % Use same x limits for comparing sites 
set(gca,'FontSize',FontSize)

%====Clean data============================================================
% Discard points that failed Gross Range and Spike Tests and do Manual Removals
sonde1_all.temperature(ind_grossRange) = NaN;
sonde1_all.temperature(ind_manual_global) = NaN;
sonde1_all.temperature(ind_spike) = NaN;
sonde1_all.temperature(ind_oow) = NaN;

% Flag all NaNs as FAIL
ind_fail = find(isnan(sonde1_all.temperature));
flags1.temperature(ind_fail) = 4;

% Create structure to save flagged temperature indices
T_flagged = struct('ind_manual',ind_manual_global,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'ind_oow',ind_oow,'d_high',d_high,'d_low',d_low,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike C ig is

%% DO
DOconc1_orig = sonde1_all.DO_conc;  % Preserve original DO concentration data for plotting

%====(1) Gross Range Test==================================================
% INPUTS
DO_low = 0;      % Lower limit (umol/L)
DO_high = 500;     % Upper limit (umol/L); based on Winkler data and O2 solubility given S & T

ind_low = find(sonde1_all.DO_conc < DO_low);
ind_high = find(sonde1_all.DO_conc > DO_high);
ind_low = setdiff(ind_low,ind_manual_global);
ind_high = setdiff(ind_high,ind_manual_global);
ind_grossRange = [ind_low;ind_high];

%====(2) Spike Test========================================================
% INPUTS
spike_threshold = 50;    % FAIL threshold (umol/L)

d_DO = zeros(length(sonde1_all.DO_conc),1);

for i = 2:length(sonde1_all.DO_conc)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (sonde1_all.DO_conc(i+1)+sonde1_all.DO_conc(i-1))/2; % Calculate the average of i+1 and i-1
        d_DO(i) = sonde1_all.DO_conc(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_DO) >= spike_threshold);
[C,ig,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of DO concentration differences
fig = figure;clf
fig.WindowState = 'maximized';
histogram(d_DO)
title('Spike Test Histogram - DO Concentration')
xlabel('d_{DO conc} (\mumol/L)')
set(gca,'YScale','log','FontSize',FontSize)

% Plot original data with flagged points
fig = figure;clf
fig.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.DO_conc,'.k','MarkerSize',dotsize);
hold on
h1 = plot(sonde1_all.datetime_utc(ind_manual_global),sonde1_all.DO_conc(ind_manual_global),'og','MarkerSize',circlesize);
h2 = plot(sonde1_all.datetime_utc(ind_grossRange),sonde1_all.DO_conc(ind_grossRange),'or','MarkerSize',circlesize);
h3 = plot(sonde1_all.datetime_utc(ind_spike),sonde1_all.DO_conc(ind_spike),'om','MarkerSize',circlesize);
h4 = plot(sonde1_all.datetime_utc(ind_oow),sonde1_all.DO_conc(ind_oow),'oc','MarkerSize',circlesize);
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
hold off
ylabel('DO Concentration (\mumol/L)')
legend([h0 h1 h2 h3 h4],'Original Data','Manual Removal','Gross Range Test','Spike Test','Out of Water','location','best')
title('QC Tests: Gull BC Sonde - DO Concentration')
xlim([dt1 dt2])                 % Use same x limits for comparing sites--09876543 
set(gca,'FontSize',FontSize)

%====Clean data============================================================
% Discard points that failed Gross Range and Spike Tests and do Manual Removals
sonde1_all.DO_conc(ind_manual_global) = NaN;
sonde1_all.DO_conc(ind_grossRange) = NaN;
sonde1_all.DO_conc(ind_spike) = NaN;
sonde1_all.DO_conc(ind_oow) = NaN;

% Flag all NaNs as FAIL
ind_fail = find(isnan(sonde1_all.DO_conc));
flags1.DO_conc(ind_fail) = 4;

% Create structure to save flagged DO concentration indices
DOconc_flagged = struct('ind_manual',ind_manual_global,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'ind_oow',ind_oow,'d_high',d_high,'d_low',d_low,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike C ig is

%%
% SALINITY
cd([rootpath,'\figures\open-water-platform-figures\gull\data-qc\bc'])

%====(1) Gross Range Test==================================================
% INPUTS
S_low = 0;      % Lower limit (psu) - flag all negative values as FAIL
S_high = 50;     % Upper limit (psu) - not sure what should be

ind_low = find(sonde1_all.salinity < S_low);
ind_high = find(sonde1_all.salinity > S_high);
% Don't duplicate points already defined to be manually removed
ind_low = setdiff(ind_low,ind_manual_global);   
ind_high = setdiff(ind_high,ind_manual_global);
ind_grossRange = [ind_low;ind_high];
ind_grossRange = setdiff(ind_grossRange,ind_oow);

ind_manual = [ind_manual_global;ind_oow];

%====(2) Spike Test========================================================
% INPUTS
spike_threshold = 10;    % FAIL threshold (psu)

d_salinity = zeros(length(sonde1_all.salinity),1);

for i = 2:length(sonde1_all.salinity)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (sonde1_all.salinity(i+1)+sonde1_all.salinity(i-1))/2; % Calculate the average of i+1 and i-1
        d_salinity(i) = sonde1_all.salinity(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_salinity) >= spike_threshold);
[C,ig,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of salinity differences
fig = figure;clf
fig.WindowState = 'maximized';
histogram(d_salinity)
title('Spike Test Histogram - Salinity')
xlabel('d_{salinity} (psu)')
set(gca,'YScale','log','FontSize',FontSize)

% Plot original data with flagged points
fig = figure;clf
fig.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.salinity,'.k','MarkerSize',dotsize);
hold on
h1 = plot(sonde1_all.datetime_utc(ind_manual_global),sonde1_all.salinity(ind_manual_global),'og','MarkerSize',circlesize);
h2 = plot(sonde1_all.datetime_utc(ind_grossRange),sonde1_all.salinity(ind_grossRange),'or','MarkerSize',circlesize);
h3 = plot(sonde1_all.datetime_utc(ind_spike),sonde1_all.salinity(ind_spike),'om','MarkerSize',circlesize);
h4 = plot(sonde1_all.datetime_utc(ind_oow),sonde1_all.salinity(ind_oow),'oc','MarkerSize',circlesize);
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label)
hold off
ylabel('Salinity (psu)')
legend([h0 h1 h2 h3 h4],'Original Data','Manual Removal','Gross Range Test','Spike Test','Out of Water','location','best')
title('QC Tests: Gull BC Sonde - Salinity')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)

%====Clean data============================================================
% Discard points that failed Gross Range and Spike Tests and do Manual Removals
sonde1_all.salinity(ind_manual_global) = NaN;
sonde1_all.salinity(ind_grossRange) = NaN;
sonde1_all.salinity(ind_spike) = NaN;
sonde1_all.salinity(ind_oow) = NaN;

% Flag all NaNs as FAILind_fail = find(isnan(sonde1_all.salinity));
flags1.salinity(ind_fail) = 4;

% Create structure to save flagged salinity indices
S_flagged = struct('ind_manual',ind_manual_global,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'ind_oow',ind_oow,'S_high',S_high,'S_low',S_low,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike C ig is

%%
%====Convert table to timetable============================================
sonde1_TT = table2timetable(sonde1_all,'RowTimes','datetime_utc');
newTimeStep = minutes(10);
sonde1_TT = retime(sonde1_TT,'regular','mean','TimeStep',newTimeStep); % Calculate mean of values in each time bin

sonde1_TT.datetime_utc = sonde1_TT.datetime_utc + newTimeStep/2;
sonde1_TT.datetime_local = sonde1_TT.datetime_utc;
sonde1_TT.datetime_local.TimeZone = 'America/New_York';

%====Plot the cleaned data=================================================
% Depth
fig = figure;clf
fig.WindowState = 'maximized';
% h0 = plot(sonde1_all.datetime_utc,depth1_orig,'.','Color',red);
% hold on
h1 = plot(sonde1_TT.datetime_utc,sonde1_TT.depth,'.k','MarkerSize',dotsize);
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
% hold off
% legend([h0 h1],'Original','Retimed','location','best')
ylabel('Depth (m)')
title('Gull BC Sonde - Cleaned & Retimed Depth Data')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylim([-.5 3])
set(gca,'FontSize',FontSize)

% Temperature
fig = figure;clf
fig.WindowState = 'maximized';
% h0 = plot(sonde1_all.datetime_utc,T1_orig,'.','Color',red);
% hold on
h1 = plot(sonde1_TT.datetime_utc,sonde1_TT.temperature,'.k','MarkerSize',dotsize);
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
% hold off
% legend([h0 h1],'Original','Retimed','location','best')
ylabel('Temperature (^oC)')
title('Gull BC Sonde - Cleaned & Retimed Temperature Data')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)

% DO
fig = figure;clf
fig.WindowState = 'maximized';
% h0 = plot(sonde1_all.datetime_utc,DOconc1_orig,'.','Color',red);
% hold on
h1 = plot(sonde1_TT.datetime_utc,sonde1_TT.DO_conc,'.k','MarkerSize',dotsize);
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
% hold off
% legend([h0 h1],'Original','Retimed','location','best')
ylabel('DO Concentration (\mumol/L)')
title('Gull BC Sonde - Cleaned & Retimed DO Concentration Data')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)

% Salinity
fig = figure;clf
fig.WindowState = 'maximized';
% h0 = plot(sonde1_all.datetime_utc,depth1_orig,'.','Color',red);
% hold on
h1 = plot(sonde1_TT.datetime_utc,sonde1_TT.salinity,'.k','MarkerSize',dotsize);
xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
% hold off
% legend([h0 h1],'Original','Retimed','location','best')
ylabel('Salinity (psu)')
title('Gull BC Sonde - Cleaned & Retimed Salinity Data')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylim([-.5 50])
set(gca,'FontSize',FontSize)

%====Save the cleaned data=================================================
sonde1_cleaned = sonde1_TT;
cd([rootpath,'open-water-platform-data\gull\cleaned'])
save('gull-cleaned.mat','sonde1_cleaned','flags1')
