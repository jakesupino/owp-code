%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataQC_gull.m
% This script performs data quality control on the adjusted merged sonde
% data for Gull.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% 11/9/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;clc

site = 'gull';

rootpath = 'G:\My Drive\Postdoc\';
cd([rootpath,'SMIIL\open-water-platform-data\',site,'\adjusted\merged'])

load(['alldeps-',site,'-adj.mat'])

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = max(sonde1_all.datetime_utc(end), sonde2_all.datetime_utc(end));
NumTicks = 13;
XTick = linspace(dt1,dt2,NumTicks);
XTickFormat = "M/yy";
XLabel = 'Month/Year';
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
FontSize = 14;
LineWidth = 1;

switch site
    case 'gull'
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
            'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
            'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
            'Deployment 15'};
    case 'north'
        label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8',...
            'Deployment 9','Deployment 10','Deployment 11','Deployment 12',...
            'Deployment 13','Deployment 14','Deployment 15'};
    case 'south'
        % BC
        label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
            'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
            'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
            'Deployment 14'};
        % ERDC
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 7',...
            'Deployment 8','Deployment 9','Deployment 10','Deployment 11',...
            'Deployment 12','Deployment 13','Deployment 14','Deployment 15'};
end

% Find indices of deployment changes
ind_dep = find(diff(sonde1_all.deployment) > 0);

%====Apply QC tests to flag data===========================================
% Make a table of flags for each sonde; start with everything as passing (1)
flags1 = ones(height(sonde1_all),width(sonde1_all));
flags1 = array2table(flags1);
flags1.Properties.VariableNames = {'deployment' 'datetime_utc' 'datetime_local' ...
    'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
    'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
    'pH' 'pH_raw' 'ORP' 'chla' 'nitrate' 'external_voltage' 'battery_capacity'};    % Or create flag columns for variables of interest only??

%% DEPTH
%====(1) Gross Range Test==================================================

% Lower limit - flag all negative values as FAIL
d_low = 0;      % [m]
ind_low = find(sonde1_all.depth < d_low);

% Upper limit - not sure what should be; set to 5 m for now
d_high = 5;     % [m]
ind_high = find(sonde1_all.depth > d_high);

ind_grossRange = [ind_low;ind_high];

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.depth,'.');
hold on
h1 = plot(sonde1_all.datetime_utc(ind_grossRange),sonde1_all.depth(ind_grossRange),'or');
x1 = xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
ylabel('Depth (m)')
title('QC Tests - BC Sonde - Depth')
legend([h0 h1],'Original Data','Gross Range Test','location','best')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)

%====(2) Spike test========================================================
d_depth = zeros(length(sonde1_all.depth),1);
fail_threshold = 0.1;    % FAIL threshold [m]

for i = 2:length(sonde1_all.depth)
    % Only check points that have not already been flagged in Gross Range Test
    if i ~= ind_grossRange
        ref = (sonde1_all.depth(i+1)+sonde1_all.depth(i-1))/2; % Calculate the average of i+1 and i-1
        d_depth(i) = sonde1_all.depth(i) - ref; % Calculate the difference between i and ref
        if abs(d_depth(i)) >= fail_threshold % If i - ref is greater than high threshold, FAIL
            flags1.depth(i) = 4;
        end
    end
end
ind_spike = find(abs(d_depth) >= fail_threshold);

fig2 = figure(2);clf
fig2.WindowState = 'maximized';
histogram(d_depth)
title('Spike Test Histogram - Depth')
xlabel('d_{depth} (m)')
set(gca,'YScale','log','FontSize',FontSize)

fig3 = figure(3);clf
fig3.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.depth,'.');
hold on
h1 = plot(sonde1_all.datetime_utc(ind_grossRange),sonde1_all.depth(ind_grossRange),'or');
h2 = plot(sonde1_all.datetime_utc(ind_spike),sonde1_all.depth(ind_spike),'om');
x1 = xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
ylabel('Depth (m)')
legend([h0 h1 h2],'Original Data','Gross Range Test','Spike Test','location','best')
title('QC Tests - BC Sonde - Depth')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)

%====Clean data============================================================
% Discard points that failed Gross Range and Spike Tests
sonde1_all.depth(ind_grossRange) = NaN;
sonde1_all.depth(ind_spike) = NaN;

%====(3) Manual removals===================================================
% See annotations in Collaborative Lab Notebook

% Discard isolated points during large gap from 7/30/21 (during Dep 1) to start of Dep 2
sonde1_all.depth(7420:12590) = NaN;
flags1.depth(7420:12590) = 4;

% Flag all NaNs as FAIL
ind_fail = find(isnan(sonde1_all.depth));
flags1.depth(ind_fail) = 4;

%====Convert table to timetable============================================
sonde1_TT = table2timetable(sonde1_all,'RowTimes','datetime_utc');
newTimeStep = minutes(10);
sonde1_TT = retime(sonde1_TT,'regular','mean','TimeStep',newTimeStep); % Calculate mean of values in each time bin

sonde1_TT.datetime_utc = sonde1_TT.datetime_utc + newTimeStep/2;
sonde1_TT.datetime_local = sonde1_TT.datetime_utc;
sonde1_TT.datetime_local.TimeZone = 'America/New_York';

fig4 = figure(4);clf
fig4.WindowState = 'maximized';
% h0 = plot(sonde1_all.datetime_utc,sonde1_all.depth,'.');
% hold on
h1 = plot(sonde1_TT.datetime_utc,sonde1_TT.depth,'.','Color',red);
x1 = xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
% legend([h0 h1],'Original','Retimed','location','best')
ylabel('Depth (m)')
title('Gull - BC Sonde - Cleaned Depth Data')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylim([0 3])
set(gca,'FontSize',FontSize)

%% TEMPERATURE
%====(1) Gross Range Test==================================================

% Lower limit - flag all negative values as FAIL
T_low = -5;    % degC
ind_low = find(sonde1_all.temperature < T_low);

% Upper limit - not sure what should be; set to 5 m for now
T_high = 35;    % degC
ind_high = find(sonde1_all.temperature > T_high);

ind_grossRange = [ind_low;ind_high];

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.temperature,'.');
hold on
h1 = plot(sonde1_all.datetime_utc(ind_grossRange),sonde1_all.temperature(ind_grossRange),'or');
x1 = xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
ylabel('Temperature (^oC)')
title('QC Tests - BC Sonde - Temperature')
legend([h0 h1],'Original Data','Gross Range Test','location','best')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)


