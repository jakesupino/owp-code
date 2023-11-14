%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data_qc-original.m
% This script performs data quality control on the adjusted merged sonde
% data.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% 11/3/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

fig = uifigure;
site = uiconfirm(fig,"Select the platform","Site selection","Options",["gull","north","south"]);
close(fig)

cd([rootpath,'open-water-platform-data\',site,'\adjusted\merged'])

load(['alldeps-',site,'-adj.mat'])

% Global plotting settings
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
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6','Deployment 7','Deployment 8','Deployment 9','Deployment 10','Deployment 11','Deployment 12','Deployment 13','Deployment 14','Deployment 15'};
    case 'north'
        label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8','Deployment 9','Deployment 10','Deployment 11','Deployment 12','Deployment 13','Deployment 14','Deployment 15'};
    case 'south'
        % BC
        label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5','Deployment 6','Deployment 7','Deployment 8','Deployment 9','Deployment 10','Deployment 11','Deployment 12','Deployment 13','Deployment 14'};
        % ERDC
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 7','Deployment 8','Deployment 9','Deployment 10','Deployment 11','Deployment 12','Deployment 13','Deployment 14','Deployment 15'};
end

% Find indices of deployment changes
ind_dep = find(diff(sonde1_all.deployment) > 0);

%====Check timestamps======================================================
timeDiff = diff(sonde1_all.datetime_utc);
ind_overlap = find(timeDiff < 0);
ind_gap = find(timeDiff > minutes(12)); % Sampling interval is 6, 10, or 12 min

%% 
% Convert table to timetable
sonde1_TT = table2timetable(sonde1_all,'RowTimes','datetime_utc');
newTimeStep = minutes(10);
% sonde1_TT = retime(sonde1_TT,'regular','fillwithmissing','TimeStep',newTimeStep);   % Fill gaps with NaN or NaT
sonde1_TT = retime(sonde1_TT,'regular','mean','TimeStep',newTimeStep);   % Calculate mean of values in each time bin

sonde1_TT.datetime_utc = sonde1_TT.datetime_utc + newTimeStep/2;

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.depth,'.');
hold on
h1 = plot(sonde1_TT.datetime_utc,sonde1_TT.depth,'.');
x1 = xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
legend([h0 h1],'Original','Retimed')
ylabel('Depth (m)')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)

% fig2 = figure(2);clf
% fig2.WindowState = 'maximized';
% h0 = plot(sonde1_all.datetime_utc,sonde1_all.depth,'.');
% hold on
% x1 = xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
% ylabel('Depth (m)')
% title('Original')
% xlim([dt1 dt2])                 % Use same x limits for comparing sites
% set(gca,'FontSize',FontSize)
%%
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.depth,'.');
hold on
h1 = plot(sonde1_all.datetime_utc(ind_overlap),sonde1_all.depth(ind_overlap),'xg','MarkerSize',10,'LineWidth',2);
h2 = plot(sonde1_all.datetime_utc(ind_gap),sonde1_all.depth(ind_gap),'xr','MarkerSize',10,'LineWidth',2);
x1 = xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
ylabel('Depth (m)')
title('QC Tests - Sonde 1 - Depth')
legend([h0 h1 h2],'Original Data','Overlap','Gap','location','best')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)

%%
% Make a table of flags for each sonde; start with everything as passing (1)
flags1 = ones(height(sonde1_all),width(sonde1_all));
flags1 = array2table(flags1);
flags1.Properties.VariableNames = {'deployment' 'datetime_utc' 'datetime_local' ...
    'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
    'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
    'pH' 'pH_raw' 'ORP' 'chla' 'nitrate' 'external_voltage' 'battery_capacity'};    % Or create flag columns for variables of interest only??


%====Gross range test======================================================
%----Depth-----------------------------------------------------------------
% Flag all NaNs as FAIL
ind1 = find(isnan(sonde1_all.depth));
flags1.depth(ind1) = 4; 

% Lower limit - flag all negative values as FAIL
ind2 = find(sonde1_all.depth < 0);
flags1.depth(ind2) = 4;    

% Upper limit - not sure what should be; set to 5 m for now
ind3 = find(sonde1_all.depth > 5);
flags1.depth(ind3) = 4;

ind_grossRange = [ind1;ind2;ind3];

fig2 = figure(1);clf
fig2.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.depth,'.');
hold on
h1 = plot(sonde1_all.datetime_utc(ind_grossRange),sonde1_all.depth(ind_grossRange),'or');
x1 = xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
ylabel('Depth (m)')
title('QC Tests - Sonde 1 - Depth')
legend([h0 h1],'Original Data','Gross Range Test','location','best')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)


%%
%====Spike test============================================================
%----Depth-----------------------------------------------------------------
d_depth = zeros(length(sonde1_all.depth),1);
fail_threshold = 0.1;    % FAIL threshold [m]
% suspect_threshold = 1; % Suspect threshold [m]

for i = 2:length(sonde1_all.depth)
    % Only check points that have not already been flagged in Gross Range Test
    if flags1.depth(i) ~= 4
        % Calculate the average of i+1 and i-1
        ref = (sonde1_all.depth(i+1)+sonde1_all.depth(i-1))/2;
        % Calculate the difference between i and ref
        d_depth(i) = sonde1_all.depth(i) - ref;
        % If i - ref is greater than high threshold, FAIL
        if abs(d_depth(i)) >= fail_threshold
            flags1.depth(i) = 4;
        end
    end
end
ind_spike = find(abs(d_depth) >= fail_threshold);  % For plotting data spikes

% Include SUSPECT threshold
% for i = 3:length(sonde1_all.depth)
%     d_depth(i-1) = sonde1_all.depth(i-1) - (sonde1_all.depth(i)+sonde1_all.depth(i-2))/2;
%     if abs(d_depth(i)) > suspect_threshold && abs(d_depth(i)) < fail_threshold
%         flags1.depth(i) = 3;
%     elseif abs(d_depth(i)) >= fail_threshold
%         flags1.depth(i) = 4;
%     end
%     if abs(d_depth(i-1)) >= fail_threshold
%         flags1.depth(i-1) = 4;
%     end
% end

fig3 = figure(2);clf
fig3.WindowState = 'maximized';
histogram(d_depth)
title('Spike Test Histogram - Depth')
xlabel('d_{depth} (m)')
set(gca,'YScale','log','FontSize',FontSize)

fig4 = figure(3);clf
fig4.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.depth,'.');
hold on
h1 = plot(sonde1_all.datetime_utc(ind_grossRange),sonde1_all.depth(ind_grossRange),'or');
h2 = plot(sonde1_all.datetime_utc(ind_spike),sonde1_all.depth(ind_spike),'om');
x1 = xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
ylabel('Depth (m)')
legend([h0 h1 h2],'Original Data','Gross Range Test','Spike Test','location','best')
title('QC Tests - Sonde 1 - Depth')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)

%%
%====Moving Mean===========================================================
window = 240;            % Sliding window length; 240 pts = 24 h for sampling frequency of 6 min

% Need to have a window with constant dt, not number of data points.
% See https://www.mathworks.com/help/matlab/matlab_prog/clean-timetable-with-missing-duplicate-or-irregular-times.html#responsive_offcanvas
% Convert table to time table
% sonde1_TT = table2timetable(sonde1_all);
% 
% dupTimes = sort(sonde1_TT.datetime_utc);
% 
% sonde1_TT = sortrows(sonde1_TT);
% sonde1_TT_reg = retime(sonde1_TT,'regular','linear','TimeStep',minutes(6));

ind_good = find(flags1.depth ~=4);
mmean1 = movmean(sonde1_all.depth(ind_good),window);
mstd1 = movstd(sonde1_all.depth(ind_good),window);
mmedian1 = movmedian(sonde1_all.depth(ind_good),window);
mmad1 = movmad(sonde1_all.depth(ind_good),window);

val1 = sonde1_all.depth(ind_good);
val2 = mmedian1-3*mmad1;

% ind_mflag = find(val1 > abs(val2));
ind_mflag = find(abs(val1) > abs(val2));

% ind_mflag = zeros(length(ind_good),1);
% for i = 1:length(ind_good)
    % ind_mflag(i) = find(sonde1_all.depth(ind_good(i)) < (mmedian1(ind_good(i))-3*mmad1(ind_good(i))));
    % val = find(sonde1_all.depth(ind_good(i)) > abs(mmedian1(i)-3*mmad1(i)));
%     ind_mflag = find(val1(i) > abs(val2(i)));
% end

figure,clf
plot(sonde1_all.datetime_utc(ind_good),val1,'.')
hold on
plot(sonde1_all.datetime_utc(ind_good),val2,'.k')
plot(sonde1_all.datetime_utc(ind_good(ind_mflag)),sonde1_all.depth(ind_good(ind_mflag)),'.r')
hold off

fig5 = figure(3);clf
fig5.WindowState = 'maximized';
h0 = plot(sonde1_all.datetime_utc,sonde1_all.depth,'.');
hold on
h1 = plot(sonde1_all.datetime_utc,mean1);
x1 = xline([sonde1_all.datetime_utc(1); sonde1_all.datetime_utc(ind_dep+1)],'--',label);
ylabel('Depth (m)')
legend([h0 h1],'Original Data','Moving Mean','location','best')
title('QC Tests - Sonde 1 - Depth')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',FontSize)

%%
% Use groupsummary to bin the data by day and calculate the daily mean and standard deviation
dailyMeans_sonde1 = groupsummary(sonde1_all,'datetime_utc',"day","mean");
dailyStd_sonde1 = groupsummary(sonde1_all,'datetime_utc',"day","std");

figure,clf;plot(dailyMeans_sonde1.day_datetime_utc,dailyMeans_sonde1.mean_depth,'.');ylabel('Daily mean depth (m)')
histogram(dailyMeans_sonde1.mean_depth);xlabel('Daily mean depth (m)')