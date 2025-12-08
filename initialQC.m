%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialQC.m
% This script performs intial quality control tests on the adjusted merged sonde
% data for the selected site and sonde.
%
% Code that requires manual input is commented with "INPUTS".
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 11/9/2023
% Last updated: 11/9/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

prompt = {'Choose which site to QC'};
site = questdlg(prompt,'Sonde Selection','Gull','North','South','Gull');

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

cd([rootpath,'open-water-platform-data\',site,'\adjusted\merged'])

load(['alldeps-',site,'-adj.mat'])

prompt = {'Choose which sonde to QC'};
sonde = questdlg(prompt,'Sonde Selection','BC','ERDC','Cancel','BC');

switch sonde
    case 'BC'
        dat = sonde1_all;
    case 'ERDC'
        dat = sonde2_all;
end

clear sonde1_all sonde2_all

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = dat.datetime_utc(end);
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 12;
circlesize = 6;

switch site
    case 'Gull'
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
            'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
            'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
            'Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
    case 'North'
        label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8',...
            'Deployment 9','Deployment 10','Deployment 11','Deployment 12',...
            'Deployment 13','Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
    case 'South'
        switch sonde
            case 'BC'
                label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
                    'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
                    'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
                    'Deployment 14','Deployment 16','Deployment 17','Deployment 18'};
            case 'ERDC'
                label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 7',...
                    'Deployment 8','Deployment 9','Deployment 10','Deployment 11',...
                    'Deployment 12','Deployment 13','Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
        end
end

% Find indices of deployment changes
ind_dep = find(diff(dat.deployment) > 0);

%====Create a table of flags for each sonde================================
% Start with everything as passing (flag = 1)
flags = ones(height(dat),width(dat));
flags = array2table(flags);
switch sonde
    case 'BC'
        flags.Properties.VariableNames = {'deployment' 'datetime_utc' 'datetime_local' ...
            'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
            'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
            'pH' 'pH_raw' 'ORP' 'chla' 'nitrate' 'external_voltage' 'battery_capacity'};
    case 'ERDC'
        flags.Properties.VariableNames = {'deployment' 'datetime_utc' 'datetime_local' ...
            'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
            'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
            'pH' 'pH_raw' 'ORP' 'turbidity' 'external_voltage' 'battery_capacity'};
end

%====(0) Manual removals (see notes in Collab Lab Notebook)================
% Plot original depth data to define points to remove manually
fig0 = figure(1);clf
fig0.WindowState = 'maximized';
plot(dat.datetime_utc,dat.depth,'.k','markersize',12);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);
ylabel('Depth (m)')
title(['Original Data - ',site,' ',sonde])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);

% INPUTS
switch site
    case 'Gull'
        switch sonde
            case 'BC'
                % Out-of-water points
                oow1 = (34474:34476)';      % Dep 5 --> 6
                oow2 = (51466:51486)';      % Dep 6 --> 7
                oow3 = (60821:60829)';      % Dep 7 --> 8
                oow4 = (67605:67643)';      % Dep 8 --> 9
                oow5 = (71503:71528)';      % Dep 9 --> 10
                oow6 = (78899);             % Dep 10 --> 11
                oow7 = (87962:87965)';      % Dep 11 --> 12
                oow8 = (104299:104326)';    % Dep 12 --> 13
                oow9 = (113499:113501)';    % Dep 13 --> 14
                oow10 = (122403:122405)';   % Dep 14 --> 15
                oow11 = (130652:130661)';   % Dep 15 --> 16
                oow12 = (143717:143719)';   % Dep 16 --> 17
                oow13 = (150787:150799)';   % Dep 17 --> 18
                oow14 = (163889:163911)';   % Dep 18 --> end
                ind_oow = [oow1;oow2;oow3;oow4;oow5;oow6;oow7;oow8;oow9;oow10;oow11;oow12;oow13;oow14];

                % Isolated points during large gap from 7/30/21 (during Dep 1) to start of Dep 2
                ind_dropout1 = (7420:12590)';
                ind_dropout2 = find(dat.depth(1:ind_dep(2)) < 0);
                ind_dropout = unique([ind_dropout1;ind_dropout2]);

            case 'ERDC'
                % Out-of-water points
                oow1 = (34473:34478)'; % Dep 5 --> 6
                oow2 = (51063:51090)'; % Dep 6 --> 7
                oow3 = (60425:60429)'; % Dep 7 --> 8
                oow4 = (67205:67233)'; % Dep 8 --> 9
                oow5 = (71093:71100)'; % Dep 9 --> 10
                oow6 = (78470:78477)'; % Dep 10 --> 11
                oow7 = (87534:87540)'; % Dep 11 --> 12
                oow8 = (103874:103903)'; % Dep 12 --> 13
                oow9 = (113075:113077)'; % Dep 13 --> 14
                oow10 = (121978:121982)'; % Dep 14 --> 15
                oow11 = (130230:130231)'; % Dep 15 --> 16
                oow12 = 143288;         % Dep 16 --> 17
                oow13 = (150357:150360)';   % Dep 17 --> 18
                oow14 = (163450:163470)';   % Dep 18 --> end
                ind_oow = [oow1;oow2;oow3;oow4;oow5;oow6;oow7;oow8;oow9;oow10;oow11;oow12;oow13;oow14];

                % Isolated points during large gap from 7/30/21 (during Dep 1) to start of Dep 2
                ind_dropout1 = (7420:12589)';
                ind_dropout2 = find(dat.depth(1:ind_dep(1)) < 0);
                ind_dropout = unique([ind_dropout1;ind_dropout2]);
        end
        
    case 'North'
        switch sonde
            case 'BC'
                % Out-of-water points
                oow1 = 8767;                % Dep 5 --> 6
                oow2 = (25818:25833)';      % Dep 6 --> 7
                oow3 = (34655:34673)';      % Dep 7 --> 8
                oow4 = (41585:41594)';      % Dep 8 --> 9
                oow5 = (45507:45509)';      % Dep 9 --> 10
                oow6 = (52681:52689)';      % Dep 10 --> 11
                oow7 = (61893:61928)';      % Dep 11 --> 12
                oow8 = (78392:78394)';      % Dep 12 --> 13
                oow9 = (87464:87466)';      % Dep 13 --> 14
                oow10 = (96663:96669)';     % Dep 14 --> 15
                oow11 = (104890:104899)';   % Dep 15 --> 16
                oow12 = (117837:117842)';   % Dep 16 --> 17
                oow13 = (125024:125025)';   % Dep 17 --> 18
                oow14 = (137984:138002)';   % Dep 18 --> end
                ind_oow = [oow1;oow2;oow3;oow4;oow5;oow6;oow7;oow8;oow9;oow10;oow11;oow12;oow13;oow14];

                % Isolated points during Deployment 2
                ind_dropout1 = (7250:7280)';
                ind_dropout2 = find(dat.depth(1:ind_dep(2)) < 0);
                ind_dropout = unique([ind_dropout1;ind_dropout2]);

            case 'ERDC'
                % Out-of-water points
                oow1 = (9255:9259)'; % Dep 6 --> 7
                oow2 = (18081:18102)'; % Dep 7 --> 8
                oow3 = (25013:25022)'; % Dep 8 --> 9
                oow4 = (28935:28938)'; % Dep 9 --> 10
                oow5 = (36111:36112)'; % Dep 10 --> 11
                oow6 = (45316:45333)'; % Dep 11 --> 12
                oow7 = (61847:61849)'; % Dep 12 --> 13
                oow8 = (70919:70919)'; % Dep 13 --> 14
                oow9 = (80118:80124)'; % Dep 14 --> 15
                oow10 = (88346:88355)'; % Dep 15 --> 16
                oow11 = (101293:101296)'; % Dep 16 --> 17
                oow12 = (108478:108480)'; % Dep 17 --> 18
                oow13 = (121439:121462)';   % Dep 18 --> end
                ind_oow = [oow1;oow2;oow3;oow4;oow5;oow6;oow7;oow8;oow9;oow10;oow11;oow12;oow13];

                % Isolated points during Deployment 2
                ind_dropout1 = (7250:7280)';
                ind_dropout2 = find(dat.depth(1:ind_dep(2)) < 0);
                ind_dropout = unique([ind_dropout1;ind_dropout2]);
        end

    case 'South'
        switch sonde
            case 'BC'
                % Out-of-water points
                oow1 = (36234:36242)';    % Dep 4 --> 5
                oow2 = (49134:49138)';    % Dep 5 --> 6
                oow3 = (52669:52688)';    % Dep 6 --> 7
                oow4 = (61341:61386)';    % Dep 7 --> 8
                oow5 = (68177:68198)';    % Dep 8 --> 9
                oow6 = (72317:72330)';    % Dep 9 --> 10
                oow7 = (79529:79530)';    % Dep 10 --> 11
                oow8 = (82141:82456)';    % Dep 11 --> 12
                oow9 = (96672:96694)';    % Dep 12 --> 13
                oow10 = (105582:105590)'; % Dep 13 --> 14
                oow11 = (114699:114713)'; % Dep 14 --> 16
                oow12 = (127949:127969)'; % Dep 16 --> 17
                oow13 = (134757:134758)'; % Dep 17 --> 18 
                oow14 = (147838:147850)';   % Dep 18 --> end
                ind_oow = [oow1;oow2;oow3;oow4;oow5;oow6;oow7;oow8;oow9;oow10;oow11;oow12;oow13;oow14];


                % Isolated points during large gap from 7/30/21 (during Dep 1) to start of Dep 2
                ind_dropout1 = (18395:18419)';
                ind_dropout2 = find(dat.depth(1:ind_dep(2)) < 0);
                ind_dropout = unique([ind_dropout1;ind_dropout2]);

            case 'ERDC'
                % Out-of-water points
                oow1 = (19871:19873)'; % Dep 2 --> 5
                oow2 = (32765:32791)'; % Dep 5 --> 7
                oow3 = (41444:41494)'; % Dep 7 --> 8
                oow4 = (48285:48299)'; % Dep 8 --> 9
                oow5 = (52422:52433)'; % Dep 9 --> 10
                oow6 = (57630:57631)'; % Dep 10 --> 11
                oow7 = (64748:65013)'; % Dep 11 --> 12
                oow8 = (79228:79235)'; % Dep 12 --> 13
                oow9 = (88098:88114)'; % Dep 13 --> 14
                oow10 = (97223:97237)'; % Dep 14 --> 15
                oow11 = (105534:105537)'; % Dep 15 --> 16
                oow12 = (118263:118268)'; % Dep 16 --> 17
                oow13 = (125056:125058)'; % Dep 17 --> 18
                oow14 = (138138:138149)';   % Dep 18 --> end
                ind_oow = [oow1;oow2;oow3;oow4;oow5;oow6;oow7;oow8;oow9;oow10;oow11;oow12;oow13;oow14];

                % Isolated points during Deployment 2
                ind_dropout1 = (18395:18419)';
                ind_dropout2 = find(dat.depth(1:ind_dep(2)) < 0);
                ind_dropout = unique([ind_dropout1;ind_dropout2]);
        end
end
        
clearvars oow1 oow2 oow3 oow4 oow5 oow6 oow7 oow8 oow9 oow10 oow11 oow12 oow13 oow14

ind_manual = [ind_dropout;ind_oow];

%====DEPTH & PRESSURE======================================================
depth_orig = dat.depth;  % Preserve original depth data for plotting
p_orig = dat.p;     % Preserve original p data for plotting

% (1) Gross Range Test
% INPUTS
low_threshold = 0;      % Lower limit (m) - flag all negative values as FAIL
high_threshold = 5;     % Upper limit (m)

ind_low = find(dat.depth < low_threshold);
ind_high = find(dat.depth > high_threshold);
ind_grossRange = [ind_low;ind_high];
% Don't duplicate points already defined to be manually removed
ind_grossRange = setdiff(ind_grossRange,ind_manual);

% (2) Spike Test
% INPUTS
spike_threshold = 0.1;    % FAIL threshold (m)

d_depth = zeros(length(dat.depth),1);

for i = 2:length(dat.depth)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (dat.depth(i+1)+dat.depth(i-1))/2; % Calculate the average of i+1 and i-1
        d_depth(i) = dat.depth(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_depth) >= spike_threshold);
[~,~,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of depth differences
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
histogram(d_depth)
title([site,' ',sonde,' - Spike Test Histogram'])
xlabel('d_{depth} (m)')
ylabel('Frequency')
set(gca,'YScale','log','FontSize',fontsize)

% Plot original data with flagged points
fig2 = figure(2);clf
fig2.WindowState = 'maximized';
yyaxis left
plot(dat.datetime_utc,depth_orig,'.k','MarkerSize',dotsize,'DisplayName','Original Data');
hold on
plot(dat.datetime_utc(ind_dropout),depth_orig(ind_dropout),'og','MarkerSize',circlesize,'DisplayName','Dropouts')
plot(dat.datetime_utc(ind_oow),depth_orig(ind_oow),'oc','MarkerSize',circlesize,'DisplayName','Out of Water')
plot(dat.datetime_utc(ind_grossRange),depth_orig(ind_grossRange),'or','MarkerSize',circlesize,'DisplayName','Gross Range Test')
plot(dat.datetime_utc(ind_spike),depth_orig(ind_spike),'om','MarkerSize',circlesize,'DisplayName','Spike Test')
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
ylabel('Depth (m)')
yyaxis right
plot(dat.datetime_utc,p_orig,'.','Color',rgb('darkgrey'),'MarkerSize',dotsize,'HandleVisibility','off')
hold on
plot(dat.datetime_utc(ind_dropout),p_orig(ind_dropout),'og','MarkerSize',circlesize,'HandleVisibility','off')
plot(dat.datetime_utc(ind_oow),p_orig(ind_oow),'oc','MarkerSize',circlesize,'HandleVisibility','off')
plot(dat.datetime_utc(ind_grossRange),p_orig(ind_grossRange),'or','MarkerSize',circlesize,'HandleVisibility','off')
plot(dat.datetime_utc(ind_spike),p_orig(ind_spike),'om','MarkerSize',circlesize,'HandleVisibility','off')
ylabel('Pressure (psi)')
xlabel('UTC')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = rgb('darkgrey');
legend('show','location','best')
title([site,' ',sonde,' - Flagged Points from Initial Data QC'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% Create structure to save flagged depth indices
depth_flags = struct('ind_dropout',ind_dropout,'ind_oow',ind_oow,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'high_threshold',high_threshold,'low_threshold',low_threshold,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike is

%====TEMPERATURE===========================================================
T_orig = dat.temperature;  % Preserve original T data for plotting

% (1) Gross Range Test
% INPUTS
low_threshold = -2;      % Lower limit (oC) - Freezing pt of SW at 35 ppt (http://www.csgnetwork.com/h2ofreezecalc.html)
high_threshold = 35;     % Upper limit (oC)

ind_low = find(T_orig < low_threshold);
ind_high = find(T_orig > high_threshold);
ind_grossRange = [ind_low;ind_high];
% Don't duplicate points already defined to be manually removed
ind_grossRange = setdiff(ind_grossRange,ind_manual);

% (2) Spike Test
% INPUTS
spike_threshold = 1;    % FAIL threshold (deg C)

d_T = zeros(length(dat.temperature),1);

for i = 2:length(dat.temperature)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (dat.temperature(i+1)+dat.temperature(i-1))/2; % Calculate the average of i+1 and i-1
        d_T(i) = dat.temperature(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_T) >= spike_threshold);
[C,ig,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of temperature differences
fig3 = figure(3);clf
fig3.WindowState = 'maximized';
histogram(d_T)
title([site,' ',sonde,' - Spike Test Histogram'])
xlabel('d_{temperature} (^oC)')
ylabel('Frequency')
set(gca,'YScale','log','FontSize',fontsize)

% Plot original data with flagged points
fig4 = figure(4);clf
fig4.WindowState = 'maximized';
plot(dat.datetime_utc,T_orig,'.k','MarkerSize',dotsize,'DisplayName','Original Data');
hold on
plot(dat.datetime_utc(ind_dropout),T_orig(ind_dropout),'og','MarkerSize',circlesize,'DisplayName','Dropouts');
plot(dat.datetime_utc(ind_oow),T_orig(ind_oow),'oc','MarkerSize',circlesize,'DisplayName','Out of Water');
plot(dat.datetime_utc(ind_grossRange),T_orig(ind_grossRange),'or','MarkerSize',circlesize,'DisplayName','Gross Range Test');
plot(dat.datetime_utc(ind_spike),T_orig(ind_spike),'om','MarkerSize',circlesize,'DisplayName','Spike Test');
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
hold off
xlabel('UTC')
ylabel('Temperature (^oC)')
legend('show','location','best')
title([site,' ',sonde,' - Flagged Points from Initial Data QC'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites 
set(gca,'FontSize',fontsize)

% Create structure to save flagged temperature indices
T_flags = struct('ind_dropout',ind_dropout,'ind_oow',ind_oow,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'high_threshold',high_threshold,'low_threshold',low_threshold,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike C ig is

%====DO CONCENTRATION======================================================
DOconc_orig = dat.DO_conc;  % Preserve original DO concentration data for plotting

% (1) Gross Range Test
% INPUTS
low_threshold = 0;      % Lower limit (umol/L)
high_threshold = 650;   % Upper limit (umol/L) - 625 umol/L = 20 mg/L

ind_low = find(DOconc_orig < low_threshold);
ind_high = find(DOconc_orig > high_threshold);
ind_grossRange = [ind_low;ind_high];
% Don't duplicate points already defined to be manually removed
ind_grossRange = setdiff(ind_grossRange,ind_manual);

% (2) Spike Test
% INPUTS
spike_threshold = 50;    % FAIL threshold (umol/L)

d_DO = zeros(length(DOconc_orig),1);

for i = 2:length(DOconc_orig)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (DOconc_orig(i+1)+DOconc_orig(i-1))/2; % Calculate the average of i+1 and i-1
        d_DO(i) = DOconc_orig(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_DO) >= spike_threshold);
[~,~,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of DO concentration differences
fig5 = figure(5);clf
fig5.WindowState = 'maximized';
histogram(d_DO)
title([site,' ',sonde,' - Spike Test Histogram'])
xlabel('d_{DO conc} (\mumol/L)')
ylabel('Frequency')
set(gca,'YScale','log','FontSize',fontsize)

% Plot original data with flagged points
fig6 = figure(6);clf
fig6.WindowState = 'maximized';
plot(dat.datetime_utc,DOconc_orig,'.k','MarkerSize',dotsize,'DisplayName','Original Data');
hold on
plot(dat.datetime_utc(ind_dropout),DOconc_orig(ind_dropout),'og','MarkerSize',circlesize,'DisplayName','Dropouts');
plot(dat.datetime_utc(ind_oow),DOconc_orig(ind_oow),'oc','MarkerSize',circlesize,'DisplayName','Out of Water');
plot(dat.datetime_utc(ind_grossRange),DOconc_orig(ind_grossRange),'or','MarkerSize',circlesize,'DisplayName','Gross Range Test');
plot(dat.datetime_utc(ind_spike),DOconc_orig(ind_spike),'om','MarkerSize',circlesize,'DisplayName','Spike Test');
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
hold off
xlabel('UTC')
ylabel('DO Concentration (\mumol/L)')
legend('show','location','best')
title([site,' ',sonde,' - Flagged Points from Initial Data QC'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites--09876543 
set(gca,'FontSize',fontsize)

% Create structure to save flagged DO concentration indices
DO_flags = struct('ind_dropout',ind_dropout,'ind_oow',ind_oow,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'high_threshold',high_threshold,'low_threshold',low_threshold,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike is

%====SALINITY==============================================================
% For this test, use limits for salinity and 
S_orig = dat.salinity;  % Preserve original S data for plotting

% (1) Gross Range Test
% INPUTS
low_threshold = 10;      % Lower limit (psu)
high_threshold = 50;     % Upper limit (psu)

ind_low = find(S_orig < low_threshold);
ind_high = find(S_orig > high_threshold);
ind_grossRange = [ind_low;ind_high];
% Don't duplicate points already defined to be manually removed
ind_grossRange = setdiff(ind_grossRange,ind_manual);

% (2) Spike Test
% INPUTS
spike_threshold = 10;    % FAIL threshold (psu)

d_salinity = zeros(length(S_orig),1);

for i = 2:length(S_orig)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (S_orig(i+1)+S_orig(i-1))/2; % Calculate the average of i+1 and i-1
        d_salinity(i) = S_orig(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_salinity) >= spike_threshold);
[~,~,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of salinity differences
fig7 = figure(7);clf
fig7.WindowState = 'maximized';
histogram(d_salinity)
title([site,' ',sonde,' - Spike Test Histogram'])
xlabel('d_{salinity} (psu)')
ylabel('Frequency')
set(gca,'YScale','log','FontSize',fontsize)

% Plot original data with flagged points
fig8 = figure(8);clf
fig8.WindowState = 'maximized';
plot(dat.datetime_utc,S_orig,'.k','MarkerSize',dotsize,'DisplayName','Original Data')
hold on
plot(dat.datetime_utc(ind_dropout),S_orig(ind_dropout),'og','MarkerSize',circlesize,'DisplayName','Dropouts')
plot(dat.datetime_utc(ind_oow),S_orig(ind_oow),'oc','MarkerSize',circlesize,'DisplayName','Out of Water')
plot(dat.datetime_utc(ind_grossRange),S_orig(ind_grossRange),'or','MarkerSize',circlesize,'DisplayName','Gross Range Test')
plot(dat.datetime_utc(ind_spike),S_orig(ind_spike),'om','MarkerSize',circlesize,'DisplayName','Spike Test')
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
ylabel('Salinity (psu)')
xlabel('UTC')
legend('show','location','best')
title([site,' ',sonde,' - Flagged Points from Initial Data QC'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% Create structure to save flagged salinity indices
S_flags = struct('ind_dropout',ind_dropout,'ind_oow',ind_oow,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'high_threshold',high_threshold,'low_threshold',low_threshold,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike is

%====pH====================================================================
pH_orig = dat.pH;  % Preserve original  data for plotting

% (1) Gross Range Test
% INPUTS
% low_threshold = 7;      % Lower limit
% high_threshold = 9.5;     % Upper limit
low_threshold = 6;      % Lower limit
high_threshold = 9;     % Upper limit

ind_low = find(pH_orig < low_threshold);
ind_high = find(pH_orig > high_threshold);
ind_grossRange = [ind_low;ind_high];
% Don't duplicate points already defined to be manually removed
ind_grossRange = setdiff(ind_grossRange,ind_manual);

% (2) Spike Test
% INPUTS
spike_threshold = 0.25;    % FAIL threshold

d_pH = zeros(length(pH_orig),1);

for i = 2:length(pH_orig)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (pH_orig(i+1)+pH_orig(i-1))/2; % Calculate the average of i+1 and i-1
        d_pH(i) = pH_orig(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_pH) >= spike_threshold);
[~,~,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of pH differences
fig9 = figure(9);clf
fig9.WindowState = 'maximized';
histogram(d_pH)
title([site,' ',sonde,' - Spike Test Histogram'])
xlabel('d_{pH}')
ylabel('Frequency')
set(gca,'YScale','log','FontSize',fontsize)

% Plot original data with flagged points
fig10 = figure(10);clf
fig10.WindowState = 'maximized';
plot(dat.datetime_utc,pH_orig,'.k','MarkerSize',dotsize,'DisplayName','Original Data');
hold on
plot(dat.datetime_utc(ind_dropout),pH_orig(ind_dropout),'og','MarkerSize',circlesize,'DisplayName','Dropouts');
plot(dat.datetime_utc(ind_oow),pH_orig(ind_oow),'oc','MarkerSize',circlesize,'DisplayName','Out of Water');
plot(dat.datetime_utc(ind_grossRange),pH_orig(ind_grossRange),'or','MarkerSize',circlesize,'DisplayName','Gross Range Test');
plot(dat.datetime_utc(ind_spike),pH_orig(ind_spike),'om','MarkerSize',circlesize,'DisplayName','Spike Test');
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
hold off
xlabel('UTC')
ylabel('pH')
legend('show','location','best')
title([site,' ',sonde,' - Flagged Points from Initial Data QC'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% Create structure to save flagged pH indices
pH_flags = struct('ind_dropout',ind_dropout,'ind_oow',ind_oow,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'high_threshold',high_threshold,'low_threshold',low_threshold,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike is

switch sonde
    case 'BC'
        %====Chl a (BC only)===============================================
        chla_orig = dat.chla;  % Preserve original  data for plotting
        
        % (1) Gross Range Test
        % INPUTS
        low_threshold = 0;      % Lower limit
        high_threshold = 100;   % Upper limit
        % high_threshold = 400;   % Upper limit

        ind_low = find(chla_orig < low_threshold);
        ind_high = find(chla_orig > high_threshold);
        ind_grossRange = [ind_low;ind_high];
        % Don't duplicate points already defined to be manually removed
        ind_grossRange = setdiff(ind_grossRange,ind_manual);

        % (2) Spike Test
        % INPUTS
        % spike_threshold = 150;    % FAIL threshold (RFU)
        spike_threshold = 10;    % FAIL threshold (RFU)

        d_chla = zeros(length(dat.chla),1);

        for i = 2:length(dat.chla)-1
            % Only check points if both they and their neighbors have not already been flagged for manual removal
            if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
                ref = (dat.chla(i+1)+dat.chla(i-1))/2; % Calculate the average of i+1 and i-1
                d_chla(i) = dat.chla(i) - ref; % Calculate the difference between i and ref
            end
        end
        % If i - ref is greater than high threshold, mark as a spike
        ind_spike = find(abs(d_chla) >= spike_threshold);
        [~,~,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
        ind_spike(is) = [];

        % Plot histogram of chla differences
        fig11 = figure(11);clf
        fig11.WindowState = 'maximized';
        histogram(d_chla)
        title([site,' ',sonde,' - Spike Test Histogram'])
        xlabel('d_{chla} (RFU)')
        ylabel('Frequency')
        set(gca,'YScale','log','FontSize',fontsize)

        % Plot original data with flagged points
        fig12 = figure(12);clf
        fig12.WindowState = 'maximized';
        plot(dat.datetime_utc,chla_orig,'.k','MarkerSize',dotsize,'DisplayName','Original Data');
        hold on
        plot(dat.datetime_utc(ind_dropout),chla_orig(ind_dropout),'og','MarkerSize',circlesize,'DisplayName','Dropouts');
        plot(dat.datetime_utc(ind_oow),chla_orig(ind_oow),'oc','MarkerSize',circlesize,'DisplayName','Out of Water');
        plot(dat.datetime_utc(ind_grossRange),chla_orig(ind_grossRange),'or','MarkerSize',circlesize,'DisplayName','Gross Range Test');
        plot(dat.datetime_utc(ind_spike),chla_orig(ind_spike),'om','MarkerSize',circlesize,'DisplayName','Spike Test');
        xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
        hold off
        xlabel('UTC')
        ylabel('Chl a (RFU)')
        legend('show','location','best')
        title([site,' ',sonde,' - Flagged Points from Initial Data QC'])
        xlim([dt1 dt2])                 % Use same x limits for comparing sites
        set(gca,'FontSize',fontsize)

        % Create structure to save flagged chla indices
        chla_flags = struct('ind_dropout',ind_dropout,'ind_oow',ind_oow,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'high_threshold',high_threshold,'low_threshold',low_threshold,'spike_threshold',spike_threshold);

        clearvars ind_high ind_low ind_grossRange ind_spike is

    case 'ERDC'
        %====Turbidity (ERDC only)===============================================
        turbidity_orig = dat.turbidity;  % Preserve original  data for plotting
        
        % (1) Gross Range Test
        % INPUTS
        low_threshold = 0;      % Lower limit
        % high_threshold = 1000;  % Upper limit
        high_threshold = 350;  % Upper limit

        ind_low = find(turbidity_orig < low_threshold);
        ind_high = find(turbidity_orig > high_threshold);
        ind_grossRange = [ind_low;ind_high];
        % Don't duplicate points already defined to be manually removed
        ind_grossRange = setdiff(ind_grossRange,ind_manual);

        % (2) Spike Test
        % INPUTS
        % spike_threshold = 2000;    % FAIL threshold (NTU)
        spike_threshold = 50;

        d_turbidity = zeros(length(dat.turbidity),1);

        for i = 2:length(dat.turbidity)-1
            % Only check points if both they and their neighbors have not already been flagged for manual removal
            if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
                ref = (dat.turbidity(i+1)+dat.turbidity(i-1))/2; % Calculate the average of i+1 and i-1
                d_turbidity(i) = dat.turbidity(i) - ref; % Calculate the difference between i and ref
            end
        end
        % If i - ref is greater than high threshold, mark as a spike
        ind_spike = find(abs(d_turbidity) >= spike_threshold);
        [~,~,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
        ind_spike(is) = [];

        % Plot histogram of turbidity differences
        fig11 = figure(11);clf
        fig11.WindowState = 'maximized';
        histogram(d_turbidity)
        title([site,' ',sonde,' - Spike Test Histogram'])
        xlabel('d_{turbidity} (NTU)')
        ylabel('Frequency')
        set(gca,'YScale','log','FontSize',fontsize)

        % Plot original data with flagged points
        fig12 = figure(12);clf
        fig12.WindowState = 'maximized';
        plot(dat.datetime_utc,turbidity_orig,'.k','MarkerSize',dotsize,'DisplayName','Original Data');
        hold on
        plot(dat.datetime_utc(ind_dropout),turbidity_orig(ind_dropout),'og','MarkerSize',circlesize,'DisplayName','Dropouts');
        plot(dat.datetime_utc(ind_oow),turbidity_orig(ind_oow),'oc','MarkerSize',circlesize,'DisplayName','Out of Water');
        plot(dat.datetime_utc(ind_grossRange),turbidity_orig(ind_grossRange),'or','MarkerSize',circlesize,'DisplayName','Gross Range Test');
        plot(dat.datetime_utc(ind_spike),turbidity_orig(ind_spike),'om','MarkerSize',circlesize,'DisplayName','Spike Test');
        xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
        hold off
        xlabel('UTC')
        ylabel('Turbidity (NTU)')
        legend('show','location','best')
        title([site,' ',sonde,' - Flagged Points from Initial Data QC'])
        xlim([dt1 dt2])                 % Use same x limits for comparing sites
        set(gca,'FontSize',fontsize)

        % Create structure to save flagged turbidity indices
        turbidity_flags = struct('ind_dropout',ind_dropout,'ind_oow',ind_oow,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'high_threshold',high_threshold,'low_threshold',low_threshold,'spike_threshold',spike_threshold);

        clearvars ind_high ind_low ind_grossRange ind_spike is
end

%====Percentage of points removed in each test=============================
% length(ind_dropout)/length(depth_orig)*100;
% length(ind_oow)/length(depth_orig)*100;

grossRange.depth = length(depth_flags.ind_grossRange)/length(depth_orig)*100;
spike.depth = length(depth_flags.ind_spike)/length(depth_orig)*100;

grossRange.T = length(T_flags.ind_grossRange)/length(T_orig)*100; 
spike.T = length(T_flags.ind_spike)/length(T_orig)*100;

grossRange.DOconc = length(DO_flags.ind_grossRange)/length(DOconc_orig)*100;
spike.DOconc = length(DO_flags.ind_spike)/length(DOconc_orig)*100;

grossRange.S = length(S_flags.ind_grossRange)/length(S_orig)*100;
spike.S = length(S_flags.ind_spike)/length(S_orig)*100;

grossRange.pH = length(pH_flags.ind_grossRange)/length(pH_orig)*100;
spike.pH = length(pH_flags.ind_spike)/length(pH_orig)*100;

switch sonde
    case 'BC'
        grossRange.chla = length(chla_flags.ind_grossRange)/length(chla_orig)*100;
        spike.chla = length(chla_flags.ind_spike)/length(chla_orig)*100;
    case 'ERDC'
        grossRange.turbidity = length(turbidity_flags.ind_grossRange)/length(turbidity_orig)*100;
        spike.turbidity = length(turbidity_flags.ind_spike)/length(turbidity_orig)*100;
end

%====Clean data============================================================
% Discard points marked for Manual Removal and that failed Gross Range and Spike Tests
dat.depth(depth_flags.ind_dropout) = NaN;
dat.depth(depth_flags.ind_oow) = NaN;
dat.depth(depth_flags.ind_grossRange) = NaN;
dat.depth(depth_flags.ind_spike) = NaN;

dat.p(depth_flags.ind_dropout) = NaN;
dat.p(depth_flags.ind_oow) = NaN;
dat.p(depth_flags.ind_grossRange) = NaN;
dat.p(depth_flags.ind_spike) = NaN;

dat.temperature(T_flags.ind_dropout) = NaN;
dat.temperature(T_flags.ind_oow) = NaN;
dat.temperature(T_flags.ind_grossRange) = NaN;
dat.temperature(T_flags.ind_spike) = NaN;

dat.DO_conc(DO_flags.ind_dropout) = NaN;
dat.DO_conc(DO_flags.ind_oow) = NaN;
dat.DO_conc(DO_flags.ind_grossRange) = NaN;
dat.DO_conc(DO_flags.ind_spike) = NaN;

dat.DO_sat(DO_flags.ind_dropout) = NaN;
dat.DO_sat(DO_flags.ind_oow) = NaN;
dat.DO_sat(DO_flags.ind_grossRange) = NaN;
dat.DO_sat(DO_flags.ind_spike) = NaN;

dat.salinity(S_flags.ind_dropout) = NaN;
dat.salinity(S_flags.ind_oow) = NaN;
dat.salinity(S_flags.ind_grossRange) = NaN;
dat.salinity(S_flags.ind_spike) = NaN;

% dat.actual_cond(S_flags.ind_dropout) = NaN;
% dat.actual_cond(S_flags.ind_oow) = NaN;
% dat.actual_cond(S_flags.ind_grossRange) = NaN;
% dat.actual_cond(S_flags.ind_spike) = NaN;

dat.pH(pH_flags.ind_dropout) = NaN;
dat.pH(pH_flags.ind_oow) = NaN;
dat.pH(pH_flags.ind_grossRange) = NaN;
dat.pH(pH_flags.ind_spike) = NaN;

switch sonde
    case 'BC'
        dat.chla(chla_flags.ind_dropout) = NaN;
        dat.chla(chla_flags.ind_oow) = NaN;
        dat.chla(chla_flags.ind_grossRange) = NaN;
        dat.chla(chla_flags.ind_spike) = NaN;
    case 'ERDC'
        dat.turbidity(turbidity_flags.ind_dropout) = NaN;
        dat.turbidity(turbidity_flags.ind_oow) = NaN;
        dat.turbidity(turbidity_flags.ind_grossRange) = NaN;
        dat.turbidity(turbidity_flags.ind_spike) = NaN;
end

% Flag all NaNs as FAIL
ind_fail = find(isnan(dat.depth));
flags.depth(ind_fail) = 4;
flags.p(ind_fail) = 4;

ind_fail = find(isnan(dat.temperature));
flags.temperature(ind_fail) = 4;

ind_fail = find(isnan(dat.DO_conc));
flags.DO_conc(ind_fail) = 4;
flags.DO_sat(ind_fail) = 4;

ind_fail = find(isnan(dat.salinity));
flags.salinity(ind_fail) = 4;
flags.actual_cond(ind_fail) = 4;

ind_fail = find(isnan(dat.pH));
flags.pH(ind_fail) = 4;

switch sonde
    case 'BC'
        ind_fail = find(isnan(dat.chla));
        flags.chla(ind_fail) = 4;
    case 'ERDC'
        ind_fail = find(isnan(dat.turbidity));
        flags.turbidity(ind_fail) = 4;
end

dat(ind_oow,:) = [];

%====Retime data to 10 min intervals=======================================
dat_TT = table2timetable(dat,'RowTimes','datetime_utc');
newTimeStep = minutes(10);
dat_TT = retime(dat_TT,'regular','mean','TimeStep',newTimeStep); % Calculate mean of values in each time bin

dat_TT.datetime_utc = dat_TT.datetime_utc + newTimeStep/2;
dat_TT.datetime_local = dat_TT.datetime_utc;
dat_TT.datetime_local.TimeZone = 'America/New_York';

dat_TT = rmmissing(dat_TT,'DataVariables',"deployment");

% Find new indices of deployment changes
ind_dep = find(diff(dat_TT.deployment) > 0);

%====Plot the cleaned data after initial QC tests=======================
% Depth & Pressure
fig13 = figure(13);clf
fig13.WindowState = 'maximized';
plot(dat_TT.datetime_utc,dat_TT.depth,'.k','MarkerSize',dotsize);
xline([dat_TT.datetime_utc(1); dat_TT.datetime_utc(ind_dep+1)],'--',label);
ylabel('Depth (m)')
yyaxis right
plot(dat_TT.datetime_utc,dat_TT.p,'.','Color',rgb('darkgrey'),'MarkerSize',dotsize)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = rgb('darkgrey');
ylabel('Pressure (psi)')
title([site,' ',sonde,' - After Initial Data QC'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
% ylim([-.5 3])
set(gca,'FontSize',fontsize)

% Temperature
fig14 = figure(14);clf
fig14.WindowState = 'maximized';
plot(dat_TT.datetime_utc,dat_TT.temperature,'.k','MarkerSize',dotsize);
xline([dat_TT.datetime_utc(1); dat_TT.datetime_utc(ind_dep+1)],'--',label);
ylabel('Temperature (^oC)')
title([site,' ',sonde,' - After Initial Data QC'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% DO Concentration
fig15 = figure(15);clf
fig15.WindowState = 'maximized';
plot(dat_TT.datetime_utc,dat_TT.DO_conc,'.k','MarkerSize',dotsize);
xline([dat_TT.datetime_utc(1); dat_TT.datetime_utc(ind_dep+1)],'--',label);
ylabel('DO Concentration (\mumol/L)')
title([site,' ',sonde,' - After Initial Data QC'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% Salinity
fig16 = figure(16);clf
fig16.WindowState = 'maximized';
plot(dat_TT.datetime_utc,dat_TT.salinity,'.k','MarkerSize',dotsize);
xline([dat_TT.datetime_utc(1); dat_TT.datetime_utc(ind_dep+1)],'--',label);
ylabel('Salinity (psu)')
title([site,' ',sonde,' - After Initial Data QC'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
% ylim([-.5 50])
set(gca,'FontSize',fontsize)

% pH
fig17 = figure(17);clf
fig17.WindowState = 'maximized';
plot(dat_TT.datetime_utc,dat_TT.pH,'.k','MarkerSize',dotsize);
xline([dat_TT.datetime_utc(1); dat_TT.datetime_utc(ind_dep+1)],'--',label);
ylabel('pH')
title([site,' ',sonde,' - After Initial Data QC'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

switch sonde
    case 'BC'
        % Chl a
        fig18 = figure(18);clf
        fig18.WindowState = 'maximized';
        plot(dat_TT.datetime_utc,dat_TT.chla,'.k','MarkerSize',dotsize);
        xline([dat_TT.datetime_utc(1); dat_TT.datetime_utc(ind_dep+1)],'--',label);
        ylabel('Chl a (RFU)')
        title([site,' ',sonde,' - After Initial Data QC'])
        xlim([dt1 dt2])                 % Use same x limits for comparing sites
        set(gca,'FontSize',fontsize)

    case 'ERDC'
        % Turbidity
        fig18 = figure(18);clf
        fig18.WindowState = 'maximized';
        plot(dat_TT.datetime_utc,dat_TT.turbidity,'.k','MarkerSize',dotsize);
        xline([dat_TT.datetime_utc(1); dat_TT.datetime_utc(ind_dep+1)],'--',label);
        ylabel('Turbidity (NTU)')
        title([site,' ',sonde,' - After Initial Data QC'])
        xlim([dt1 dt2])                 % Use same x limits for comparing sites
        set(gca,'FontSize',fontsize)
end

%====Save the cleaned data=================================================
option = questdlg('Save cleaned data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\',site,'\cleaned\initial-qc'])
        switch sonde
            case 'BC'
                sonde1_cleaned = dat_TT;
                save([site,'-bc-cleaned.mat'],'sonde1_cleaned')
            case 'ERDC'
                sonde2_cleaned = dat_TT;
                save([site,'-erdc-cleaned.mat'],'sonde2_cleaned')
        end
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end

%====Save the flagged and cleaned plots====================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        switch sonde
            case 'BC'
                saveFilePath = ['figures\open-water-platform\',site,'\data-qc\bc\initial-qc'];
            case 'ERDC'
                saveFilePath = ['figures\open-water-platform\',site,'\data-qc\erdc\initial-qc'];
        end
        cd([rootpath,saveFilePath])
        saveas(fig1,'depth_spike-histogram.png')
        saveas(fig1,'depth_spike-histogram.fig')
        saveas(fig2,'depth&p_flagged.png')
        saveas(fig2,'depth&p_flagged.fig')
        saveas(fig3,'T_spike-histogram.png')
        saveas(fig3,'T_spike-histogram.fig')
        saveas(fig4,'T_flagged.png')
        saveas(fig4,'T_flagged.fig')
        saveas(fig5,'DOconc_spike-histogram.png')
        saveas(fig5,'DOconc_spike-histogram.fig')
        saveas(fig6,'DOconc_flagged.png')
        saveas(fig6,'DOconc_flagged.fig')
        saveas(fig7,'S_spike-histogram.png')
        saveas(fig7,'S_spike-histogram.fig')
        saveas(fig8,'S_flagged.png')
        saveas(fig8,'S_flagged.fig')
        saveas(fig9,'pH_spike-histogram.png')
        saveas(fig9,'pH_spike-histogram.fig')
        saveas(fig10,'pH_flagged.png')
        saveas(fig10,'pH_flagged.fig')
        switch sonde
            case 'BC'
                saveas(fig11,'chla_spike-histogram.png')
                saveas(fig11,'chla_spike-histogram.fig')
                saveas(fig12,'chla_flagged.png')
                saveas(fig12,'chla_flagged.fig')
            case 'ERDC'
                saveas(fig11,'turbidity_spike-histogram.png')
                saveas(fig11,'turbidity_spike-histogram.fig')
                saveas(fig12,'turbidity_flagged.png')
                saveas(fig12,'turbidity_flagged.fig')
        end
        saveas(fig13,'depth&p_cleaned&retimed.png')
        saveas(fig13,'depth&p_cleaned&retimed.fig')
        saveas(fig14,'T_cleaned&retimed.png')
        saveas(fig14,'T_cleaned&retimed.fig')
        saveas(fig15,'DOconc_cleaned&retimed.png')
        saveas(fig15,'DOconc_cleaned&retimed.fig')
        saveas(fig16,'S_cleaned&retimed.png')
        saveas(fig16,'S_cleaned&retimed.fig')
        saveas(fig17,'pH_cleaned&retimed.png')
        saveas(fig17,'pH_cleaned&retimed.fig')
        switch sonde
            case 'BC'
                saveas(fig18,'chla_cleaned&retimed.png')
                saveas(fig18,'chla_cleaned&retimed.fig')
            case 'ERDC'
                saveas(fig18,'turbidity_cleaned&retimed.png')
                saveas(fig18,'turbidity_cleaned&retimed.fig')
        end
        
        disp('Plots saved!')
    
    case 'No'
        disp('Plots not saved.')
end