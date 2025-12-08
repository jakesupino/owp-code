%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalCuts.m
% This script looks at the turbidity and chla data from all platforms
% together and periods to cut are defined based on visual assessment.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 11/12/2024
% Last updated: 11/20/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

% cd([rootpath,'open-water-platform-data\gull\cleaned\movmed'])
% load('gull-erdc-cleaned.mat');
% gull = sonde2_cleaned;
% 
% cd([rootpath,'open-water-platform-data\north\cleaned\movmed'])
% load('north-erdc-cleaned.mat');
% north = sonde2_cleaned;
% 
% cd([rootpath,'open-water-platform-data\south\cleaned\movmed'])
% load('south-erdc-cleaned.mat');
% south = sonde2_cleaned;
% 
% clearvars sonde2_cleaned;

cd([rootpath,'open-water-platform-data\north\cleaned\combined'])
load('north-combined.mat');
north = sondes_comb;

cd([rootpath,'open-water-platform-data\gull\cleaned\combined'])
load('gull-combined.mat');
gull = sondes_comb;

cd([rootpath,'open-water-platform-data\south\cleaned\combined'])
load('south-combined.mat');
south = sondes_comb;

clearvars sondes_comb

%====Import windspeed data==============================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

% % Find daily tidal range for each site
% gull_range = groupsummary(gull,"datetime_utc","day","range","depth");
% gull_range.day_datetime_utc = datetime(cellstr(gull_range.day_datetime_utc),'TimeZone','UTC');
% gull_range = table2timetable(gull_range);
% gull_range.GroupCount = [];
% gull_range.Properties.VariableNames = {'tidal'};

% Find indices of deployment changes
ind_dep = find(diff(south.deployment) > 0);
dep = table([1;south.deployment(ind_dep+1)],[1;ind_dep+1]);
dep.Properties.VariableNames = {'depNum','ind'};

label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
    'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
    'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
    'Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};

% See Wong (2011) and https://www.rapidtables.com/convert/color/rgb-to-hex.html?r=86&g=180&b=233
north_clr = '#56B4E9';
gull_clr = '#019E73';
south_clr = '#D55E00';

%====TURBIDITY=============================================================
fig = figure;clf
fig.WindowState = 'maximized';
ax = axes;
yyaxis left
plot(era5Dat.datetime,era5Dat.wspd,'-k','HandleVisibility','off')
ylabel('Wind speed (m/s)')
ax.YAxis(1).Color = 'b';
yyaxis right
plot(gull.datetime_utc,gull.turbidity,'.','Color',gull_clr,'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.turbidity,'.','Color',north_clr,'DisplayName','North')
plot(south.datetime_utc,south.turbidity,'.','Color',south_clr,'DisplayName','South')
yline(100,'--r','LineWidth',2,'DisplayName','100 NTU')
yline(300,'--r','LineWidth',2,'DisplayName','300 NTU')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Turbidity (NTU)')
legend('show','location','best')
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
title('After Moving Median Test')

% Retain original turbidity values for comparison
turbidity_orig.north = north.turbidity; 
turbidity_orig.gull = gull.turbidity;
turbidity_orig.south = south.turbidity;

% Remove values that are exactly zero
north.turbidity(find(ismember(north.turbidity,0))) = NaN;
gull.turbidity(find(ismember(gull.turbidity,0))) = NaN;
south.turbidity(find(ismember(south.turbidity,0))) = NaN;

% Remove periods based on visual assessment
%----North-----------------------------------------------------------------
% Find start indices of partial deployment periods to cut
s{1} = '7/21/22 8:45';
s{2} = '8/13/22 16:25';
s{3} = '11/18/22 11:45';
s{4} = '5/26/23 10:55';
s{5} = '7/2/23 12:15';
s{6} = '4/26/24 09:45';
for i = 1:length(s)
    t = datetime(s{i},'InputFormat','MM/dd/yy HH:mm','TimeZone','UTC');
    [~, ind_start(i,1)] = ismember(t,north.datetime_utc);
end
% Find start indices of whole deployments to cut
ind_start(end+1) = find(ismember(north.deployment,2),1,'first');
ind_start(end+1) = find(ismember(north.deployment,15),1,'first');
ind_start(end+1) = find(ismember(north.deployment,16),1,'first');
ind_start(end+1) = find(ismember(north.deployment,17),1,'first');
ind_start = sort(ind_start);

% Find end indices of periods to cut (last index of current deployment)
dep_rm = [2;7;8;10;13;14;15;16;17;18];
for i = 1:length(dep_rm)
    ind_end(i,1) = find(ismember(north.deployment,dep_rm(i)),1,'last');
end

ind_rm = table(ind_start,ind_end);

for i = 1:height(ind_rm)
    north.turbidity(ind_rm.ind_start(i):ind_rm.ind_end(i)) = NaN;
end

clearvars ind_start ind_end ind_rm

%----Gull------------------------------------------------------------------
% Find start indices of partial deployment periods to cut
s{1} = '7/21/21 9:05';
s{2} = '10/23/21 18:45';
s{3} = '5/3/22 8:05';
s{4} = '7/12/22 13:25';
s{5} = '11/11/22 18:35';
s{6} = '6/28/23 8:15';
s{7} = '9/20/23 17:25';
s{8} = '1/8/24 22:55';
s{9} = '3/21/24 12:45';
for i = 1:length(s)
    t = datetime(s{i},'InputFormat','MM/dd/yy HH:mm','TimeZone','UTC');
    [~, ind_start(i,1)] = ismember(t,gull.datetime_utc);
end
% Find start indices of whole deployments to cut
ind_start(end+1) = find(ismember(gull.deployment,5),1,'first');
ind_start(end+1) = find(ismember(gull.deployment,8),1,'first');
ind_start(end+1) = find(ismember(gull.deployment,13),1,'first');
ind_start(end+1) = find(ismember(gull.deployment,17),1,'first');
ind_start = sort(ind_start);

% Find end indices of periods to cut (last index of current deployment)
dep_rm = [1;2;5;6;7;8;10;13;14;15;16;17;18];
for i = 1:length(dep_rm)
    ind_end(i,1) = find(ismember(gull.deployment,dep_rm(i)),1,'last');
end

ind_rm = table(ind_start,ind_end);

for i = 1:height(ind_rm)
    gull.turbidity(ind_rm.ind_start(i):ind_rm.ind_end(i)) = NaN;
end

% Special case: Translate Dep 16 down
ind1 = find(ismember(gull.deployment,16),1,'first');
ind2 = find(ismember(gull.deployment,16),1,'last');
min_val = min(gull.turbidity(ind1:ind2),[],1,'omitmissing');
gull.turbidity(ind1:ind2) = gull.turbidity(ind1:ind2) - min_val;

clearvars s t ind_start ind_end ind_rm

%----South-----------------------------------------------------------------
% Find start indices of partial deployment periods to cut
s{1} = '7/18/21 20:15';
s{2} = '10/22/21 20:45';
s{3} = '7/21/22 6:35';
s{4} = '8/11/22 3:15';
s{5} = '11/07/22 14:35';
s{6} = '3/28/23 7:05';
s{7} = '7/2/23 12:55';
s{8} = '4/26/24 9:05';
for i = 1:length(s)
    t = datetime(s{i},'InputFormat','MM/dd/yy HH:mm','TimeZone','UTC');
    [~, ind_start(i,1)] = ismember(t,south.datetime_utc);
end
% Find start indices of whole deployments to cut
ind_start(end+1) = find(ismember(south.deployment,5),1,'first');
ind_start(end+1) = find(ismember(south.deployment,13),1,'first');
ind_start(end+1) = find(ismember(south.deployment,15),1,'first');
ind_start(end+1) = find(ismember(south.deployment,16),1,'first');
ind_start(end+1) = find(ismember(south.deployment,17),1,'first');
ind_start = sort(ind_start);

% Find end indices of periods to cut (last index of current deployment)
dep_rm = [1;2;5;7;8;10;12;13;14;15;16;17;18];
for i = 1:length(dep_rm)
    ind_end(i,1) = find(ismember(south.deployment,dep_rm(i)),1,'last');
end

ind_rm = table(ind_start,ind_end);

for i = 1:height(ind_rm)
    south.turbidity(ind_rm.ind_start(i):ind_rm.ind_end(i)) = NaN;
end

clearvars s t ind_start ind_end ind_rm

fig = figure;clf
fig.WindowState = 'maximized';
ax = axes;
yyaxis left
plot(era5Dat.datetime,era5Dat.wspd,':k','HandleVisibility','off')
ylabel('Wind speed (m/s)')
ax.YAxis(1).Color = 'b';
yyaxis right
plot(north.datetime_utc,turbidity_orig.north,'o','MarkerSize',5,'LineWidth',.5,'Color',rgb('LightSteelBlue'),'DisplayName','North orig')
hold on
plot(gull.datetime_utc,turbidity_orig.gull,'o','MarkerSize',5,'LineWidth',.5,'Color',rgb('DarkSeaGreen'),'DisplayName','Gull orig')
plot(south.datetime_utc,turbidity_orig.south,'o','MarkerSize',5,'LineWidth',.5,'Color',rgb('Tan'),'DisplayName','South orig')
plot(north.datetime_utc,north.turbidity,'.','Color',north_clr,'DisplayName','North final')
plot(gull.datetime_utc,gull.turbidity,'.','Color',gull_clr,'DisplayName','Gull final')
plot(south.datetime_utc,south.turbidity,'.','Color',south_clr,'DisplayName','South final')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Turbidity (NTU)')
legend('show','location','best')
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
title('After Final Cuts')
%%
%====CHL A=================================================================
fig = figure;clf
fig.WindowState = 'maximized';
ax = axes;
yyaxis left
plot(era5Dat.datetime,era5Dat.wspd,':k','HandleVisibility','off')
ylabel('Wind speed (m/s)')
ax.YAxis(1).Color = 'b';
yyaxis right
plot(gull.datetime_utc,gull.chla,'.','Color',gull_clr,'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.chla,'.','Color',north_clr,'DisplayName','North')
plot(south.datetime_utc,south.chla,'.','Color',south_clr,'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Chl a (RFU)')
legend('show','location','best')
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
title('After Moving Median Test')

% Retain original chl a values for comparison
chla_orig.north = north.chla; 
chla_orig.gull = gull.chla;
chla_orig.south = south.chla;

% Remove values that are exactly zero
north.chla(find(ismember(north.chla,0))) = NaN;
gull.chla(find(ismember(gull.chla,0))) = NaN;
south.chla(find(ismember(south.chla,0))) = NaN;

% Remove periods based on visual assessment
%----North-----------------------------------------------------------------
% Find start indices of partial deployment periods to cut
s{1} = '10/16/22 4:15';
s{2} = '11/2/22 22:05';
for i = 1:length(s)
    t = datetime(s{i},'InputFormat','MM/dd/yy HH:mm','TimeZone','UTC');
    [~, ind_start(i,1)] = ismember(t,north.datetime_utc);
end
ind_start = sort(ind_start);

% Find end indices of periods to cut (last index of current deployment)
dep_rm = [9;10];
for i = 1:length(dep_rm)
    ind_end(i,1) = find(ismember(north.deployment,dep_rm(i)),1,'last');
end

ind_rm = table(ind_start,ind_end);

for i = 1:height(ind_rm)
    north.chla(ind_rm.ind_start(i):ind_rm.ind_end(i)) = NaN;
end

clearvars s t ind_start ind_end ind_rm

%----Gull-----------------------------------------------------------------
% Find start indices of partial deployment periods to cut
s{1} = '6/13/23 14:05';
for i = 1:length(s)
    t = datetime(s{i},'InputFormat','MM/dd/yy HH:mm','TimeZone','UTC');
    [~, ind_start(i,1)] = ismember(t,gull.datetime_utc);
end
ind_start = sort(ind_start);

% Find end indices of periods to cut (last index of current deployment)
dep_rm = [13];
for i = 1:length(dep_rm)
    ind_end(i,1) = find(ismember(gull.deployment,dep_rm(i)),1,'last');
end

ind_rm = table(ind_start,ind_end);

for i = 1:height(ind_rm)
    gull.chla(ind_rm.ind_start(i):ind_rm.ind_end(i)) = NaN;
end

clearvars ind_start ind_end ind_rm

%----South-----------------------------------------------------------------
% Find start indices of partial deployment periods to cut
s{1} = '11/2/22 22:05';
for i = 1:length(s)
    t = datetime(s{i},'InputFormat','MM/dd/yy HH:mm','TimeZone','UTC');
    [~, ind_start(i,1)] = ismember(t,south.datetime_utc);
end
ind_start = sort(ind_start);

% Find end indices of periods to cut (last index of current deployment)
dep_rm = [10];
for i = 1:length(dep_rm)
    ind_end(i,1) = find(ismember(south.deployment,dep_rm(i)),1,'last');
end

ind_rm = table(ind_start,ind_end);

for i = 1:height(ind_rm)
    south.chla(ind_rm.ind_start(i):ind_rm.ind_end(i)) = NaN;
end

clearvars s t ind_start ind_end ind_rm

fig = figure;clf
fig.WindowState = 'maximized';
ax = axes;
yyaxis left
plot(era5Dat.datetime,era5Dat.wspd,':k','HandleVisibility','off')
ylabel('Wind speed (m/s)')
ax.YAxis(1).Color = 'b';
yyaxis right
plot(north.datetime_utc,chla_orig.north,'o','MarkerSize',5,'LineWidth',.5,'Color',rgb('LightSteelBlue'),'DisplayName','North orig')
hold on
plot(gull.datetime_utc,chla_orig.gull,'o','MarkerSize',5,'LineWidth',.5,'Color',rgb('DarkSeaGreen'),'DisplayName','Gull orig')
plot(south.datetime_utc,chla_orig.south,'o','MarkerSize',5,'LineWidth',.5,'Color',rgb('Tan'),'DisplayName','South orig')
plot(north.datetime_utc,north.chla,'.','Color',north_clr,'DisplayName','North final')
plot(gull.datetime_utc,gull.chla,'.','Color',gull_clr,'DisplayName','Gull final')
plot(south.datetime_utc,south.chla,'.','Color',south_clr,'DisplayName','South final')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Chl a (RFU)')
legend('show','location','best')
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
title('After Final Cuts')
%%
%====Save data=============================================================
option = questdlg('Save data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\north\cleaned\final-qc'])
        finalQC = north;
        save('finalQC.mat','finalQC')
        
        cd([rootpath,'open-water-platform-data\gull\cleaned\final-qc'])
        finalQC = gull;
        save('finalQC.mat','finalQC')
     
        cd([rootpath,'open-water-platform-data\south\cleaned\final-qc'])
        finalQC = south;
        save('finalQC.mat','finalQC')

        disp('Files saved!')
    case 'No'
        disp('File not saved.')
end