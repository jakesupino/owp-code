%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotTimeseries.m
% This script plots daily averages of the measured parameters and
% calculated metabolic rates by site and year, to compare site and 
% interannual variability. 
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 11/25/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import final QC'd data & diel analysis results========================
fn = {'north','gull','south'};

for i = 1:length(fn)
    cd([rootpath,'open-water-platform-data\',fn{i},'\cleaned\final-qc'])
    load('finalQC.mat');
    params.(fn{i}) = finalQC;
    cd([rootpath,'diel-method\uncertainty-analysis\',fn{i}])
    load('MonteCarloResults')
    metab.(fn{i}) = diel_dtd_MC;
end
clearvars finalQC diel_dtd diel_obs

% Find daily tidal range for each site
gull_range = groupsummary(params.gull,"datetime_utc","day","range","depth");
gull_range(find(gull_range.GroupCount < 144),:) = []; % Remove days that don't have all depth measurements (full day = 144 points)
gull_range.day_datetime_utc = datetime(cellstr(gull_range.day_datetime_utc),'TimeZone','UTC');
gull_range = table2timetable(gull_range);
gull_range.GroupCount = [];
gull_range.Properties.VariableNames = {'tidal'};

north_range = groupsummary(params.north,"datetime_utc","day","range","depth");
north_range(find(north_range.GroupCount < 144),:) = []; % Remove days that don't have all depth measurements (full day = 144 points)
north_range.day_datetime_utc = datetime(cellstr(north_range.day_datetime_utc),'TimeZone','UTC');
north_range = table2timetable(north_range);
north_range.GroupCount = [];
north_range.Properties.VariableNames = {'tidal'};

south_range = groupsummary(params.south,"datetime_utc","day","range","depth");
south_range(find(south_range.GroupCount < 144),:) = []; % Remove days that don't have all depth measurements (full day = 144 points)
south_range.day_datetime_utc = datetime(cellstr(south_range.day_datetime_utc),'TimeZone','UTC');
south_range = table2timetable(south_range); 
south_range.GroupCount = [];
south_range.Properties.VariableNames = {'tidal'};

%====Import windspeed and PAR data=========================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

% cd([rootpath,'physical-data\final-dataset'])
% load('par.mat')

%====Calculate daily means for each site===================================
gull_dailyAvg = retime(params.gull,'daily','mean');
north_dailyAvg = retime(params.north,'daily','mean');
south_dailyAvg = retime(params.south,'daily','mean');

wspd_dailyAvg = retime(era5Dat,'daily','mean');

%====Create time tables for each site======================================
%----Gull------------------------------------------------------------------
dt2 = dateshift(metab.gull.daystart_dt,'start','day');
metab.gull.daystart_dt = dt2;
daily.gull = synchronize(gull_dailyAvg,gull_range,wspd_dailyAvg,metab.gull);
daily.gull = removevars(daily.gull,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});
daily.gull(daily.gull.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily.gull(daily.gull.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Remove timepoints with anomalous GPP and ER values
anomER = find(daily.gull.ER_avg > 0);
anomGPP = find(daily.gull.GPP_avg < 0);
daily.gull([anomER;anomGPP],:) = [];

%----North-----------------------------------------------------------------
dt2 = dateshift(metab.north.daystart_dt,'start','day');
metab.north.daystart_dt = dt2;
daily.north = synchronize(north_dailyAvg,north_range,wspd_dailyAvg,metab.north);
daily.north = removevars(daily.north,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});
daily.north(daily.north.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily.north(daily.north.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Remove timepoints with anomalous GPP and ER values
anomER = find(daily.north.ER_avg > 0);
anomGPP = find(daily.north.GPP_avg < 0);
daily.north([anomER;anomGPP],:) = [];

%----South-----------------------------------------------------------------
dt2 = dateshift(metab.south.daystart_dt,'start','day');
metab.south.daystart_dt = dt2;
daily.south = synchronize(south_dailyAvg,south_range,wspd_dailyAvg,metab.south);
daily.south = removevars(daily.south,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});
daily.south(daily.south.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily.south(daily.south.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Remove timepoints with anomalous GPP and ER values
anomER = find(daily.south.ER_avg > 0);
anomGPP = find(daily.south.GPP_avg < 0);
daily.south([anomER;anomGPP],:) = [];

% See triad scheme in https://files-aje-com.s3.amazonaws.com/www/row/_assets/docs/Using_Color_In_Your_Manuscript_Figures.pdf
% and https://www.rapidtables.com/convert/color/rgb-to-hex.html?r=86&g=180&b=233
north_clr = '#5D3A9B';
gull_clr = '#019E73';
south_clr = '#D55E00';

t = datetime(2021,01,01):calyears(1):datetime(2023,01,01);
yr = year(t);

xaxisLbl = t;
xtk = 1:3;

cd([rootpath,'figures'])

%% ====Plot daily means of the parameters==================================
fig1 = figure;clf
fig1.WindowState = 'maximized';
% t1 = tiledlayout(6,1,'tilespacing','tight');
t1 = tiledlayout(8,1,'Padding','compact','TileSpacing','tight');
 
nexttile;
ymin = -5;
ymax = 35;
plot(daily.north.datetime_utc,daily.north.temperature,'-','LineWidth',1.5,'Color',north_clr,'DisplayName','North')
hold on
plot(daily.gull.datetime_utc,daily.gull.temperature,'-.','LineWidth',1.5,'Color',gull_clr,'DisplayName','Gull')
plot(daily.south.datetime_utc,daily.south.temperature,':','LineWidth',1.5,'Color',south_clr,'DisplayName','South')
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'a) Temperature','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('^oC','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on
legend('show','Location','northoutside','Orientation','horizontal')

nexttile;
ymin = 25;
ymax = 40;
plot(daily.north.datetime_utc,daily.north.salinity,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.salinity,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.salinity,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'b) Salinity','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('psu','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = 25;
ymax = 125;
plot(daily.north.datetime_utc,daily.north.DOsat,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.DOsat,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.DOsat,':','LineWidth',1.5,'Color',south_clr)
yline(100,'k','LineWidth',1.5)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'c) DO','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('% sat','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = 7.3;
ymax = 8.7;
plot(daily.north.datetime_utc,daily.north.pH,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.pH,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.pH,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'d) pH','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
% ylabel('-')
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = .5;
ymax = 2.5;
plot(daily.north.datetime_utc,daily.north.tidal,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.tidal,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.tidal,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'e) Tidal range','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('m','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = 0;
ymax = 15;
plot(daily.south.datetime_utc,daily.south.wspd,'-','Color',rgb('darkslategray'),'LineWidth',1)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'f) Wind speed','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('m/s','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = -700;
ymax = 700;
plot(daily.north.datetime_utc,daily.north.GPP_avg,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.GPP_avg,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.GPP_avg,':','LineWidth',1.5,'Color',south_clr)
plot(daily.north.datetime_utc,daily.north.ER_avg,'-','LineWidth',1.5,'Color',north_clr)
plot(daily.gull.datetime_utc,daily.gull.ER_avg,'--','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.ER_avg,':','LineWidth',1.5,'Color',south_clr)
yline(0,'k','LineWidth',1.5)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'g) GPP','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
text(0.02,0.05,'h) ER','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('mmol m^{-2} d^{-1}','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = -500;
ymax = 500;
plot(daily.north.datetime_utc,daily.north.NEM_avg,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.NEM_avg,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.NEM_avg,':','LineWidth',1.5,'Color',south_clr)
yline(0,'k','LineWidth',1.5)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'i) NEM','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('mmol m^{-2} d^{-1}','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

xlabel('Calendar Year')

%====SAVE THE FIGURE WITH THE DESIRED SIZE MANUALLY========================
cd([rootpath,'figures\open-water-platform\all-sites'])

disp('Go to Export Setup, then hit Enter to save .fig and .jpg')
pause

% STEP 1:
% File --> Export Setup --> Width = 7 in; Height = 9 in --> Check "Expand axes to fill figure" --> Apply to Figure --> Export...  

% Step 2:
saveas(fig1,'daily_params&rates.fig')
exportgraphics(fig1,'daily_params&rates.jpg','Resolution',600)

%% ===Supplemental plot w/ chla and turbidity==============================
fig2 = figure(2);clf
fig2.WindowState = 'maximized';
t1 = tiledlayout(8,1,'TileSpacing','tight');
 
nexttile;
ymin = -5;
ymax = 35;
plot(daily.north.datetime_utc,daily.north.temperature,'-','LineWidth',1.5,'Color',north_clr,'DisplayName','North')
hold on
plot(daily.gull.datetime_utc,daily.gull.temperature,'-.','LineWidth',1.5,'Color',gull_clr,'DisplayName','Gull')
plot(daily.south.datetime_utc,daily.south.temperature,':','LineWidth',1.5,'Color',south_clr,'DisplayName','South')
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'a) Temperature','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('^oC','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on
legend('show','Location','northoutside','Orientation','horizontal')

nexttile;
ymin = 25;
ymax = 40;
plot(daily.north.datetime_utc,daily.north.salinity,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.salinity,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.salinity,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'b) Salinity','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('psu','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = 25;
ymax = 125;
plot(daily.north.datetime_utc,daily.north.DOsat,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.DOsat,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.DOsat,':','LineWidth',1.5,'Color',south_clr)
yline(100,'k','LineWidth',1.5)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'c) DO','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('% sat','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = 7.3;
ymax = 8.7;
plot(daily.north.datetime_utc,daily.north.pH,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.pH,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.pH,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'d) pH','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
% ylabel('-')
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = 0;
ymax = 60;
plot(daily.north.datetime_utc,daily.north.chla,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.chla,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.chla,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'e) Chl a','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('RFU','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
plot(daily.north.datetime_utc,daily.north.turbidity,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.turbidity,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.turbidity,':','LineWidth',1.5,'Color',south_clr)
set(gca,'XTickLabel',[])
text(0.02,0.8,'f) Turbidity','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('NTU','FontSize',20)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = .5;
ymax = 2.5;
plot(daily.north.datetime_utc,daily.north.tidal,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.tidal,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.tidal,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'g) Tidal range','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('m','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

nexttile;
ymin = 0;
ymax = 15;
plot(daily.south.datetime_utc,daily.south.wspd,'-','Color',rgb('darkslategray'),'LineWidth',1)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.8,'h) Wind speed','Units','normalized','VerticalAlignment','bottom','FontSize',20,'FontWeight','bold')
ylabel('m/s','FontSize',20)
set(gca,'box','off','LineWidth',1)
grid on

xlabel('Calendar Year')

%====SAVE THE FIGURE WITH THE DESIRED SIZE MANUALLY========================
cd([rootpath,'figures\open-water-platform\all-sites'])

disp('Go to Export Setup, then hit Enter to save .fig and .jpg')
pause

% STEP 1:
% File --> Export Setup --> Width = 7 in; Height = 9 in --> Check "Expand axes to fill figure" --> Apply to Figure

% Step 2:
saveas(fig2,'daily_params-SI.fig')
exportgraphics(fig2,'daily_params-SI.jpg','Resolution',600)