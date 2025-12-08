%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotSeasonalTimescales.m
% This script takes seasonal averages of the measured parameters and
% calculated metabolic rates, and plots all three sites together.
% Within-season variability is plotted as errorbars.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 9/16/2024
% Last updated: 10/5/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

site = {'Gull','North','South'};

for i = 1:length(site)
    %====Import data===========================================================
    % Load the sonde data
    cd([rootpath,'\diel-method\owp-data\final-qc'])
    load([site{i},'_obs.mat']);
    % dat = rmmissing(dat);
    dat_10min.(site{i}) = dat;
    dat_daily.(site{i}) = retime(dat,'daily','mean');

    clear dat

    % Load the MATLAB metabolism results
    cd([rootpath,'diel-method\uncertainty-analysis\',site{i}])
    load('MonteCarloResults');
    metab.(site{i}) = diel_dtd_MC;

    % dt2 = dateshift(metab.(site{i}).daystart_dt,'start','day');
    % metab.(site{i}).daystart_dt = dt2;
    metab.(site{i}) = removevars(metab.(site{i}),{'daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});

    % Remove timepoints with anomalous GPP and ER values
    anomER = find(metab.(site{i}).ER_avg > 0);
    anomGPP = find(metab.(site{i}).GPP_avg < 0);
    metab.(site{i})([anomER;anomGPP],:) = [];

    clear diel_dtd_MC

    %====Combine daily sonde data and metabolic rates into one table=======
    daily.(site{i}) = synchronize(dat_daily.(site{i}),metab.(site{i}));
    daily.(site{i}) = removevars(daily.(site{i}),{'deployment'});
    
    daily = rmmissing(daily);
    monthlyMean.(site{i}) = retime(daily.(site{i}),'monthly','mean');
    monthlyStd.(site{i}) = retime(daily.(site{i}),'monthly',@std);
    
    %====Calculate seasonal means==========================================
    % Metabolism data
    % Get month numbers
    mo = month(daily.(site{i}).datetime_utc);
    % Get year numbers
    yr = year(daily.(site{i}).datetime_utc);

    % Define the seasons by month
    isWinter = ismember(mo,[12 1 2]);
    isSpring = ismember(mo,[3 4 5]);
    isSummer = ismember(mo,[6 7 8]);
    isFall = ismember(mo,[9 10 11]);

    indWinter = find(isWinter);
    indSpring = find(isSpring);
    indSummer = find(isSummer);
    indFall = find(isFall);

    % Add a column for month
    daily.(site{i}).month = mo;

    % Add a column for season number
    daily.(site{i}).season(indWinter) = 1;
    daily.(site{i}).season(indSpring) = 2;
    daily.(site{i}).season(indSummer) = 3;
    daily.(site{i}).season(indFall) = 4;
    daily.(site{i}).season = categorical(daily.(site{i}).season);

    % Add a column for year
    % IF statement so December is included in following winter
    for j = 1:height(daily.(site{i}))
        if daily.(site{i}).month(j) == 12
            daily.(site{i}).year(j) = categorical(yr(j) + 1);
        else
            daily.(site{i}).year(j) = categorical(yr(j));
        end
    end

    % Find rows belonging to each year
    ind_2021 = find(daily.(site{i}).year == '2021');
    ind_2022 = find(daily.(site{i}).year == '2022');
    ind_2023 = find(daily.(site{i}).year == '2023');
    ind_2024 = find(daily.(site{i}).year == '2024');

    % Find seasonal means
    % G = groupsummary(T,groupvars,method,datavars)
    seasonal_means_2021 = groupsummary(daily.(site{i})(ind_2021,:),"season","mean",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP_avg" "ER_avg" "NEM_avg"]);
    seasonal_means_2022 = groupsummary(daily.(site{i})(ind_2022,:),"season","mean",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP_avg" "ER_avg" "NEM_avg"]);
    seasonal_means_2023 = groupsummary(daily.(site{i})(ind_2023,:),"season","mean",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP_avg" "ER_avg" "NEM_avg"]);
    seasonal_means_2024 = groupsummary(daily.(site{i})(ind_2024,:),"season","mean",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP_avg" "ER_avg" "NEM_avg"]);

    seasonal_sd_2021 = groupsummary(daily.(site{i})(ind_2021,:),"season","std",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP_avg" "ER_avg" "NEM_avg"]);
    seasonal_sd_2022 = groupsummary(daily.(site{i})(ind_2022,:),"season","std",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP_avg" "ER_avg" "NEM_avg"]);
    seasonal_sd_2023 = groupsummary(daily.(site{i})(ind_2023,:),"season","std",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP_avg" "ER_avg" "NEM_avg"]);
    seasonal_sd_2024 = groupsummary(daily.(site{i})(ind_2024,:),"season","std",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP_avg" "ER_avg" "NEM_avg"]);

    seasonal_sd_2021(:,1:2) = [];
    seasonal_sd_2022(:,1:2) = [];
    seasonal_sd_2023(:,1:2) = [];
    seasonal_sd_2024(:,1:2) = [];

    seasonal_means_2021.year = repelem(categorical(2021), height(seasonal_means_2021))';
    seasonal_means_2022.year = repelem(categorical(2022), height(seasonal_means_2022))';
    seasonal_means_2023.year = repelem(categorical(2023), height(seasonal_means_2023))';
    seasonal_means_2024.year = repelem(categorical(2024), height(seasonal_means_2024))';

    seasonal.(site{i}) =  [seasonal_means_2021 seasonal_sd_2021; seasonal_means_2022 seasonal_sd_2022; seasonal_means_2023 seasonal_sd_2023; seasonal_means_2024 seasonal_sd_2024];

    seasonal.(site{i}) = movevars(seasonal.(site{i}),"year",'After','season');
end

% For North, need to make a dummy first row since it's missing Summer 2021 data
dummy_row = NaN(1,width(seasonal.North));
dummy_row = array2table(dummy_row);
dummy_row.Properties.VariableNames = seasonal.North.Properties.VariableNames;
dummy_row.season = 3;
dummy_row.year = 2021;
dummy_row.season = categorical(dummy_row.season);
dummy_row.year = categorical(dummy_row.year);
seasonal.North = [dummy_row; seasonal.North];

%====Plot the seasonal means===============================================
cd([rootpath,'figures\open-water-platform\all-sites'])

xaxisLbl = {'Su ''21','F ''21',...
    'W ''22','Sp ''22','Su ''22','F ''22',...
    'W ''23','Sp ''23','Su ''23','F ''23',...
    'W ''24','Sp ''24','Su ''24'};
xtk = 1:13;

% See Wong (2011) and https://www.rapidtables.com/convert/color/rgb-to-hex.html?r=86&g=180&b=233
north_clr = '#5D3A9B';
gull_clr = '#019E73';
south_clr = '#D55E00';

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
t1 = tiledlayout(4,2,'TileSpacing','tight');

ax1 = nexttile;
errorbar(seasonal.North.mean_temperature,seasonal.North.std_temperature,'o-','MarkerSize',8,'LineWidth',2,'Color',north_clr)
hold on
errorbar(seasonal.Gull.mean_temperature,seasonal.Gull.std_temperature,'o-.','MarkerSize',8,'LineWidth',2,'Color',gull_clr)
errorbar(seasonal.South.mean_temperature,seasonal.South.std_temperature,'o:','MarkerSize',8,'LineWidth',2,'Color',south_clr)
set(gca,'XTick',xtk,'XTickLabel',[])
text(0.02,0.8,'Temp','Units','normalized','VerticalAlignment','bottom','FontSize',14,'FontWeight','bold')
ylabel('^oC')

ax2 = nexttile;
errorbar(seasonal.North.mean_pH,seasonal.North.std_pH,'o-','MarkerSize',8,'LineWidth',2,'Color',north_clr)
hold on
errorbar(seasonal.Gull.mean_pH,seasonal.Gull.std_pH,'o-.','MarkerSize',8,'LineWidth',2,'Color',gull_clr)
errorbar(seasonal.South.mean_pH,seasonal.South.std_pH,'o:','MarkerSize',8,'LineWidth',2,'Color',south_clr)
set(gca,'XTick',xtk,'XTickLabel',[])
text(0.02,0.8,'pH','Units','normalized','VerticalAlignment','bottom','FontSize',14,'FontWeight','bold')
% ylabel('-')

ax3 = nexttile;
errorbar(seasonal.North.mean_DOsat,seasonal.North.std_DOsat,'o-','MarkerSize',8,'LineWidth',2,'Color',north_clr)
hold on
errorbar(seasonal.Gull.mean_DOsat,seasonal.Gull.std_DOsat,'o-.','MarkerSize',8,'LineWidth',2,'Color',gull_clr)
errorbar(seasonal.South.mean_DOsat,seasonal.South.std_DOsat,'o:','MarkerSize',8,'LineWidth',2,'Color',south_clr)
set(gca,'XTick',xtk,'XTickLabel',[])
text(0.02,0.8,'DO','Units','normalized','VerticalAlignment','bottom','FontSize',14,'FontWeight','bold')
ylim([50 120])
ylabel('% sat')

ax4 = nexttile;
errorbar(seasonal.North.mean_salinity,seasonal.North.std_salinity,'o-','MarkerSize',8,'LineWidth',2,'Color',north_clr,'DisplayName','North');
hold on
errorbar(seasonal.Gull.mean_salinity,seasonal.Gull.std_salinity,'o-.','MarkerSize',8,'LineWidth',2,'Color',gull_clr,'DisplayName','Gull');
errorbar(seasonal.South.mean_salinity,seasonal.South.std_salinity,'o:','MarkerSize',8,'LineWidth',2,'Color',south_clr,'DisplayName','South');
set(gca,'XTick',xtk,'XTickLabel',[])
legend('show','location','best')
text(0.02,0.8,'Sal','Units','normalized','VerticalAlignment','bottom','FontSize',14,'FontWeight','bold')
ylim([27 38])
ylabel('psu')

ax5 = nexttile;
errorbar(seasonal.North.mean_chla,seasonal.North.std_chla,'o-','MarkerSize',8,'LineWidth',2,'Color',north_clr)
hold on
errorbar(seasonal.Gull.mean_chla,seasonal.Gull.std_chla,'o-.','MarkerSize',8,'LineWidth',2,'Color',gull_clr)
errorbar(seasonal.South.mean_chla,seasonal.South.std_chla,'o:','MarkerSize',8,'LineWidth',2,'Color',south_clr)
set(gca,'XTick',xtk,'XTickLabel',[])
text(0.02,0.8,'Chl a','Units','normalized','VerticalAlignment','bottom','FontSize',14,'FontWeight','bold')
ylim([0 30])
ylabel('RFU')

ax6 = nexttile;
errorbar(seasonal.North.mean_turbidity,seasonal.North.std_turbidity,'o-','MarkerSize',8,'LineWidth',2,'Color',north_clr)
hold on
errorbar(seasonal.Gull.mean_turbidity,seasonal.Gull.std_turbidity,'o-.','MarkerSize',8,'LineWidth',2,'Color',gull_clr)
errorbar(seasonal.South.mean_turbidity,seasonal.South.std_turbidity,'o:','MarkerSize',8,'LineWidth',2,'Color',south_clr)
set(gca,'XTick',xtk,'XTickLabel',[])
text(0.02,0.8,'Turbidity','Units','normalized','VerticalAlignment','bottom','FontSize',14,'FontWeight','bold')
ylim([0 100])
ylabel('NTU')

ax7 = nexttile;
errorbar(seasonal.North.mean_GPP_avg,seasonal.North.std_GPP_avg,'o-','MarkerSize',8,'LineWidth',2,'Color',north_clr)
hold on
errorbar(seasonal.Gull.mean_GPP_avg,seasonal.Gull.std_GPP_avg,'o-.','MarkerSize',8,'LineWidth',2,'Color',gull_clr)
errorbar(seasonal.South.mean_GPP_avg,seasonal.South.std_GPP_avg,'o:','MarkerSize',8,'LineWidth',2,'Color',south_clr)
hold on
errorbar(seasonal.North.mean_ER_avg,seasonal.North.std_ER_avg,'o-','MarkerSize',8,'LineWidth',2,'Color',north_clr)
hold on
errorbar(seasonal.Gull.mean_ER_avg,seasonal.Gull.std_ER_avg,'o-.','MarkerSize',8,'LineWidth',2,'Color',gull_clr)
errorbar(seasonal.South.mean_ER_avg,seasonal.South.std_ER_avg,'o:','MarkerSize',8,'LineWidth',2,'Color',south_clr)
set(gca,'XTick',xtk,'XTickLabel',xaxisLbl)
yline(0,'k','LineWidth',2)
text(0.02,0.8,'ER and GPP','Units','normalized','VerticalAlignment','bottom','FontSize',14,'FontWeight','bold')
ylim([-400 400])
ylabel('mmol O_2 m^{-2} d^{-1}')

ax8 = nexttile;
errorbar(seasonal.North.mean_NEM_avg,seasonal.North.std_NEM_avg,'o-','MarkerSize',8,'LineWidth',2,'Color',north_clr)
hold on
errorbar(seasonal.Gull.mean_NEM_avg,seasonal.Gull.std_NEM_avg,'o-.','MarkerSize',8,'LineWidth',2,'Color',gull_clr)
errorbar(seasonal.South.mean_NEM_avg,seasonal.South.std_NEM_avg,'o:','MarkerSize',8,'LineWidth',2,'Color',south_clr)
set(gca,'XTick',xtk,'XTickLabel',xaxisLbl)
yline(0,'k','LineWidth',2)
text(0.02,0.8,'NEM','Units','normalized','VerticalAlignment','bottom','FontSize',14,'FontWeight','bold')
ylim([-250 250])
ylabel('mmol O_2 m^{-2} d^{-1}')

% title(t1,'Seasonal Means','FontSize',18,'FontWeight','bold')

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'x')
%%
%====Plot the monthly means================================================
fig2 = figure(2);clf
fig2.WindowState = 'maximized';
t2 = tiledlayout(4,2,'TileSpacing','tight');

ax1 = nexttile;
errorbar(monthlyMean.North.datetime_utc,monthlyMean.North.temperature,monthlyStd.North.temperature,'.:','MarkerSize',12,'LineWidth',2,'Color',north_clr)
hold on
errorbar(monthlyMean.Gull.datetime_utc,monthlyMean.Gull.temperature,monthlyStd.Gull.temperature,'.:','MarkerSize',12,'LineWidth',2,'Color',gull_clr)
errorbar(monthlyMean.South.datetime_utc,monthlyMean.South.temperature,monthlyStd.South.temperature,'.:','MarkerSize',12,'LineWidth',2,'Color',south_clr)
set(gca,'XTickLabel',[])
ylabel('T (^oC)')

ax2 = nexttile;
errorbar(monthlyMean.North.datetime_utc,monthlyMean.North.pH,monthlyStd.North.pH,'.:','MarkerSize',12,'LineWidth',2,'Color',north_clr)
hold on
errorbar(monthlyMean.Gull.datetime_utc,monthlyMean.Gull.pH,monthlyStd.Gull.pH,'.:','MarkerSize',12,'LineWidth',2,'Color',gull_clr)
errorbar(monthlyMean.South.datetime_utc,monthlyMean.South.pH,monthlyStd.South.pH,'.:','MarkerSize',12,'LineWidth',2,'Color',south_clr)
set(gca,'XTickLabel',[])
ylabel('pH')

ax3 = nexttile;
errorbar(monthlyMean.North.datetime_utc,monthlyMean.North.DOsat,monthlyStd.North.DOsat,'.:','MarkerSize',12,'LineWidth',2,'Color',north_clr)
hold on
errorbar(monthlyMean.Gull.datetime_utc,monthlyMean.Gull.DOsat,monthlyStd.Gull.DOsat,'.:','MarkerSize',12,'LineWidth',2,'Color',gull_clr)
errorbar(monthlyMean.South.datetime_utc,monthlyMean.South.DOsat,monthlyStd.South.DOsat,'.:','MarkerSize',12,'LineWidth',2,'Color',south_clr)
set(gca,'XTickLabel',[])
ylabel('DO (% sat)')

ax4 = nexttile;
errorbar(monthlyMean.North.datetime_utc,monthlyMean.North.salinity,monthlyStd.North.salinity,'.:','MarkerSize',12,'LineWidth',2,'Color',north_clr,'DisplayName','North')
hold on
errorbar(monthlyMean.Gull.datetime_utc,monthlyMean.Gull.salinity,monthlyStd.Gull.salinity,'.:','MarkerSize',12,'LineWidth',2,'Color',gull_clr,'DisplayName','Gull')
errorbar(monthlyMean.South.datetime_utc,monthlyMean.South.salinity,monthlyStd.South.salinity,'.:','MarkerSize',12,'LineWidth',2,'Color',south_clr,'DisplayName','South')
set(gca,'XTickLabel',[])
ylabel('S (psu)')
legend('show','location','best')

ax5 = nexttile;
errorbar(monthlyMean.North.datetime_utc,monthlyMean.North.chla,monthlyStd.North.chla,'.:','MarkerSize',12,'LineWidth',2,'Color',north_clr)
hold on
errorbar(monthlyMean.Gull.datetime_utc,monthlyMean.Gull.chla,monthlyStd.Gull.chla,'.:','MarkerSize',12,'LineWidth',2,'Color',gull_clr)
errorbar(monthlyMean.South.datetime_utc,monthlyMean.South.chla,monthlyStd.South.chla,'.:','MarkerSize',12,'LineWidth',2,'Color',south_clr)
set(gca,'XTickLabel',[])
ylabel('Chl a (RFU)')

ax6 = nexttile;
errorbar(monthlyMean.North.datetime_utc,monthlyMean.North.turbidity,monthlyStd.North.turbidity,'.:','MarkerSize',12,'LineWidth',2,'Color',north_clr)
hold on
errorbar(monthlyMean.Gull.datetime_utc,monthlyMean.Gull.turbidity,monthlyStd.Gull.turbidity,'.:','MarkerSize',12,'LineWidth',2,'Color',gull_clr)
errorbar(monthlyMean.South.datetime_utc,monthlyMean.South.turbidity,monthlyStd.South.turbidity,'.:','MarkerSize',12,'LineWidth',2,'Color',south_clr)
set(gca,'XTickLabel',[])
ylim([0 100])
ylabel('Turbidity (NTU)')

ax7 = nexttile;
errorbar(monthlyMean.North.datetime_utc,monthlyMean.North.GPP_avg,monthlyStd.North.GPP_avg,'.:','MarkerSize',12,'LineWidth',2,'Color',north_clr)
hold on
errorbar(monthlyMean.Gull.datetime_utc,monthlyMean.Gull.GPP_avg,monthlyStd.Gull.GPP_avg,'.:','MarkerSize',12,'LineWidth',2,'Color',gull_clr)
errorbar(monthlyMean.South.datetime_utc,monthlyMean.South.GPP_avg,monthlyStd.South.GPP_avg,'.:','MarkerSize',12,'LineWidth',2,'Color',south_clr)
hold on
errorbar(monthlyMean.North.datetime_utc,monthlyMean.North.ER_avg,monthlyStd.North.ER_avg,'.:','MarkerSize',12,'LineWidth',2,'Color',north_clr)
errorbar(monthlyMean.Gull.datetime_utc,monthlyMean.Gull.ER_avg,monthlyStd.Gull.ER_avg,'.:','MarkerSize',12,'LineWidth',2,'Color',gull_clr)
errorbar(monthlyMean.South.datetime_utc,monthlyMean.South.ER_avg,monthlyStd.South.ER_avg,'.:','MarkerSize',12,'LineWidth',2,'Color',south_clr)
ylabel('ER and GPP')
ylim([-400 400])
yline(0,'k','LineWidth',2)

ax8 = nexttile;
errorbar(monthlyMean.North.datetime_utc,monthlyMean.North.NEM_avg,monthlyStd.North.NEM_avg,'.:','MarkerSize',12,'LineWidth',2,'Color',north_clr)
hold on
errorbar(monthlyMean.Gull.datetime_utc,monthlyMean.Gull.NEM_avg,monthlyStd.Gull.NEM_avg,'.:','MarkerSize',12,'LineWidth',2,'Color',gull_clr)
errorbar(monthlyMean.South.datetime_utc,monthlyMean.South.NEM_avg,monthlyStd.South.NEM_avg,'.:','MarkerSize',12,'LineWidth',2,'Color',south_clr)
ylabel('NEM')
ylim([-250 250])
yline(0,'k','LineWidth',2)

title(t2,'Monthly Means','FontSize',18,'FontWeight','bold')

%% ====Plots of continuous data============================================
fig = figure(1);clf
fig.WindowState = 'maximized';
t1 = tiledlayout(8,1,'TileSpacing','compact');

ax1 = nexttile;
plot(dat_10min.North.datetime_utc,dat_10min.North.temperature,'-','Color',north_clr)
hold on
plot(dat_10min.Gull.datetime_utc,dat_10min.Gull.temperature,'-','Color',gull_clr)
plot(dat_10min.South.datetime_utc,dat_10min.South.temperature,'-','Color',south_clr)
set(gca,'XTickLabel',[])
ylabel('T (^oC)')

ax2 = nexttile;
plot(dat_10min.North.datetime_utc,dat_10min.North.salinity,'-','Color',north_clr)
hold on
plot(dat_10min.Gull.datetime_utc,dat_10min.Gull.salinity,'-','Color',gull_clr)
plot(dat_10min.South.datetime_utc,dat_10min.South.salinity,'-','Color',south_clr)
set(gca,'XTickLabel',[])
ylabel('S (psu)')

ax3 = nexttile;
plot(dat_10min.North.datetime_utc,dat_10min.North.DOsat,'-','Color',north_clr)
hold on
plot(dat_10min.Gull.datetime_utc,dat_10min.Gull.DOsat,'-','Color',gull_clr)
plot(dat_10min.South.datetime_utc,dat_10min.South.DOsat,'-','Color',south_clr)
set(gca,'XTickLabel',[])
ylabel('DO (% sat)')

ax4 = nexttile;
plot(dat_10min.North.datetime_utc,dat_10min.North.pH,'-','Color',north_clr)
hold on
plot(dat_10min.Gull.datetime_utc,dat_10min.Gull.pH,'-','Color',gull_clr)
plot(dat_10min.South.datetime_utc,dat_10min.South.pH,'-','Color',south_clr)
set(gca,'XTickLabel',[])
ylabel('pH')

ax5 = nexttile;
plot(north_tidal.date,north_tidal.tidal,'-','Color',north_clr)
hold on
plot(gull_tidal.date,gull_tidal.tidal,'-','Color',gull_clr)
plot(south_tidal.date,south_tidal.tidal,'-','Color',south_clr)
set(gca,'XTickLabel',[])
ylabel('Tidal (m)')

ax6 = nexttile;
plot(era5Dat.datetime,era5Dat.wspd,'-')
set(gca,'XTickLabel',[])
ylabel('Wind speed (m/s)')

% ax5 = nexttile;
% plot(dat_10min.North.datetime_utc,dat_10min.North.chla,'-','Color',north_clr)
% hold on
% plot(dat_10min.Gull.datetime_utc,dat_10min.Gull.chla,'-','Color',gull_clr)
% plot(dat_10min.South.datetime_utc,dat_10min.South.chla,'-','Color',south_clr)
% set(gca,'XTickLabel',[])
% ylabel('Chl a (RFU)')
% 
% ax6 = nexttile;
% plot(dat_10min.North.datetime_utc,dat_10min.North.turbidity,'-','Color',north_clr)
% hold on
% plot(dat_10min.Gull.datetime_utc,dat_10min.Gull.turbidity,'-','Color',gull_clr)
% plot(dat_10min.South.datetime_utc,dat_10min.South.turbidity,'-','Color',south_clr)
% set(gca,'XTickLabel',[])
% ylabel('Turbidity (NTU)')

ax7 = nexttile;
plot(metab.North.date,metab.North.GPP_avg,'-','Color',north_clr)
hold on
plot(metab.Gull.date,metab.Gull.GPP_avg,'-','Color',gull_clr)
plot(metab.South.date,metab.South.GPP_avg,'-','Color',south_clr)
set(gca,'XTickLabel',[])
ylabel('Turbidity (NTU)')
ylabel('ER and GPP')

set(gcf, 'Units', 'Inches', 'Position', [-9.9000   -2.6583    8.5000   10.6500], 'PaperUnits', 'Inches', 'PaperSize', [8.5, 11])

xlabel('Calendar Year')

%%
yyaxis right
% plot(daily.north.datetime_utc,daily.north.NEM)
% hold on
plot(daily_gull.datetime_utc,daily_gull.NEM)
% plot(daily.south.datetime_utc,daily.south.NEM)
ylabel('NEM')

nexttile
yyaxis left
plot(daily_gull.datetime_utc,daily_gull.depth)
ylabel('Depth')
yyaxis right
plot(daily_gull.datetime_utc,daily_gull.par)
ylabel('PAR')

title(t1,'Daily Means','FontSize',18,'FontWeight','bold')