%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dielAnalysis.m
% This script uses the diel oxygen method to calculate GPP, ER, and NEM
% following Beck et al. (2015)
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 1/3/2024
% Last updated: 2/2/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
cd([rootpath,'open-water-platform-data\gull\cleaned'])
prompt = {'Choose the sonde to analyze'};
answer = questdlg(prompt,'Sonde Selection','Sonde 1','Sonde 2','Cancel','Cancel');

switch answer
    case 'Sonde 1'
        load gull-bc-cleaned.mat
        dat = sonde1_cleaned;
        sondename = 'BC';
    case 'Sonde 2'
        load gull-erdc-cleaned.mat
        dat = sonde2_cleaned;
        sondename = 'ERDC';
end

cd([rootpath,'physical-data\wind-speed'])
load windSpeed.mat

cd([rootpath,'physical-data\par'])
load par.mat

cd([rootpath,'physical-data\baro-pressure'])
load baroPress.mat

% Retime cleaned Gull met data to same datetimes as sonde data 
newTimes = dat.datetime_utc(1):minutes(10):dat.datetime_utc(end);
metDat_rt = retime(metDat_cleaned,newTimes,'mean'); 

%====Gap fill Gull met station wind speed data using NOAA daily means======
% Retime NOAA data to same datetimes as sonde data 
noaaDat_rt = retime(noaaDat,newTimes,'previous');

% Find where there are gaps in Gull met station wind speed data
for i = 1:height(metDat_rt.datetime_utc)
    nid = find(isnan(metDat_rt.wspd));
    id(1:length(nid),1) = nid;
end

figure(1),clf
plot(metDat_rt.datetime_utc,metDat_rt.wspd,'.-')
hold on
plot(noaaDat_rt.date(id),noaaDat_rt.wspd_avg(id),'or','MarkerSize',4)
ylabel('Wind Speed (m/s)')
legend('Gull Met Station','NOAA Daily Mean')
title('Gap-Filled Wind Speed Data')

% Replace missing Gull windspeed data with NOAA windspeed data
ind_nan = find(isnan(metDat_rt.wspd));
metDat_rt.wspd(ind_nan) = noaaDat_rt.wspd_avg(ind_nan);

%====Gap filling Gull met station air T data===============================
% Retime PAR data to same datetimes as sonde data 
parDat_rt = retime(parDat,newTimes,'previous');
% Cut off data once original dataset ends
endDate = parDat.datetime_utc(end);
ind_end = find(ismember(parDat_rt.datetime_utc,endDate));
parDat_rt(ind_end:end,{'Tair' 'light_lux' 'par'}) = {NaN};

% Retime Baro Pressure HOBO data to same datetimes as sonde data
bpDat_rt = retime(bpDat,newTimes,'previous');
% Cut off data once original dataset ends
endDate = bpDat.datetime_utc(end);
ind_end = find(ismember(bpDat_rt.datetime_utc,endDate));
bpDat_rt(ind_end:end,{'patm' 'Tair'}) = {NaN};

% Find where there are gaps in Gull met station air T data
for i = 1:height(metDat_rt.datetime_utc)
    nid = find(isnan(metDat_rt.Tair));
    id(1:length(nid),1) = nid;
end

% See if the PAR or the BP sensor do a better job of matching the met station data
figure(2),clf
plot(metDat_rt.datetime_utc,metDat_rt.Tair,'.','MarkerSize',10)
hold on
plot(parDat_rt.datetime_utc,parDat_rt.Tair,'.r','MarkerSize',4)
plot(bpDat_rt.datetime_utc,bpDat_rt.Tair,'.g','MarkerSize',4)
ylabel('Air Temperature (^oC)')
legend('Gull Met Station','PAR Dataset','Baro Pressure Dataset')
title('Assessing Air Temperature Datasets')

% Replace missing Gull T_air data with HOBO Baro Pressure T_air data
ind_nan = find(isnan(metDat_rt.Tair));
metDat_rt.Tair(ind_nan) = bpDat_rt.Tair(ind_nan);

figure(3),clf
plot(metDat_rt.datetime_utc,metDat_rt.Tair,'.','MarkerSize',4)
hold on
plot(metDat_rt.datetime_utc(ind_nan),metDat_rt.Tair(ind_nan),'og','MarkerSize',6,'LineWidth',1)
ylabel('Air Temperature (^oC)')
legend('Gull Met Station','Baro Pressure Dataset','location','southeast')
title('Gap-Filled Air Temperature Data')

%====Gap filling Gull met station atmos p data============================
figure(4),clf
plot(metDat_rt.datetime_utc,metDat_rt.patm,'.','MarkerSize',4)
hold on
plot(bpDat_rt.datetime_utc,bpDat_rt.patm,'.g','MarkerSize',4)
ylabel('p_{atm} (hPa)')
legend('Gull Met Station','Baro Pressure Dataset')
title('Assessing Atmospheric Pressure Datasets')

% Replace missing Gull p_atm data with HOBO Baro Pressure p_atm data
ind_nan = find(isnan(metDat_rt.patm));
metDat_rt.patm(ind_nan) = bpDat_rt.patm(ind_nan);

figure(5),clf
plot(metDat_rt.datetime_utc,metDat_rt.patm,'.','MarkerSize',4)
hold on
plot(metDat_rt.datetime_utc(ind_nan),metDat_rt.patm(ind_nan),'og','MarkerSize',6,'LineWidth',1)
ylabel('p_{atm} (hPa)')
legend('Gull Met Station','Baro Pressure Dataset')
title('Gap-Filled Atmospheric Pressure Data')

% Horizontally concatenate sonde and met data (already have common time vector)
dat = synchronize(dat,metDat_rt);
% dat = rmmissing(dat,'DataVariables',{'DO_conc','wspd'});

% Check how many days of data remain in synchronized data
dataPoints = daysact(dat.datetime_utc(1),dat.datetime_utc);
wholeDays = length(unique(floor(dataPoints)));

%% STEP 1
%====Calculate DO concentration at equilibrium (DO_sat)=================
DO_sat = O2sol(dat.salinity,dat.temperature);     % [umol kg-1]
% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(dat.salinity,dat.temperature,dat.p); % [kg m-3]
% Convert DO_sat units
DO_sat = DO_sat.*rho_sw/1000;       % [mmol m-3]
% Convert DO_conc units
DO_conc = dat.DO_conc*1000/1000;    % [mmol m-3]
% Calculate the percent oxygen saturation
DO_per_sat = dat.DO_conc./DO_sat*100;   % [%]

%% STEP 2
%====Calculate gas exchange at air-water interface=======================
% Caffrey/Beck (uses "reaeration coefficient", ka)

H = mean(dat.depth,'omitnan');      % Mean water depth [m] -- add a constant (e.g., 0.5 m) to account for height of sonde above bottom?
Tkair = dat.Tair + 273.15;          % Absolute air temperature [K]
R_specific = 287.0500676;           % Specific gas constant for dry air [J kg-1 K-1]
rho_a = dat.patm*100 ./ (R_specific*Tkair); % Density of dry air [kg m-3]
T = dat.temperature;                % Water temperature [degC]
S = dat.salinity;
u = dat.wspd;

% Calculate Vw
for i = 1:length(T)
    Uw(i,1) = swp('m',S(i),T(i));  % Dynamic viscosity [kg m-1 s-1]
end
Vw = Uw./rho_sw;    % Kinematic viscosity [m2 s-1]
kB = 1.3806503E-23; % Boltzmann's constant [m2 kg s-2 K-1]
R0 = 1.72E-10;      % Radius of O2 molecule [m]
Tk = T + 273.15;    % Absolute temperature [K]
ka = 1/H * 0.1706 .* (kB*Tk./(4*rho_sw.*Vw.^2.*R0)).^0.5 .* (rho_a./rho_sw).^0.5 .* u.^1.81;

% Air-water flux
D = ka.*(DO_sat - DO_conc); % [mmol m-3 h-1]

%% STEP 3
%====Determine day/night times=============================================
% Get day/night indices based on lat and lon and Hilary's indexDayNight function
lat = 39.08;
lon = -74.78;
time_in = dat.datetime_utc;
tol = 0;
UTCoffset = 0;  % Input datetime vector is in UTC

[dayind,nightind] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);

% Manually define DST start and end times
start1 = datenum('03/14/2021 02:00','mm/dd/yyyy HH:MM');
end1 = datenum('11/07/2021 02:00','mm/dd/yyyy HH:MM');
start2 = datenum('03/13/2022 02:00','mm/dd/yyyy HH:MM');
end2 = datenum('11/6/2022 02:00','mm/dd/yyyy HH:MM');
start3 = datenum('03/12/2023 02:00','mm/dd/yyyy HH:MM');
end3 = datenum('11/5/2023 02:00','mm/dd/yyyy HH:MM');

% For plotting, create datetimes for when EDT begins/ends each year during data record
end1_dt = datetime(end1,'ConvertFrom','datenum','TimeZone','UTC');
start2_dt = datetime(start2,'ConvertFrom','datenum','TimeZone','UTC');
end2_dt = datetime(end2,'ConvertFrom','datenum','TimeZone','UTC');
start3_dt = datetime(start3,'ConvertFrom','datenum','TimeZone','UTC');
end3_dt = datetime(end3,'ConvertFrom','datenum','TimeZone','UTC');

figure(3),clf
plot(dat.datetime_utc(dayind),dat.DO_conc(dayind),'.b','MarkerSize',12)
hold on
plot(dat.datetime_utc(nightind),dat.DO_conc(nightind),'.k','MarkerSize',12)
xline(end1_dt,'--','label','DST Ends')
xline(start2_dt,'--','label','DST Starts')
xline(end2_dt,'--','label','DST Ends')
xline(start3_dt,'--','label','DST Starts')
xline(end3_dt,'--','label','DST Ends')
ylabel('DO conc (\mumol/L)')
xlabel('UTC')

figure(6),clf
plot(dat.datetime_local(dayind),dat.DO_conc(dayind),'.b','MarkerSize',12)
hold on
plot(dat.datetime_local(nightind),dat.DO_conc(nightind),'.k','MarkerSize',12)
xline(end1_dt,'--','label','DST Ends')
xline(start2_dt,'--','label','DST Starts')
xline(end2_dt,'--','label','DST Ends')
xline(start3_dt,'--','label','DST Starts')
xline(end3_dt,'--','label','DST Ends')
ylabel('DO conc (\mumol/L)')
xlabel('Local')

% tiledlayout(2,1)
% ax1 = nexttile;
% plot(dat.datetime_local(dayind),dat.DO_conc(dayind),'.b','MarkerSize',12)
% hold on
% plot(dat.datetime_local(nightind),dat.DO_conc(nightind),'.k','MarkerSize',12)
% xline(end1_dt,'--','label','DST Ends')
% xline(start2_dt,'--','label','DST Starts')
% xline(end2_dt,'--','label','DST Ends')
% xline(start3_dt,'--','label','DST Starts')
% xline(end3_dt,'--','label','DST Ends')
% ylabel('DO conc (\mumol/L)')
% xlabel('Local')
% ax2 = nexttile;
% plot(dat.datetime_utc(dayind),dat.DO_conc(dayind),'.b','MarkerSize',12)
% hold on
% plot(dat.datetime_utc(nightind),dat.DO_conc(nightind),'.k','MarkerSize',12)
% xline(end1_dt,'--','label','DST Ends')
% xline(start2_dt,'--','label','DST Starts')
% xline(end2_dt,'--','label','DST Ends')
% xline(start3_dt,'--','label','DST Starts')
% xline(end3_dt,'--','label','DST Ends')
% ylabel('DO conc (\mumol/L)')
% xlabel('UTC')
% linkaxes([ax1 ax2],'x')

% Find the indices for when each day starts and stops
daystart = dayind(find(diff(dayind) > 1) + 1);
dayend = dayind(find(diff(dayind) > 1));

% Length of each day
% daylength = dat.datetime_local(dayend(2:end)) - dat.datetime_local(daystart(1:end-1));
daylength = dat.datetime_local(dayend(2:end)) - dat.datetime_local(daystart(1:end-1) - 1);
daylength = hours(daylength);

daystart_dt = dat.datetime_local(daystart(1:end-1));
dayend_dt = dat.datetime_local(dayend(2:end));

figure(7),clf;
plot(dat.datetime_local(daystart(1:end-1)),daylength,'.')
ylabel('Day length (h)')

% Plot day/night indices with PAR as a check
figure(8),clf
yyaxis left
plot(dat.datetime_local(dayind),dat.DO_conc(dayind),'.b','MarkerSize',6)
hold on
plot(dat.datetime_local(nightind),dat.DO_conc(nightind),'.k','MarkerSize',6)
set(gca,'ycolor','k')
ylabel('DO conc (\mumol/L)')
yyaxis right
plot(parDat_rt.datetime_local,parDat_rt.par,'.g','MarkerSize',6)
set(gca,'ycolor','g')
ylabel('PAR')
xlabel('Local')
xline(end1_dt,'--','label','DST Ends')
xline(start2_dt,'--','label','DST Starts')
xline(end2_dt,'--','label','DST Ends')
xline(start3_dt,'--','label','DST Starts')
xline(end3_dt,'--','label','DST Ends')
title(['Gull - ',sondename,' Sonde'])

%====Calculate rates=======================================================
dCdt = nan(length(DO_conc),1);
dCdt(2:end,1) = diff(DO_conc) ./ hours(diff(dat.datetime_utc));  % [mmol m-3 h-1]

% Hourly rates of respiration and net production
for i = 1:length(daylength)
    % During night hours, P = 0
    R_hourly(i,1) = mean(dCdt(dayend(i+1):daystart(i+1)) - D(dayend(i+1):daystart(i+1)),'omitnan'); % Hourly rate of respiration; [mmol m-3 h-1]
    % During day hours, P != 0
    P_hourly(i,1) = mean(dCdt(daystart(i):dayend(i+1)) - D(daystart(i):dayend(i+1)),'omitnan');     % Hourly rate of net production; [mmol m-3 h-1]
end

figure(6),clf
yyaxis left;plot(dat.datetime_utc,dCdt,'.');ylabel('dC/dt (mmol m^{-3} h^{-1})')
yyaxis right;plot(dat.datetime_utc,D,'.');ylabel('D (mmol m^{-3} h^{-1})')
title(['Mass Balance - Gull ',sondename,' Sonde'])

% Daily rates
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

% Convert volumetric rates to depth-integrated (areal) estimates
GPP = P_daily * H;      % [mmol O2 m-2 d-1]
ER = R_daily * H;       % [mmol O2 m-2 d-1]
NEM = GPP + ER;         % [mmol O2 m-2 d-1]

metab_tbl = table(daystart_dt,dayend_dt,daylength,R_hourly,P_hourly,R_daily,P_daily,GPP,ER,NEM);

% Plot metabolism results
figure(7),clf
plot(dat.datetime_utc(daystart(1:end-1)),GPP,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(dat.datetime_utc(daystart(1:end-1)),ER,'k.-','MarkerSize',12,'LineWidth',1)
plot(dat.datetime_utc(daystart(1:end-1)),NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title(['Gull - ',sondename,' Sonde'])

% Plot DO concentration and metabolism results 
figure(8),clf
tiledlayout(2,1)
ax1 = nexttile;
plot(dat.datetime_utc,dat.DO_conc,'.-','MarkerSize',6,'LineWidth',1)
ylabel('DO conc (mmol m^{-3})','FontSize',14)
title(['Gull - ',sondename,' Sonde'])
ax2 = nexttile;
plot(dat.datetime_utc(daystart(1:end-1)),GPP,'.-','MarkerSize',12,'LineWidth',1);
hold on
plot(dat.datetime_utc(daystart(1:end-1)),ER,'k.-','MarkerSize',12,'LineWidth',1)
plot(dat.datetime_utc(daystart(1:end-1)),NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
linkaxes([ax1 ax2],'x')

% Plot DO concentration and tidal level (c.f. Beck Fig. 6)
figure(9),clf
tiledlayout(2,1)
ax1 = nexttile;
plot(dat.datetime_utc,DO_conc,'.-','MarkerSize',6,'Linewidth',1)
ylabel('DO conc (mmol m^{-3})','FontSize',14)
title(['Gull - ',sondename,' Sonde'])
ax2 = nexttile;
plot(dat.datetime_utc,dat.depth,'.-','MarkerSize',6,'Linewidth',1)
ylabel('Water Depth (m)')
linkaxes([ax1 ax2],'x')

cd('G:\My Drive\Postdoc\Work\SMIIL\figures\diel-analysis-figures')

%%
%====Export data in proper format to run Beck's R code=====================
DO_mgL = dat.DO_conc / 1000 * 31.999;
varNames = ["DateTimeStamp","Temp","Sal","DO_obs","ATemp","BP","WSpd","Tide"];
dat_tbl = table(dat.datetime_utc,dat.temperature,dat.salinity,DO_mgL,dat.Tair,dat.patm,dat.wspd,dat.depth,'VariableNames',varNames);
dat_tbl = rmmissing(dat_tbl);

switch answer
    case 'Sonde 1'
        writetable(dat_tbl,'G:\My Drive\Postdoc\Work\SMIIL\diel-method\owp-data\gull-bc.csv')
    case 'Sonde 2'
        writetable(dat_tbl,'G:\My Drive\Postdoc\Work\SMIIL\diel-method\owp-data\gull-erdc.csv')
end

%%
%====Import metabolism results from running Beck's R code==================

% See help("ecometab") in R for variable names, units, etc.
varNames = ["datetime_local","Pg","Rt","NEM","Pg_vol","Rt_vol"];
varUnits = ["","mmol O2 m-2 d-1","mmol O2 m-2 d-1","mmol O2 m-2 d-1","mmol O2 m-3 d-1","mmol O2 m-3 d-1"];

cd([rootpath,'diel-method\owp-results\gull-erdc\observed'])
metab_obs = readtable('metab_obs.csv');
metab_obs.Properties.VariableNames = varNames;
metab_obs.Properties.VariableUnits = varUnits;

cd([rootpath,'diel-method\owp-results\gull-erdc\sapelo'])
metab_dtd = readtable('metab_dtd.csv');
metab_dtd.Properties.VariableNames = varNames;
metab_dtd.Properties.VariableUnits = varUnits;

cd('G:\My Drive\Postdoc\Work\SMIIL\figures\diel-analysis-figures')

% Plot metabolism results from Beck's R code
figure(10),clf
plot(metab_obs.datetime_local,metab_obs.Pg,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(metab_obs.datetime_local,metab_obs.Rt,'k.-','MarkerSize',12,'LineWidth',1)
plot(metab_obs.datetime_local,metab_obs.NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
title(['Beck Metabolism Results for Gull ',sondename,' Sonde - Not Detided'])
legend({'GPP','ER','NEM'},'FontSize',14)

figure(11),clf
plot(metab_dtd.datetime_local,metab_dtd.Pg,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(metab_dtd.datetime_local,metab_dtd.Rt,'k.-','MarkerSize',12,'LineWidth',1)
plot(metab_dtd.datetime_local,metab_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
title(['Beck Metabolism Results for Gull ',sondename,' Sonde - Detided'])
legend({'GPP','ER','NEM'},'FontSize',14)

% Plot in same figure with linked axes
figure(12),clf
tiledlayout(2,1)
ax1 = nexttile;
plot(metab_obs.datetime_local,metab_obs.Pg,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(metab_obs.datetime_local,metab_obs.Rt,'k.-','MarkerSize',12,'LineWidth',1)
plot(metab_obs.datetime_local,metab_obs.NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
title(['Beck Metabolism Results for Gull ',sondename,' Sonde - Not Detided'])
ax2 = nexttile;
plot(metab_dtd.datetime_local,metab_dtd.Pg,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(metab_dtd.datetime_local,metab_dtd.Rt,'k.-','MarkerSize',12,'LineWidth',1)
plot(metab_dtd.datetime_local,metab_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
title(['Beck Metabolism Results for Gull ',sondename,' Sonde - Detided'])
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
linkaxes([ax1 ax2],'x')

%% Compare observed vs. detided results by month

window = 'Elkhorn';

% Grouping by 12 months
cd([rootpath,'diel-method\owp-results\gull-bc\observed'])
month_avg_obs = readtable('meteval_mos_obs.csv');
cd([rootpath,'diel-method\owp-results\gull-bc\',window])
month_avg_dtd = readtable('meteval_mos_dtd.csv');

cd([rootpath,'figures\diel-analysis-figures'])
% Plot means with errorbars (SD) and bars (anomalous values)
month = month_avg_obs.month';
anomPg_all = [month_avg_obs.anomPg'; month_avg_dtd.anomPg'];
anomRt_all = [month_avg_obs.anomRt'; month_avg_dtd.anomRt'];

figure(13),clf
errorbar(month_avg_obs.month,month_avg_obs.meanPg,month_avg_obs.sdPg,'LineWidth',2)
hold on
errorbar(month_avg_dtd.month,month_avg_dtd.meanPg,month_avg_dtd.sdPg,'LineWidth',2)
b = bar(month,anomPg_all);
b(1,1).FaceColor = [0 0.4470 0.7410];
b(1,2).FaceColor = [0.8500 0.3250 0.0980];
legend('Observed',['Detided (',window,' Window Widths)'],'Anomalous (Observed)','Anomalous (Detided)','location','best')
xlabel('Month')
ylabel('P_g (mmol O_2 m^{-2} d^{-1})','FontSize',14)
title(['Beck Metabolism Results for Gull ',sondename,' Sonde'])

figure(14),clf
errorbar(month_avg_obs.month,month_avg_obs.meanRt,month_avg_obs.sdPg,'LineWidth',2)
hold on
errorbar(month_avg_dtd.month,month_avg_dtd.meanRt,month_avg_dtd.sdPg,'LineWidth',2)
b = bar(month,anomRt_all);
b(1,1).FaceColor = [0 0.4470 0.7410];
b(1,2).FaceColor = [0.8500 0.3250 0.0980];
legend('Observed',['Detided (',window,' Window Widths)'],'Anomalous (Observed)','Anomalous (Detided)','location','best')
xlabel('Month')
ylabel('R_t (mmol O_2 m^{-2} d^{-1})','FontSize',14)
title(['Beck Metabolism Results for Gull ',sondename,' Sonde'])

% Plot anomalous values only
anomPg_all = [month_avg_obs.anomPg'; month_avg_dtd.anomPg'];
figure(15),clf
bar(month,anomPg_all)
xlabel('Month')
ylabel('Anomalous P_g (mmol O_2 m^{-2} d^{-1})')
legend('Observed',['Detided (',window,' Window Widths)'])
title(['Beck Metabolism Results for Gull ',sondename,' Sonde'])

anomRt_all = [month_avg_obs.anomRt'; month_avg_dtd.anomRt'];
figure(16),clf
bar(month,anomRt_all)
xlabel('Month')
ylabel('Anomalous R_t (mmol O_2 m^{-2} d^{-1})')
legend('Observed',['Detided (',window,' Window Widths)'])
title(['Beck Metabolism Results for Gull ',sondename,' Sonde'])

% Grouping by month-year
Pg_monthmeans_obs = groupsummary(metab_obs,'datetime_local','month','mean');
Pg_monthmeans_dtd = groupsummary(metab_dtd,'datetime_local','month','mean');
Pg_monthmeans_all = [Pg_monthmeans_obs.mean_Pg'; Pg_monthmeans_dtd.mean_Pg'];
month_yr = Pg_monthmeans_obs.month_datetime_local';

figure(17),clf
bar(month_yr,Pg_monthmeans_all)
ylabel('P_g (mmol O_2 m^{-2} d^{-1})','FontSize',14)
legend('Observed',['Detided (',window,' Window Widths)'],'location','northwest')
title(['Beck Metabolism Results for Gull ',sondename,' Sonde'])

Rt_monthmeans_obs = groupsummary(metab_obs,'datetime_local','month','mean');
Rt_monthmeans_dtd = groupsummary(metab_dtd,'datetime_local','month','mean');
Rt_monthmeans_all = [Rt_monthmeans_obs.mean_Rt'; Rt_monthmeans_dtd.mean_Rt'];
month_yr = Rt_monthmeans_obs.month_datetime_local';

figure(18),clf
bar(month_yr,Rt_monthmeans_all)
ylabel('R_t (mmol O_2 m^{-2} d^{-1})','FontSize',14)
legend('Observed',['Detided (',window,' Window Widths)'],'location','southwest')
title(['Beck Metabolism Results for Gull ',sondename,' Sonde'])

%% Take the mean of the metabolism results from the BC and ERDC datasets

window = 'Elkhorn';

varNames = ["Pg","Rt","NEM","Pg_vol","Rt_vol"];
varUnits = ["mmol O2 m-2 d-1","mmol O2 m-2 d-1","mmol O2 m-2 d-1","mmol O2 m-3 d-1","mmol O2 m-3 d-1"];

% Load BC data
cd([rootpath,'diel-method\owp-results\gull-bc\observed'])
metab_obs = readtable('metab_obs.csv');
metab_obs = table2timetable(metab_obs);
metab_obs.Properties.VariableNames = varNames;
metab_obs.Properties.VariableUnits = varUnits;

cd([rootpath,'diel-method\owp-results\gull-bc\',window])
metab_dtd = readtable('metab_dtd.csv');
metab_dtd = table2timetable(metab_dtd);
metab_dtd.Properties.VariableNames = varNames;
metab_dtd.Properties.VariableUnits = varUnits;

bc_metab = struct('obs',metab_obs,'dtd',metab_dtd);

% Load ERDC data
cd([rootpath,'diel-method\owp-results\gull-erdc\observed'])
metab_obs = readtable('metab_obs.csv');
metab_obs = table2timetable(metab_obs);
metab_obs.Properties.VariableNames = varNames;
metab_obs.Properties.VariableUnits = varUnits;

cd([rootpath,'diel-method\owp-results\gull-erdc\',window])
metab_dtd = readtable('metab_dtd.csv');
metab_dtd = table2timetable(metab_dtd);
metab_dtd.Properties.VariableNames = varNames;
metab_dtd.Properties.VariableUnits = varUnits;

erdc_metab = struct('obs',metab_obs,'dtd',metab_dtd);
%%
% Average the observed results and plot
cd([rootpath,'figures\diel-analysis-figures'])

all_metab_obs = synchronize(bc_metab.obs,erdc_metab.obs);

Pg_obs = mean(all_metab_obs{:,["Pg_1","Pg_2"]},2,"includemissing");
Rt_obs= mean(all_metab_obs{:,["Rt_1","Rt_2"]},2,"includemissing");
NEM_obs = mean(all_metab_obs{:,["NEM_1","NEM_2"]},2,"includemissing");

varNames = ["Pg","Rt","NEM"];
varUnits = ["mmol O2 m-2 d-1","mmol O2 m-2 d-1","mmol O2 m-2 d-1"];

mean_metab_obs = timetable(all_metab_obs.Date,Pg_obs,Rt_obs,NEM_obs);
mean_metab_obs.Properties.VariableNames = varNames;
mean_metab_obs.Properties.VariableUnits = varUnits;

figure(19),clf
plot(mean_metab_obs.Time,mean_metab_obs.Pg,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(mean_metab_obs.Time,mean_metab_obs.Rt,'k.-','MarkerSize',12,'LineWidth',1)
plot(mean_metab_obs.Time,mean_metab_obs.NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
title('Beck Metabolism Results for Gull - Mean of Both Sondes - Not Detided')
legend({'GPP','ER','NEM'},'FontSize',14)

% figure(20),clf
% plot(bc_metab.obs.Date,bc_metab.obs.Pg,'.-','MarkerSize',12,'LineWidth',1)
% hold on
% plot(erdc_metab.obs.Date,erdc_metab.obs.Pg,'.-','MarkerSize',12,'LineWidth',1)
% plot(mean_metab_obs.Time,mean_metab_obs.Pg,'.-','MarkerSize',12,'LineWidth',1)
% legend('BC','ERDC','Mean')
%%
% Average the detided results and plot
cd([rootpath,'figures\diel-analysis-figures'])

all_metab_dtd = synchronize(bc_metab.dtd,erdc_metab.dtd);

Pg_dtd = mean(all_metab_dtd{:,["Pg_1","Pg_2"]},2,"includemissing");
Rt_dtd= mean(all_metab_dtd{:,["Rt_1","Rt_2"]},2,"includemissing");
NEM_dtd = mean(all_metab_dtd{:,["NEM_1","NEM_2"]},2,"includemissing");

varNames = ["Pg","Rt","NEM"];
varUnits = ["mmol O2 m-2 d-1","mmol O2 m-2 d-1","mmol O2 m-2 d-1"];

mean_metab_dtd = timetable(all_metab_dtd.Date,Pg_dtd,Rt_dtd,NEM_dtd);
mean_metab_dtd.Properties.VariableNames = varNames;
mean_metab_dtd.Properties.VariableUnits = varUnits;

figure(21),clf
plot(mean_metab_dtd.Time,mean_metab_dtd.Pg,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(mean_metab_dtd.Time,mean_metab_dtd.Rt,'k.-','MarkerSize',12,'LineWidth',1)
plot(mean_metab_dtd.Time,mean_metab_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
title('Beck Metabolism Results for Gull - Mean of Both Sondes - Detided')
legend({'GPP','ER','NEM'},'FontSize',14)

% figure(22),clf
% plot(bc_metab.dtd.Date,bc_metab.dtd.Pg,'.-','MarkerSize',12,'LineWidth',1)
% hold on
% plot(erdc_metab.dtd.Date,erdc_metab.dtd.Pg,'.-','MarkerSize',12,'LineWidth',1)
% plot(mean_metab_dtd.Time,mean_metab_dtd.Pg,'.-','MarkerSize',12,'LineWidth',1)
% legend('BC','ERDC','Mean')