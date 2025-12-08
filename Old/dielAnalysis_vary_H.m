%%%%%%%%%%%%%%%%%/%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dielAnalysis_vary_H.m
% This script uses the diel oxygen method to calculate GPP, ER, and NEM
% following Caffrey et al. (2014) and Beck et al. (2015).
% 
% Computes summary statistics for both the observed and detided data.
%
% Output diel analysis results from the detided data (diel_dtd) include
% anomalous values, but have endpoints around big gaps of data removed.
%
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 9/4/2024
% Last updated: 9/6/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

cd([rootpath,'\diel-method\owp-data\final-qc'])

load([site,'_obs.mat'])   % Load the table

% Load R results
cd([rootpath,'diel-method\R-results\final-qc\',site])
wtreg_res = readtable('wtreg_res.csv');
wtreg_res.Properties.DimensionNames{1} = 'datetime_utc';
wtreg_res.DateTimeStamp.TimeZone = "UTC";
wtreg_res = table2timetable(wtreg_res);

% Retime weighted regression data to same datetimes as sonde/physical data
wtreg_res.solar_period = [];
newTimes = dat.datetime_utc(1):minutes(10):dat.datetime_utc(end);
wtreg_res_rt = retime(wtreg_res,newTimes,'mean');

metab_obs = table2timetable(readtable('metab_obs.csv'));
metab_dtd = table2timetable(readtable('metab_dtd.csv'));

%====Define input variables================================================
dt_utc = wtreg_res_rt.DateTimeStamp;
dt_local = dt_utc;
dt_local.TimeZone = 'America/New_York';

% Water
S = wtreg_res_rt.Sal;
T = wtreg_res_rt.Temp;
% p = dat.p;
p = wtreg_res_rt.Tide;  % Pressure in dbar and depth in meters are approx. equal
d = wtreg_res_rt.Tide;
% Mean water column depth, H = d + D [m]
% See Collab Lab Notebook - Table 2 for manual measurements for D for each site
switch site
    case 'Gull'
        H = mean(dat.depth,'omitnan') + 0.42;
    case 'North'
        H = mean(dat.depth,'omitnan') + 0.47;
    case 'South'
        H = mean(dat.depth,'omitnan') + 0.80;
end
% Air
Tair = wtreg_res_rt.ATemp;  % [deg C]
patm = wtreg_res_rt.BP;     % [hPa]
u = wtreg_res_rt.WSpd;      % [m/s]

% DO concentration conversions
DO_obs = wtreg_res_rt.DO_obs*1000/32;   % Observed DO concentration [mmol m-3]
DO_nrm = wtreg_res_rt.DO_nrm*1000/32;   % Detided DO concentration [mmol m-3]

% Plot results of R detiding (c.f. Beck Fig. 6)
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
t = tiledlayout(2,1,'TileSpacing','compact');
ax1 = nexttile;
plot(dt_local,DO_obs,'.-','MarkerSize',6,'Linewidth',1,'DisplayName','Observed')
hold on
plot(dt_local,DO_nrm,'.-','MarkerSize',6,'Linewidth',1,'DisplayName','Detided')
legend('show','location','best')
ylabel('DO conc. (mmol m^{-3})','FontSize',14)
title(site)
ax2 = nexttile;
plot(dt_local,d,'.-','MarkerSize',6,'Linewidth',1)
xlabel('Local Time')
ylabel('Water depth (m)')
linkaxes([ax1 ax2],'x')

%==========================================================================
%   Conduct diel analysis with observed DO data
%==========================================================================
DO_conc = DO_obs;

%====STEP 1: Calculate DO concentration at equilibrium (DO_sat)============
DO_sat = O2sol(S,T);                % [umol kg-1]

% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(S,T,p); % [kg m-3]

% Convert DO_sat units
DO_sat = DO_sat.*rho_sw/1000;       % [mmol m-3]

% Calculate the percent oxygen saturation
DO_per_sat = DO_conc./DO_sat*100;   % [%]

%====STEP 2: Calculate gas exchange at air-water interface=================
% Caffrey/Beck (uses "reaeration coefficient", ka)
Tkair = Tair + 273.15;                  % Absolute air temperature [K]
R_specific = 287.0500676;               % Specific gas constant for dry air [J kg-1 K-1]
rho_a = patm*100 ./ (R_specific*Tkair); % Density of dry air [kg m-3]

% Calculate Vw (kinematic viscosity of sw)
for i = 1:length(T)
    Uw(i,1) = swp('m',S(i),T(i));  % Dynamic viscosity of sw [kg m-1 s-1]
end
Vw = Uw./rho_sw;    % Kinematic viscosity [m2 s-1]
kB = 1.3806503E-23; % Boltzmann's constant [m2 kg s-2 K-1]
R0 = 1.72E-10;      % Radius of O2 molecule [m]
Tk = T + 273.15;    % Absolute temperature [K]
ka = 1/H * 0.1706 .* (kB*Tk./(4*rho_sw.*Vw.^2.*R0)).^0.5 .* (rho_a./rho_sw).^0.5 .* u.^1.81;

% Air-water flux
D = ka.*(DO_sat - DO_conc); % [mmol m-3 h-1]

%====STEP 3: Determine day/night times=====================================
% Get day/night indices based on lat and lon and Hilary's indexDayNight function
lat = 39.08;
lon = -74.78;
time_in = dt_utc;
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

fig2 = figure(2);clf
fig2.WindowState = 'maximized';
% plot(dt_utc(dayind),DO_conc(dayind),'.b','MarkerSize',12)
plot(dt_local(dayind),DO_conc(dayind),'.b','MarkerSize',12)
hold on
% plot(dt_utc(nightind),DO_conc(nightind),'.k','MarkerSize',12)
plot(dt_local(nightind),DO_conc(nightind),'.k','MarkerSize',12)
xline(end1_dt,'--','label','DST Ends')
xline(start2_dt,'--','label','DST Starts')
xline(end2_dt,'--','label','DST Ends')
xline(start3_dt,'--','label','DST Starts')
xline(end3_dt,'--','label','DST Ends')
ylabel('DO conc (mmol m^{-3})')
% xlabel('UTC')
xlabel('Local Time')
title([site,': Observed DO Concentration'])

% Find the indices for when each day starts and stops
daystart = dayind(find(diff(dayind) > 1) + 1);
dayend = dayind(find(diff(dayind) > 1));

% Length of each day
daylength = dt_utc(dayend(2:end)) - dt_utc(daystart(1:end-1) - 1);
daylength = hours(daylength);

daystart_dt = dt_utc(daystart(1:end-1));
dayend_dt = dt_utc(dayend(2:end));

fig3 = figure(3);clf
fig3.WindowState = 'maximized';
plot(dt_local(daystart(1:end-1)),daylength,'.')
% xlabel('Local')
ylabel('Day length (h)')
title([site,': Observed DO Concentration'])

%====STEP 4: Calculate rates===============================================
dCdt = nan(length(DO_conc),1);
dCdt(2:end,1) = diff(DO_conc) ./ hours(diff(dt_utc));  % [mmol m-3 h-1]

% Mean hourly rates
for i = 1:length(daylength)
    % During night hours, P = 0
    R_hourly(i,1) = mean(dCdt(dayend(i+1):daystart(i+1)) - D(dayend(i+1):daystart(i+1)),'omitnan'); % Mean hourly rate of nighttime respiration; [mmol m-3 h-1]
    % During day hours, P != 0
    P_hourly(i,1) = mean(dCdt(daystart(i):dayend(i+1)) - D(daystart(i):dayend(i+1)),'omitnan');     % Mean hourly rate of apparent/net production; [mmol m-3 h-1]
end

fig4 = figure(4);clf
fig4.WindowState = 'maximized';
yyaxis left;plot(dt_utc,dCdt,'.','markersize',6)
ylabel('dC/dt (mmol m^{-3} h^{-1})')
yyaxis right;plot(dt_utc,D,'.','markersize',6)
ylabel('D (mmol m^{-3} h^{-1})')
xlabel('UTC')
title([site,': Observed DO Concentration'])

% Daily rates
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

% Convert volumetric rates to depth-integrated (areal) estimates
GPP = P_daily * H;      % [mmol O2 m-2 d-1]
ER = R_daily * H;       % [mmol O2 m-2 d-1]
NEM = GPP + ER;         % [mmol O2 m-2 d-1]

diel_obs = table(daystart_dt,dayend_dt,daylength,R_hourly,P_hourly,R_daily,P_daily,GPP,ER,NEM);
diel_obs = table2timetable(diel_obs);

%==========================================================================
%   Conduct diel analysis with detided DO data
%==========================================================================
DO_conc = DO_nrm;

%====STEP 1: Calculate DO concentration at equilibrium (DO_sat)============
DO_sat = O2sol(S,T);              % [umol kg-1]

% Use Gibbs Seawater toolbox to calculate seawater density
rho_sw = gsw_rho(S,T,p);          % [kg m-3]

% Convert DO_sat units
DO_sat = DO_sat.*rho_sw/1000;     % [mmol m-3]

% Calculate the percent oxygen saturation
DO_per_sat = DO_conc./DO_sat*100; % [%]

%====STEP 2: Calculate gas exchange at air-water interface=================
% Caffrey/Beck (uses "reaeration coefficient", ka)
Tkair = Tair + 273.15;                  % Absolute air temperature [K]
R_specific = 287.0500676;               % Specific gas constant for dry air [J kg-1 K-1]
rho_a = patm*100 ./ (R_specific*Tkair); % Density of dry air [kg m-3]

% Calculate Vw (kinematic viscosity of sw)
for i = 1:length(T)
    Uw(i,1) = swp('m',S(i),T(i));  % Dynamic viscosity of sw [kg m-1 s-1]
end
Vw = Uw./rho_sw;    % Kinematic viscosity [m2 s-1]
kB = 1.3806503E-23; % Boltzmann's constant [m2 kg s-2 K-1]
R0 = 1.72E-10;      % Radius of O2 molecule [m]
Tk = T + 273.15;    % Absolute temperature [K]
ka = 1/H * 0.1706 .* (kB*Tk./(4*rho_sw.*Vw.^2.*R0)).^0.5 .* (rho_a./rho_sw).^0.5 .* u.^1.81;

% Air-water flux
D = ka.*(DO_sat - DO_conc); % [mmol m-3 h-1]

%====STEP 3: Determine day/night times=====================================
% Get day/night indices based on lat and lon and Hilary's indexDayNight function
lat = 39.08;
lon = -74.78;
time_in = dt_utc;
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

fig5 = figure(5);clf
fig5.WindowState = 'maximized';
plot(dt_local(dayind),DO_conc(dayind),'.b','MarkerSize',12)
hold on
plot(dt_local(nightind),DO_conc(nightind),'.k','MarkerSize',12)
xline(end1_dt,'--','label','DST Ends')
xline(start2_dt,'--','label','DST Starts')
xline(end2_dt,'--','label','DST Ends')
xline(start3_dt,'--','label','DST Starts')
xline(end3_dt,'--','label','DST Ends')
ylabel('DO conc (mmol m^{-3})')
xlabel('Local Time')
title([site,': Detided DO Concentration'])

% Find the indices for when each day starts and stops
daystart = dayind(find(diff(dayind) > 1) + 1);
dayend = dayind(find(diff(dayind) > 1));

% Find corresponding datetimes
daystart_dt = dt_utc(daystart(1:end-1));
dayend_dt = dt_utc(dayend(2:end));

% Length of each day
daylength = dt_utc(dayend(2:end)) - dt_utc(daystart(1:end-1) - 1);
daylength = hours(daylength);

fig6 = figure(6);clf
fig6.WindowState = 'maximized';
plot(dt_local(daystart(1:end-1)),daylength,'.')
% xlabel('UTC')
ylabel('Day length (h)')
title([site,': Detided DO Concentration'])

%====STEP 4: Calculate rates===============================================
dCdt = nan(length(DO_conc),1);
dCdt(2:end,1) = diff(DO_conc) ./ hours(diff(dt_utc));  % [mmol m-3 h-1]

% Hourly rates of respiration and net production
for i = 1:length(daylength)
    % During night hours, P = 0
    R_hourly(i,1) = mean(dCdt(dayend(i+1):daystart(i+1)) - D(dayend(i+1):daystart(i+1)),'omitnan'); % Hourly rate of respiration; [mmol m-3 h-1]
    % During day hours, P != 0
    P_hourly(i,1) = mean(dCdt(daystart(i):dayend(i+1)) - D(daystart(i):dayend(i+1)),'omitnan');     % Hourly rate of net production; [mmol m-3 h-1]
end

fig7 = figure(7);clf
fig7.WindowState = 'maximized';
yyaxis left;plot(dt_utc,dCdt,'.','markersize',6)
ylabel('dC/dt (mmol m^{-3} h^{-1})')
yyaxis right;plot(dt_utc,D,'.','markersize',6)
xlabel('UTC')
ylabel('D (mmol m^{-3} h^{-1})')
title([site,': Detided DO Concentration'])

% Daily rates
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

% Convert volumetric rates to depth-integrated (areal) estimates
GPP = P_daily * H;      % [mmol O2 m-2 d-1]
ER = R_daily * H;       % [mmol O2 m-2 d-1]
NEM = GPP + ER;         % [mmol O2 m-2 d-1]

diel_dtd = table(daystart_dt,dayend_dt,daylength,R_hourly,P_hourly,R_daily,P_daily,GPP,ER,NEM);
diel_dtd = table2timetable(diel_dtd);

%==========================================================================
%   Compare using time-varying H with constant mean H -- this part is
%   different from dielAnalysis.m
%==========================================================================
daily_depth = groupsummary(dat,"datetime_utc","day","mean","depth");
switch site
    case 'Gull'
        H = daily_depth.mean_depth + 0.42;
    case 'North'
        H = daily_depth.mean_depth + 0.47;
    case 'South'
        H = daily_depth.mean_depth + 0.80;
end

% Convert volumetric rates to depth-integrated (areal) estimates using daily H
GPP_vary = P_daily .* H(2:end-1);
ER_vary = R_daily .* H(2:end-1);
NEM_vary = GPP_vary + ER_vary;
diel_dtd_vary = table(daystart_dt,dayend_dt,daylength,R_hourly,P_hourly,R_daily,P_daily,GPP_vary,ER_vary,NEM_vary);

figure,clf
plot(diel_dtd.daystart_dt,diel_dtd.GPP,'DisplayName','Constant H')
hold on
plot(diel_dtd_vary.daystart_dt,diel_dtd_vary.GPP_vary,'DisplayName','Time-varying H')
ylabel('GPP (mmol O_2 m^{-2} d^{-1})','FontSize',14)
legend('show')
title(site)

figure,clf
plot(diel_dtd.daystart_dt,diel_dtd.ER,'DisplayName','Constant H')
hold on
plot(diel_dtd_vary.daystart_dt,diel_dtd_vary.ER_vary,'DisplayName','Time-varying H')
ylabel('ER (mmol O_2 m^{-2} d^{-1})','FontSize',14)
legend('show')
title(site)

figure,clf
plot(diel_dtd.daystart_dt,diel_dtd.NEM,'DisplayName','Constant H')
hold on
plot(diel_dtd_vary.daystart_dt,diel_dtd_vary.NEM_vary,'DisplayName','Time-varying H')
ylabel('NEM (mmol O_2 m^{-2} d^{-1})','FontSize',14)
legend('show')
title(site)

%==========================================================================
%   Make summary statistics tables
%==========================================================================
% Remove rows with missing data
diel_obs = rmmissing(diel_obs);
diel_dtd = rmmissing(diel_dtd);

% Observed data
anomGPP = length(find(diel_obs.GPP < 0)) / length(diel_obs.GPP) * 100;
anomER = length(find(diel_obs.ER > 0)) / length(diel_obs.ER) * 100;
meanGPP = mean(diel_obs.GPP);
sdGPP = std(diel_obs.GPP);
meanER = mean(diel_obs.ER);
sdER = std(diel_obs.ER);
summStats.obs = table(meanGPP,sdGPP,anomGPP,meanER,sdER,anomER);

% Detided data
anomGPP = length(find(diel_dtd.GPP < 0)) / length(diel_dtd.GPP) * 100;
anomER = length(find(diel_dtd.ER > 0)) / length(diel_dtd.ER) * 100;
meanGPP = mean(diel_dtd.GPP);
sdGPP = std(diel_dtd.GPP);
meanER = mean(diel_dtd.ER);
sdER = std(diel_dtd.ER);
summStats.dtd = table(meanGPP,sdGPP,anomGPP,meanER,sdER,anomER);

%==========================================================================
%   Plot results from R ecometab and MATLAB diel analysis
%==========================================================================

% Observed
fig8 = figure(8);clf
fig8.WindowState = 'maximized';
plot(metab_obs.Date,metab_obs.Pg,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(metab_obs.Date,metab_obs.Rt,'k.-','MarkerSize',12,'LineWidth',1)
plot(metab_obs.Date,metab_obs.NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
title([site,': R "ecometab" Results Using Observed Data'])
legend({'GPP','ER','NEM'},'FontSize',14)
% ylim([-1000 800])

fig9 = figure(9);clf
fig9.WindowState = 'maximized';
plot(diel_obs.daystart_dt,diel_obs.GPP,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(diel_obs.daystart_dt,diel_obs.ER,'k.-','MarkerSize',12,'LineWidth',1)
plot(diel_obs.daystart_dt,diel_obs.NEM,'r.-','MarkerSize',12,'LineWidth',1)
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,': MATLAB Results Using Observed Data'])
% ylim([-1000 800])

% Detided
fig10 = figure(10);clf
fig10.WindowState = 'maximized';
plot(metab_dtd.Date,metab_dtd.Pg,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(metab_dtd.Date,metab_dtd.Rt,'k.-','MarkerSize',12,'LineWidth',1)
plot(metab_dtd.Date,metab_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1)
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,': R "ecometab" Results Using Detided Data'])
% ylim([-500 500])

fig11 = figure(11);clf
fig11.WindowState = 'maximized';
plot(diel_dtd.daystart_dt,diel_dtd.GPP,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(diel_dtd.daystart_dt,diel_dtd.ER,'k.-','MarkerSize',12,'LineWidth',1)
plot(diel_dtd.daystart_dt,diel_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1)
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,': MATLAB Results Using Detided Data'])
% ylim([-500 500])

% Delete endpoint values of GPP, ER, and NEM around big gaps in data 
% (happens for North and South -- gaps accompanied by large NEM absolute values)
threshold = duration(days(7));
gap = diff(diel_dtd.daystart_dt);
idx = find(gap > threshold);
idx_delete = [idx;idx+1];

diel_dtd(idx_delete,:) = [];

fig11 = figure(11);clf
fig11.WindowState = 'maximized';
plot(diel_dtd.daystart_dt,diel_dtd.GPP,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(diel_dtd.daystart_dt,diel_dtd.ER,'k.-','MarkerSize',12,'LineWidth',1)
plot(diel_dtd.daystart_dt,diel_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1)
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,': MATLAB Results Using Detided Data'])
% ylim([-500 500])

%==========================================================================
%   Option to save the results and figures
%==========================================================================

%====Save the data=========================================================
option = questdlg('Save diel analysis results?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'diel-method\matlab-results\final-qc\',site])
        save('diel_res.mat','diel_obs','diel_dtd','summStats')
        disp('Files saved!')
    case 'No'
        disp('Files not saved.')
end

%====Save the plots========================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\diel-analysis\R-results\final-qc\',site])
        saveas(fig1,'detiding_results.fig')
        saveas(fig1,'detiding_results.png')

        % Save plots of R results
        cd([rootpath,'figures\diel-analysis\R-results\final-qc\',site])
        saveas(fig8,'ecometab_obs.png')
        saveas(fig8,'ecometab_obs.fig')
        saveas(fig10,'ecometab_dtd.png')
        saveas(fig10,'ecometab_dtd.fig')

        % Save plots of MATLAB results
        cd([rootpath,'figures\diel-analysis\matlab-results\final-qc\',site,'\sanity-checks'])
        saveas(fig2,'day-night_obs.fig')
        saveas(fig2,'day-night_obs.png')
        saveas(fig3,'daylength_obs.fig')
        saveas(fig3,'daylength_obs.png')
        saveas(fig4,'mass-balance_obs.fig')
        saveas(fig4,'mass-balance_obs.png')
        saveas(fig5,'day-night_dtd.fig')
        saveas(fig5,'day-night_dtd.png')
        saveas(fig6,'daylength_dtd.fig')
        saveas(fig6,'daylength_dtd.png')
        saveas(fig7,'mass-balance_dtd.fig')
        saveas(fig7,'mass-balance_dtd.png')

        cd([rootpath,'figures\diel-analysis\matlab-results\final-qc\',site])
        saveas(fig9,'dielAnalysis_obs.png')
        saveas(fig9,'dielAnalysis_obs.fig')
        saveas(fig11,'dielAnalysis_dtd.png')
        saveas(fig11,'dielAnalysis_dtd.fig')

        disp('Plots saved!')

    case 'No'
        disp('Plots not saved.')
end