clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
cd([rootpath,'diel-method\pond-data'])
pondNum = '1';
filename = ['pond',num2str(pondNum),'_obs.csv'];
opts = detectImportOptions(filename);
opts = setvaropts(opts,'DateTimeStamp','DatetimeFormat','MM/dd/yy HH:mm:ss');
dat_in = readtable(filename,opts);

%====Export data in proper format to run Beck's R code=====================
% Remove rows with any missing data
% indNaN = find(isnan(dat_in.Tide));
% indNaN = find(isnan(dat(:,2:end)));
% dat(indNaN,:) = [];

dat = rmmissing(dat_in);

% Remove duplicate rows
[dat2,I,J] = unique(dat(:,2:end),'rows');
indDups = setdiff(1:size(dat,1),I);
dat(indDups,:) = []; 

% Export .csv file to run WtRegDO in R
% writetable(dat,['pond',num2str(pondNum),'.csv'])

%% ====Run Beck's R code=====================================================
% Go to R and run WtRegDO

%====Import R results======================================================
pondNum = '5';

cd([rootpath,'diel-method\R-results\ponds\pond',num2str(pondNum)])
output = readtable('wtreg_res.csv');

output.DateTimeStamp.TimeZone = "UTC";
dt_local = output.DateTimeStamp;
dt_local.TimeZone = "America/New_York";
DO_obs = output.DO_obs;
DO_nrm = output.DO_nrm;
d = output.Tide;

% Plot observed & detided DO concentration and tidal level (c.f. Beck Fig. 6)
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
t = tiledlayout(2,1,'TileSpacing','compact');
ax1 = nexttile;
plot(dt_local,DO_obs,'.-','MarkerSize',6,'Linewidth',1)
hold on
plot(dt_local,DO_nrm,'.-','MarkerSize',6,'Linewidth',1)
legend('Observed','Detided')
ylabel('DO conc. (mmol m^{-3})','FontSize',14)
title(['Pond ',num2str(pondNum)])
ax2 = nexttile;
plot(dt_local,d,'.-','MarkerSize',6,'Linewidth',1)
xlabel('Local Time')
ylabel('Water depth (m)')
linkaxes([ax1 ax2],'x')

filename2 = 'metab_obs.csv';
varNames = {'Pg','Rt','NEM','Pg_vol','Rt_vol'};
varTypes = {'double','double','double','double','double'};
opts = setvartype(opt,varNames,varTypes);
metab_obs = readtable('metab_obs.csv',opts);
metab_dtd = readtable('metab_dtd.csv',opts);

% Plot metabolism results from observed DO data
fig2 = figure(2);clf
fig2.WindowState = 'maximized';
plot(metab_obs.Date,metab_obs.Pg,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(metab_obs.Date,metab_obs.Rt,'k.-','MarkerSize',12,'LineWidth',1)
plot(metab_obs.Date,metab_obs.NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
title(['Pond ',num2str(pondNum),': R "ecometab" Results Using Observed Data'])
legend({'GPP','ER','NEM'},'FontSize',14)
% ylim([-1000 800])

% Plot metabolism results from detided DO data
fig3 = figure(3);clf
fig3.WindowState = 'maximized';
plot(metab_dtd.Date,metab_dtd.Pg,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(metab_dtd.Date,metab_dtd.Rt,'k.-','MarkerSize',12,'LineWidth',1)
plot(metab_dtd.Date,metab_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
title(['Pond ',num2str(pondNum),': R "ecometab" Results Using Detided Data'])
legend({'GPP','ER','NEM'},'FontSize',14)

% Save plots of R results
cd([rootpath,'figures\diel-analysis-figures\ponds\pond',num2str(pondNum)])
saveas(fig1,['pond',num2str(pondNum),'_DOcomparison.fig'])
saveas(fig1,['pond',num2str(pondNum),'_DOcomparison.png'])
saveas(fig2,['pond',num2str(pondNum),'_observed.fig'])
saveas(fig2,['pond',num2str(pondNum),'_observed.png'])
saveas(fig3,['pond',num2str(pondNum),'_detided.fig'])
saveas(fig3,['pond',num2str(pondNum),'_detided.png'])

%% Diel Analysis
pondNum = 6;
cd([rootpath,'diel-method\R-results\ponds\pond',num2str(pondNum)])

wtreg_res = readtable('wtreg_res.csv');
wtreg_res.Properties.DimensionNames{1} = 'datetime_utc';
wtreg_res.DateTimeStamp.TimeZone = "UTC";
wtreg_res = table2timetable(wtreg_res);

%====Define input variables================================================
dt_utc = wtreg_res.DateTimeStamp;
dt_local = dt_utc;
dt_local.TimeZone = 'America/New_York';

% Water
S = wtreg_res.Sal;
T = wtreg_res.Temp;
p = wtreg_res.Tide; % Pressure in dbar and depth in meters are approx. equal
d = wtreg_res.Tide; 
% H = mean(d,'omitnan') + 0.1524;   % Mean water depth [m] plus 6 in to account for height of sonde above bottom
H = 1.5 + 0.1524; % Mean water depth [m] plus 6 in to account for height of sonde above bottom

% Air
Tair = wtreg_res.ATemp;
patm = wtreg_res.BP;
u = wtreg_res.WSpd;

% DO concentration conversions
DO_obs = wtreg_res.DO_obs*1000/32;   % Observed DO concentration [mmol m-3]
DO_nrm = wtreg_res.DO_nrm*1000/32;   % Detided DO concentration [mmol m-3]

% Plot observed & detided DO concentration and tidal level (c.f. Beck Fig. 6)
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
t = tiledlayout(2,1,'TileSpacing','compact');
ax1 = nexttile;
plot(dt_local,DO_obs,'.-','MarkerSize',6,'Linewidth',1)
hold on
plot(dt_local,DO_nrm,'.-','MarkerSize',6,'Linewidth',1)
legend('Observed','Detided')
ylabel('DO conc. (mmol m^{-3})','FontSize',14)
title(['Pond ',num2str(pondNum)])
ax2 = nexttile;
plot(dt_local,d,'.-','MarkerSize',6,'Linewidth',1)
xlabel('Local Time')
ylabel('Water depth (m)')
linkaxes([ax1 ax2],'x')

% Conduct diel analysis with detided DO data
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

fig5 = figure(5);clf
fig5.WindowState = 'maximized';
plot(dt_local(dayind),DO_conc(dayind),'.b','MarkerSize',12)
hold on
plot(dt_local(nightind),DO_conc(nightind),'.k','MarkerSize',12)
ylabel('DO conc (mmol m^{-3})')
xlabel('Local Time')
title(['Pond ',num2str(pondNum),': Detided DO Concentration'])

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
title(['Pond ',num2str(pondNum),': Detided DO Concentration'])

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
title(['Pond ',num2str(pondNum),': Detided DO Concentration'])

% Daily rates
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

% Convert volumetric rates to depth-integrated (areal) estimates
GPP = P_daily * H;      % [mmol O2 m-2 d-1]
ER = R_daily * H;       % [mmol O2 m-2 d-1]
NEM = GPP + ER;         % [mmol O2 m-2 d-1]

diel_dtd = table(daystart_dt,dayend_dt,daylength,R_hourly,P_hourly,R_daily,P_daily,GPP,ER,NEM);
diel_dtd = table2timetable(diel_dtd);

fig8 = figure(8);clf
fig8.WindowState = 'maximized';
plot(diel_dtd.daystart_dt,diel_dtd.GPP,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(diel_dtd.daystart_dt,diel_dtd.ER,'k.-','MarkerSize',12,'LineWidth',1)
plot(diel_dtd.daystart_dt,diel_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1)
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title(['Pond ',num2str(pondNum),': MATLAB Results Using Detided Data'])

