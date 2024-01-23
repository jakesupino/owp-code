%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dielAnalysis.m
% This script uses the diel oxygen method to calculate GPP, ER, and NEM
% following Beck et al. (2015)
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% January 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
cd([rootpath,'open-water-platform-data\gull\cleaned'])
load gull-cleaned.mat

cd([rootpath,'physical-data\wind-speed'])
load windSpeed.mat

cd([rootpath,'physical-data\par'])
load par.mat

% Retime cleaned Gull met data to same datetimes as sonde data 
newTimes = sonde1_cleaned.datetime_utc(1):minutes(10):sonde1_cleaned.datetime_utc(end);
metDat_rt = retime(metDat_cleaned,newTimes,'mean'); 

%====Find where there are gaps in Gull met station wind speed data=========
% Retime NOAA data to same datetimes as sonde data 
noaaDat_rt = retime(noaaDat,newTimes,'previous');

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

%====Gap fill Gull met station wind speed data using NOAA daily means======
ind_nan = find(isnan(metDat_rt.wspd));
metDat_rt.wspd(ind_nan) = noaaDat_rt.wspd_avg(ind_nan);

%====Find where there are gaps in Gull met station air T data==============
% Retime PAR data to same datetimes as sonde data 
parDat_rt = retime(parDat,newTimes,'previous');

for i = 1:height(metDat_rt.datetime_utc)
    nid = find(isnan(metDat_rt.Tair));
    id(1:length(nid),1) = nid;
end

figure(2),clf
plot(metDat_rt.datetime_utc,metDat_rt.Tair,'.-')
hold on
plot(parDat_rt.datetime_utc(id),parDat_rt.Tair(id),'or','MarkerSize',4)
ylabel('Air T (^oC)')
legend('Gull Met Station','PAR Dataset')
title('Gap-Filled Air Temperature Data')

%====Gap fill Gull met station air T data using PAR data===================
ind_nan = find(isnan(metDat_rt.Tair));
metDat_rt.Tair(ind_nan) = parDat_rt.Tair(ind_nan);

%====Gap fill Gull met station atm p data using global mean value==========
patm_mean = mean(metDat_cleaned.patm,'omitmissing');
ind_nan = find(isnan(metDat_rt.patm));
metDat_rt.patm(ind_nan) = patm_mean;

figure(3),clf
plot(metDat_rt.datetime_utc,metDat_rt.patm,'.-')
hold on
plot(metDat_rt.datetime_utc(ind_nan),metDat_rt.patm(ind_nan),'or','MarkerSize',4)
ylabel('p_{atm} (hPa)')
legend('Gull Met Station','Gull Met Station Mean')
title('Gap-Filled Atmospheric Pressure Data')

% Horizontally concatenate sonde and met data (already have common time vector)
dat = synchronize(sonde1_cleaned,metDat_rt);
dat = rmmissing(dat,'DataVariables',{'DO_conc','wspd'});

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
D = ka.*(DO_sat - dat.DO_conc); % [mmol m-3 h-1]

%% STEP 3
%====Apply filter========================================================
% Needoba uses HBF
% Caffrey et al. (2014) doesn't bother with filtering
% Beck et al. (2015) uses their R package to remove tidal advection

%% STEP 4
%====Calculate rates=====================================================
dCdt = nan(length(dat.DO_conc),1);
dCdt(2:end,1) = diff(dat.DO_conc) ./ hours(diff(dat.datetime_utc));  % [mmol m-3 h-1]

% Get day/night indices based on lat and lon and Hilary's indexDayNight function
% Or use PAR data for SMIIL analysis??
lat = 39.08;
lon = -74.78;
time_in = dat.datetime_utc;
tol = 0;

% Manually define DST start and end times
start1 = datenum('03/14/2021 02:00','mm/dd/yyyy HH:MM');
end1 = datenum('11/07/2021 02:00','mm/dd/yyyy HH:MM');
start2 = datenum('03/13/2022 02:00','mm/dd/yyyy HH:MM');
end2 = datenum('11/6/2022 02:00','mm/dd/yyyy HH:MM');
start3 = datenum('03/12/2023 02:00','mm/dd/yyyy HH:MM');
end3 = datenum('11/5/2023 02:00','mm/dd/yyyy HH:MM');

% Calculate all day/night indices assuming EDT
UTCoffset = -4; % [h]
[dayind_edt,nightind_edt] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);
ind_edt = [dayind_edt;nightind_edt];

% Calculate all day/night indices assuming EST
UTCoffset = -5; % [h]
[dayind_est,nightind_est] = indexDayNight_EJC(lat,lon,UTCoffset,time_in,tol);

% Find closest indices in datetime_utc that match DST start and end times
[~,edt_end1] = min(abs(datenum(dat.datetime_utc)-end1));    % Note: This date isn't close due to the data gap in fall 2021/winter 2022
[~,edt_start2] = min(abs(datenum(dat.datetime_utc)-start2));
[~,edt_end2] = min(abs(datenum(dat.datetime_utc)-end2));
[~,edt_start3] = min(abs(datenum(dat.datetime_utc)-start3));
[~,edt_end3] = min(abs(datenum(dat.datetime_utc)-end3));

dayind_edt = table(dayind_edt,repmat({'day'},[length(dayind_edt),1]));
dayind_edt.Properties.VariableNames = ["index","day/night"];
nightind_edt = table(nightind_edt,repmat({'night'},[length(nightind_edt),1]));
nightind_edt.Properties.VariableNames = ["index","day/night"];
ind_edt = [dayind_edt; nightind_edt];
ind_edt = sortrows(ind_edt);

dayind_est = table(dayind_est,repmat({'day'},[length(dayind_est),1]));
dayind_est.Properties.VariableNames = ["index","day/night"];
nightind_est = table(nightind_est,repmat({'night'},[length(nightind_est),1]));
nightind_est.Properties.VariableNames = ["index","day/night"];
ind_est = [dayind_est; nightind_est];
ind_est = sortrows(ind_est);

% Define periods for whether it is DST or not
t1 = ind_edt(1:edt_end1,:);              % EDT period
t2 = ind_est(edt_end1+1:edt_start2-1,:); % EST period
t3 = ind_edt(edt_start2:edt_end2,:);     % EDT period
t4 = ind_est(edt_end2+1:edt_start3-1,:);   % EST period
t5 = ind_edt(edt_start3:end,:);          % EDT period

allind = [t1;t2;t3;t4;t5];

dayind = strcmp(allind.("day/night"),'day');
nightind = strcmp(allind.("day/night"),'night');

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
% xline(end3_dt,'--','label','Daylight Saving Time Ends','labelhorizontalalignment','left')
ylabel('DO conc (\mumol/L)')

% Find the indices for when each day starts and stops
daystart = find(diff(dayind) == 1) + 1;
dayend = find(diff(dayind) == -1);

% Length of each day
daylength = dat.datetime_utc(dayend(2:end)) - dat.datetime_utc(daystart(1:end));
daylength = hours(daylength);

% figure,clf;plot(dat.datetime_utc(daystart(2:end)),daylength,'.')

% Hourly rates of respiration and net production
for i = 1:length(daylength)
    % During night hours, P = 0
    R_hourly(i,1) = mean(dCdt(dayend(i):daystart(i)) - D(dayend(i):daystart(i)),'omitnan');     % Hourly rate of respiration; [mmol m-3 h-1]
    % During day hours, P != 0
    P_hourly(i,1) = mean(dCdt(daystart(i):dayend(i+1)) - D(daystart(i):dayend(i+1)),'omitnan'); % Hourly rate of net production; [mmol m-3 h-1]
end

% Daily rates
R_daily = R_hourly .* 24;                      % Daily rate of respiration; [mmol m-3 d-1]
P_daily = (P_hourly - R_hourly) .* daylength;  % Daily rate of gross production; [mmol m-3 d-1]

% Convert volumetric rates to depth-integrated (aereal) estimates
GPP = P_daily * H;      % [mmol O2 m-2 d-1]
ER = R_daily * H;       % [mmol O2 m-2 d-1]
NEM = GPP + ER;         % [mmol O2 m-2 d-1]

% Plot results
figure(1),clf
plot(dat.datetime_utc(daystart),GPP,'.-','MarkerSize',12,'LineWidth',2)
hold on
plot(dat.datetime_utc(daystart),ER,'k.-','MarkerSize',12,'LineWidth',2)
plot(dat.datetime_utc(daystart),NEM,'r.-','MarkerSize',12,'LineWidth',2)
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title('Gull - BC Sonde')

cd('G:\My Drive\Postdoc\Work\SMIIL\figures\diel-analysis-figures')