%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotWinklers.m
% This script...
%
% AUTHOR:
% Emily Chua 
% 
% First created: 1/2024
% Last updated: 4/16/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import Winkler data===================================================
cd([rootpath,'discrete-samples'])
wink = readtable('winklers_owp.csv');
varNames = ["datetime_utc","datetime_local","platform","S_lab","DO_mean","DO_std","DO_%err"];
varUnits = ["","","","psu","umol/L","umol/L","%"];
wink.Properties.VariableNames = varNames;
wink.Properties.VariableUnits = varUnits;
wink.datetime_local.TimeZone = 'America/New_York';
wink.datetime_utc.TimeZone = 'UTC';
wink = table2timetable(wink);

%====Import sonde data=====================================================
% Gull
cd([rootpath,'open-water-platform-data\gull\cleaned\movmed'])
load gull-bc-cleaned.mat
load gull-erdc-cleaned.mat
gull.sonde1 = sonde1_cleaned;
gull.sonde2 = sonde2_cleaned;

clearvars sonde1_cleaned sonde2_cleaned

% North
% cd([rootpath,'open-water-platform-data\north\cleaned\movmed'])
cd([rootpath,'open-water-platform-data\north\cleaned\initial-qc'])
load north-bc-cleaned.mat
load north-erdc-cleaned.mat
north.sonde1 = sonde1_cleaned;
north.sonde2 = sonde2_cleaned;

clearvars sonde1_cleaned sonde2_cleaned

% South
% cd([rootpath,'open-water-platform-data\south\cleaned\movmed'])
cd([rootpath,'open-water-platform-data\south\cleaned\initial-qc'])
load south-bc-cleaned.mat
load south-erdc-cleaned.mat
south.sonde1 = sonde1_cleaned;
south.sonde2 = sonde2_cleaned;

clearvars sonde1_cleaned sonde2_cleaned

%% Find closest matching times between Winkler and sonde data
%====Find closest matching indices=========================================
winkler = find(ismember(wink.platform,'Gull'));
t1 = wink.datetime_utc(winkler);
t2 = gull.sonde1.datetime_utc;
bc = interp1(t2,1:length(t2),t1,'nearest','extrap');

% If closest matching samples are more than 20 minutes off, do not include
for i = 1:length(bc(~isnan(bc)))
    dt = between(gull.sonde1.datetime_utc(bc(i)),wink.datetime_utc(winkler(i)));
    if abs(time(dt)) > minutes(20)
        bc(i) = NaN;
        % winkler(i) = NaN;
    end
end

ind_bc = table(winkler,bc);

% Gull ERDC
winkler = find(ismember(wink.platform,'Gull'));
t1 = wink.datetime_utc(winkler);
t2 = gull.sonde2.datetime_utc;
erdc = interp1(t2,1:length(t2),t1,'nearest');

% If closest matching samples are more than 20 minutes off, do not include
for i = 1:length(erdc(~isnan(erdc)))
    dt = between(gull.sonde2.datetime_utc(erdc(i)),wink.datetime_utc(winkler(i)));
    if abs(time(dt)) > minutes(20)
        erdc(i) = NaN;
        % winkler(i) = NaN;
    end
end

ind_erdc = table(winkler,erdc);

ind_gull = join(ind_bc,ind_erdc);

%%
%====Global plotting settings==============================================
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 12;
circlesize = 6;

% Plot DO timeseries data for both sondes, with Winkler samples on top
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(gull.sonde1.datetime_utc,gull.sonde1.DO_conc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC')
hold on
plot(gull.sonde2.datetime_utc,gull.sonde2.DO_conc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
errorbar(wink.datetime_utc(ind_gull.winkler),wink.DO_mean(ind_gull.winkler),wink.DO_std(ind_gull.winkler),'.','MarkerSize',dotsize,'LineWidth',2,'DisplayName','Winkler')
ylabel('DO Concentration (\mumol/L)')
legend('show')
title('Gull Dissolved Oxygen Validation')

% fig2 = figure(2);clf
% fig2.WindowState = 'maximized';
% plot(north.sonde1.datetime_utc,north.sonde1.DO_conc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC')
% hold on
% plot(north.sonde2.datetime_utc,north.sonde2.DO_conc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
% errorbar(wink.datetime_utc(ind_north),wink.DO_mean(ind_north),wink.DO_std(ind_north),'.','MarkerSize',dotsize,'LineWidth',2,'DisplayName','Winkler')
% ylabel('DO Concentration (\mumol/L)')
% legend('show')
% title('North Dissolved Oxygen Validation')
% 
% fig3 = figure(3);clf
% fig3.WindowState = 'maximized';
% plot(south.sonde1.datetime_utc,south.sonde1.DO_conc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC')
% hold on
% plot(south.sonde2.datetime_utc,south.sonde2.DO_conc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
% errorbar(wink.datetime_utc(ind_south),wink.DO_mean(ind_south),wink.DO_std(ind_south),'.','MarkerSize',dotsize,'LineWidth',2,'DisplayName','Winkler')
% ylabel('DO Concentration (\mumol/L)')
% legend('show')
% title('South Dissolved Oxygen Validation')

%% %====Linear Regressions=================================================
% Gull BC linear regression
include = ~isnan(ind_gull.bc);
x = gull.sonde1.DO_conc(ind_gull.bc(include));
y = wink.DO_mean(ind_gull.winkler(include));
mdl_bc = fitlm(x,y,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl_bc.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl_bc.Rsquared.Ordinary,2);

figure(5),clf
h = plot(mdl_bc,'marker','.');
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([0 ylim], [0 ylim],'--k')
errorbar(gull.sonde1.DO_conc(ind_gull.bc(include)),wink.DO_mean(ind_gull.winkler(include)),wink.DO_std(ind_gull.winkler(include)),'.b','MarkerSize',dotsize,'LineWidth',2)
legend('Data',[eqn,newline,'R^2 = ',R2],'1:1 line')
xlabel('BC DO Concentration (\mumol/L)')                                                                                                                                   
ylabel('Winkler DO Concentration (\mumol/L)')
title('Gull BC - DO Validation')
daspect([1 1 1])

% Gull ERDC linear regression
include = ~isnan(ind_gull.erdc);
x = gull.sonde2.DO_conc(ind_gull.erdc(include));
y = wink.DO_mean(ind_gull.winkler(include));
mdl_erdc = fitlm(x,y,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl_erdc.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl_erdc.Rsquared.Ordinary,2);

figure(6),clf
h = plot(mdl_erdc,'marker','.');
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([0 400], [0 400],'--k')
errorbar(gull.sonde2.DO_conc(ind_gull.erdc(include)),wink.DO_mean(ind_gull.winkler(include)),wink.DO_std(ind_gull.winkler(include)),'.b','MarkerSize',dotsize,'LineWidth',2)
legend('Data',[eqn,newline,'R^2 = ',R2],'1:1 line')
xlabel('ERDC DO Concentration (\mumol/L)')                                                                                                                                   
ylabel('Winkler DO Concentration (\mumol/L)')
title('Gull ERDC - DO Validation')
daspect([1 1 1])
