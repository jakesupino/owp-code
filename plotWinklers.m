%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotWinklers.m
% This script...
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% January 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
cd([rootpath,'open-water-platform-data\gull\cleaned'])
load gull-cleaned.mat

cd([rootpath,'discrete-samples'])
wink = readtable('winklers_owp.csv');
varNames = ["datetime_utc","datetime_local","platform","S_lab","DO_mean","DO_std","DO_%err"];
varUnits = ["","","","psu","umol/L","umol/L","%"];
wink.Properties.VariableNames = varNames;
wink.Properties.VariableUnits = varUnits;
wink.datetime_local.TimeZone = 'America/New_York';
wink.datetime_utc.TimeZone = 'UTC';
wink = table2timetable(wink);

ind_gull = find(ismember(wink.platform,'Gull'));
ind_north = find(ismember(wink.platform,'North'));
ind_south = find(ismember(wink.platform,'South'));

%%
figure(1),clf
plot(sonde1_cleaned.datetime_utc,sonde1_cleaned.DO_conc,'.','MarkerSize',12)
hold on
errorbar(wink.datetime_utc(ind_gull),wink.DO_mean(ind_gull),wink.DO_std(ind_gull),'.','MarkerSize',12,'LineWidth',2)
ylabel('DO Concentration (\mumol/L)')
legend('Sonde Data','Winkler Data')
title('Gull - BC')

%%
t1 = wink.datetime_utc(ind_gull);
t2 = sonde1_cleaned.datetime_utc;
ind = interp1(t2,1:length(t2),t1,'nearest');
ind = rmmissing(ind);
%%
% Plot linear regression between Gull sonde and Winkler DO data
tbl = table(sonde1_cleaned.DO_conc(ind),wink.DO_mean(ind_gull(1:end-1)));
tbl.Properties.VariableNames = ["sonde","winkler"];
mdl = fitlm(tbl.sonde,tbl.winkler,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl.Rsquared.Ordinary,2);

figure(5),clf
h = plot(mdl,'marker','.');
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([0 ylim], [0 ylim],'--k')
errorbar(tbl.sonde,tbl.winkler,wink.DO_std(ind_gull(1:end-1)),'.b','MarkerSize',12,'LineWidth',2)
legend('Data',[eqn,newline,'R^2 = ',R2],'1:1 line')
xlabel('Sonde (Cleaned)')                                                                                                                                   
ylabel('Winkler')
title('Gull BC - DO Concentration')
daspect([1 1 1])
