clear;close all;clc

cd('G:\My Drive\Postdoc\Work\SMIIL\diel-method\matlab-results\gull-bc\movmed\')
load('diel_res.mat');

site = 'Gull';
sondename = 'BC';

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(diel_obs.daystart_dt,diel_obs.GPP,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(diel_obs.daystart_dt,diel_obs.ER,'k.-','MarkerSize',12,'LineWidth',1)
plot(diel_obs.daystart_dt,diel_obs.NEM,'r.-','MarkerSize',12,'LineWidth',1)
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,' ',sondename,' Sonde: MATLAB Results Using Observed Data'])
ylim([-1000 800])

fig2 = figure(2);clf
fig2.WindowState = 'maximized';
plot(diel_dtd.daystart_dt,diel_dtd.GPP,'.-','MarkerSize',12,'LineWidth',1)
hold on
plot(diel_dtd.daystart_dt,diel_dtd.ER,'k.-','MarkerSize',12,'LineWidth',1)
plot(diel_dtd.daystart_dt,diel_dtd.NEM,'r.-','MarkerSize',12,'LineWidth',1)
xlabel('UTC')
ylabel('mmol O_2 m^{-2} d^{-1}','FontSize',14)
legend({'GPP','ER','NEM'},'FontSize',14)
set(gca,'FontSize',14,'LineWidth',2)
title([site,' ',sondename,' Sonde: MATLAB Results Using Detided Data'])
ylim([-500 500])

cd('G:\My Drive\Postdoc\Work\SMIIL\figures\diel-analysis-figures\gull-bc\matlab-results\Moving Median')

%%
% Compute summary statistics
% Observed
meanGPP = mean(diel_obs.GPP,'omitmissing');
stdGPP = std(diel_obs.GPP,'omitmissing');
anomGPP = length(find(diel_obs.GPP < 0))/length(diel_obs.GPP)*100; % Percent of production estimates that are anomalous
meanER = mean(diel_obs.ER,'omitmissing');
stdER = std(diel_obs.ER,'omitmissing');
anomER = length(find(diel_obs.ER > 0))/length(diel_obs.ER)*100; % Percent of respiration estimates that are anomalous
meanNEM = mean(diel_obs.NEM,'omitmissing');
stdNEM = std(diel_obs.NEM,'omitmissing');
obs = table(meanGPP,stdGPP,anomGPP,meanER,stdER,anomER,meanNEM,stdNEM);

meanGPP = mean(diel_dtd.GPP,'omitmissing');
stdGPP = std(diel_dtd.GPP,'omitmissing');
anomGPP = length(find(diel_dtd.GPP < 0))/length(diel_dtd.GPP)*100; % Percent of production estimates that are anomalous
meanER = mean(diel_dtd.ER,'omitmissing');
stdER = std(diel_dtd.ER,'omitmissing');
anomER = length(find(diel_dtd.ER > 0))/length(diel_dtd.ER)*100; % Percent of respiration estimates that are anomalous
meanNEM = mean(diel_dtd.NEM,'omitmissing');
stdNEM = std(diel_dtd.NEM,'omitmissing');
dtd = table(meanGPP,stdGPP,anomGPP,meanER,stdER,anomER,meanNEM,stdNEM);

stats = struct("obs",obs,"dtd",dtd);
%%
cd('G:\My Drive\Postdoc\Work\SMIIL\open-water-platform-data\gull\cleaned\movmed')
save('summary-stats.mat','stats')