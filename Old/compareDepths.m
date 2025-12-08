clear; close all; clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

figure(1),clf

cd([rootpath,'open-water-platform-data\gull\cleaned\final-qc'])
load('Gull-cleaned.mat')
datetime_local = finalQC.datetime_utc;
datetime_local.TimeZone = 'America/New_York';
plot(datetime_local,finalQC.depth,'.','Color',rgb('goldenrod'),'DisplayName','Gull')

hold on

cd([rootpath,'open-water-platform-data\north\cleaned\final-qc'])
load('north-cleaned.mat')
datetime_local = finalQC.datetime_utc;
datetime_local.TimeZone = 'America/New_York';
plot(datetime_local,finalQC.depth,'.','Color',rgb('lightseagreen'),'DisplayName','North')

cd([rootpath,'open-water-platform-data\south\cleaned\final-qc'])
load('south-cleaned.mat')
datetime_local = finalQC.datetime_utc;
datetime_local.TimeZone = 'America/New_York';
plot(datetime_local,finalQC.depth,'.','Color',rgb('mediumpurple'),'DisplayName','South')

legend('show')

%%
cd('G:\My Drive\Postdoc\Work\SMIIL\open-water-platform-data\south\original\merged')
load('alldeps-south-raw.mat')
bc_orig = sonde1_all;
erdc_orig = sonde2_all;

cd('G:\My Drive\Postdoc\Work\SMIIL\open-water-platform-data\south\adjusted\merged')
load('alldeps-south-adj.mat')
bc_adj = sonde1_all;
erdc_adj = sonde2_all;

clear sonde1_all sonde2_all

figure,clf
plot(bc_orig.datetime_local,bc_orig.depth,'.','DisplayName','Original')
hold on
plot(bc_adj.datetime_local,bc_adj.depth,'.k','DisplayName','Adjusted')
ylabel('Depth (m)')
xlabel('Local time')
legend('show')
