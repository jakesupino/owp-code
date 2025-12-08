%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% statsAnalysis_pH_DO_coupling.m
% This script...
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 7/22/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import final QC'd data================================================
site = 'gull';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
gull_orig = finalQC;

site = 'north';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
north_orig = finalQC;

site = 'south';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
south_orig = finalQC;

clearvars finalQC
%%
%====Predict monthly mean of pH from monthly means of DO and S=============
data.north = north_orig;
data.gull = gull_orig;
data.south = south_orig;

fn = fieldnames(data);
siteName = {'North','Gull','South'};

cd([rootpath,'figures\stats-analyses'])
fig = figure(1);clf
fig.WindowState = 'maximized';
tl = tiledlayout(1,3,'TileSpacing','compact');

for i = 1:numel(fn)
    monAvg = retime(data.(fn{i})(:,2:end),'monthly','mean');
    monAvg = timetable2table(monAvg);
    
    % Fit model to monthly means of DO%sat and salinity
    mdl.(fn{i}) = fitlm(monAvg,'pH ~ DOsat + salinity');
    % gull_mdl = fitlm(gull_monAvg,'pH ~ DOsat');

    xobs = table(monAvg.salinity,monAvg.DOsat,'VariableNames',{'salinity','DOsat'});
    pHpred = predict(mdl.(fn{i}),xobs);

    % Plot observed mean pH vs. predicted mean pH
    nexttile
    plot(pHpred,monAvg.pH,'.k','markersize',12)
    set(gca,'xlim',[7.6 8.5]);
    set(gca,'ylim',[7.6 8.5]);
    xl = get(gca,'xlim');
    yl = xl;
    hold on
    plot(xl,yl,'--','color',rgb('lightgrey'))
    daspect([1 1 1])
    t1 = text(xl(1)+.05,yl(2)-.03,['R_{adj}^2 = ',num2str(mdl.(fn{i}).Rsquared.Adjusted,'%.3f')],'fontsize',12);
    t1.HorizontalAlignment = 'left';
    t1.VerticalAlignment = 'top';

    % Only include salinity if p < 0.05
    if mdl.(fn{i}).Coefficients.pValue(2) < 0.05
        t2 = text(xl(1)+.05,yl(2)-.1,['pH = ',num2str(mdl.(fn{i}).Coefficients.Estimate(1),'%.3f'),' + ',...
            num2str(mdl.(fn{i}).Coefficients.Estimate(2),'%.3f'),'*S + ',...
            num2str(mdl.(fn{i}).Coefficients.Estimate(3),'%.3f'),'*DO%sat'],'fontsize',12);
    else
        t2 = text(xl(1)+.05,yl(2)-.1,['pH = ',num2str(mdl.(fn{i}).Coefficients.Estimate(1),'%.3f'),' + ',...
            num2str(mdl.(fn{i}).Coefficients.Estimate(3),'%.3f'),'*DO%sat'],'fontsize',12);
    end
    t2.HorizontalAlignment = 'left';
    t2.VerticalAlignment = 'top';
    title(siteName(i))

end
xlabel(tl,'Predicted monthly mean pH','fontsize',16)
ylabel(tl,'Observed monthly mean pH','fontsize',16)

% Combined data
all_data = synchronize(data.north,data.gull,data.south);

pH_mat = [all_data.pH_1,all_data.pH_2,all_data.pH_3];
S_mat = [all_data.salinity_1,all_data.salinity_2,all_data.salinity_3];
DO_mat = [all_data.DOsat_1,all_data.DOsat_2,all_data.DOsat_3];
% Calculate means for each row
all_data.mean_pH = mean(pH_mat,2,'omitmissing');
all_data.mean_S = mean(S_mat,2,'omitmissing');
all_data.mean_DOsat = mean(DO_mat,2,'omitmissing');
% Calculate monthly means
monAvg = retime(all_data,'monthly','mean');
monAvg = timetable2table(monAvg);
% Fit model to monthly means of DO%sat and salinity
mdl.all = fitlm(monAvg,'mean_pH ~ mean_DOsat + mean_S');

% Method 1: Predict pH from monthly across-site means
% xobs = table(monAvg.mean_S,monAvg.mean_DOsat,'VariableNames',{'mean_S','mean_DOsat'});
% pHpred = predict(mdl,xobs);

% Method 2: Predict pH from monthly site-specific means
% North
monAvg_north = retime(data.north(:,2:end),'monthly','mean');
northobs = table(monAvg_north.salinity,monAvg_north.DOsat,'VariableNames',{'mean_S','mean_DOsat'});
pHpred_north = predict(mdl.all,northobs);
% Gull
monAvg_gull = retime(data.gull(:,2:end),'monthly','mean');
gullobs = table(monAvg_gull.salinity,monAvg_gull.DOsat,'VariableNames',{'mean_S','mean_DOsat'});
pHpred_gull = predict(mdl.all,gullobs);
% South
monAvg_south = retime(data.south(:,2:end),'monthly','mean');
southobs = table(monAvg_south.salinity,monAvg_south.DOsat,'VariableNames',{'mean_S','mean_DOsat'});
pHpred_south = predict(mdl.all,southobs);

monAvg_gull = retime(data.gull(:,2:end),'monthly','mean');
monAvg_south = retime(data.south(:,2:end),'monthly','mean');

cd([rootpath,'figures\stats-analyses'])
fig = figure;clf
% Plot observed mean pH vs. predicted mean pH
% plot(pHpred,monAvg.mean_pH,'.k','markersize',12)  % Method 1: Plot across-site means
plot(pHpred_north,monAvg_north.pH,'.','Color',rgb('lightseagreen'),'DisplayName','North')
hold on
plot(pHpred_gull,monAvg_gull.pH,'.','Color',rgb('goldenrod'),'DisplayName','Gull')
plot(pHpred_south,monAvg_south.pH,'.','Color',rgb('mediumpurple'),'DisplayName','South')
set(gca,'xlim',[7.6 8.5]);
set(gca,'ylim',[7.6 8.5]);
xl = get(gca,'xlim');
yl = xl;
plot(xl,yl,'--','color',rgb('lightgrey'),'DisplayName','1:1 line') 
daspect([1 1 1])
t1 = text(xl(1)+.02,yl(2)-.01,['R_{adj}^2 = ',num2str(mdl.all.Rsquared.Adjusted,'%.3f')],'fontsize',12);
t1.HorizontalAlignment = 'left';
t1.VerticalAlignment = 'top';

% Only include salinity in regression if p < 0.05
if mdl.all.Coefficients.pValue(2) < 0.05
    t2 = text(xl(1)+.02,yl(2)-.053,['pH = ',num2str(mdl.all.Coefficients.Estimate(1),'%.3f'),' + ',...
        num2str(mdl.all.Coefficients.Estimate(2),'%.3f'),'*S + ',...
        num2str(mdl.all.Coefficients.Estimate(3),'%.3f'),'*DO%sat'],'fontsize',12);
else
    t2 = text(xl(1)+.02,yl(2)-.05,['pH = ',num2str(mdl.(fn{i}).Coefficients.Estimate(1),'%.3f'),' + ',...
        num2str(mdl.all.(fn{i}).Coefficients.Estimate(3),'%.3f'),'*DO%sat'],'fontsize',12);
end
t2.HorizontalAlignment = 'left';
t2.VerticalAlignment = 'top';
title('All sites')
xlabel('Predicted monthly mean pH','fontsize',16)
ylabel('Observed monthly mean pH','fontsize',16)
legend('show','location','southeast')
daspect([1 1 1])