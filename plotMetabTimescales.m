%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotMetabTimescales.m
% This script plots the metabolic rates averaged across all sites and years
% by day of year and month of year.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 11/25/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import diel analysis results==========================================
fn = {'north','gull','south'};

for i = 1:length(fn)
    cd([rootpath,'diel-method\uncertainty-analysis\',fn{i}])
    load('MonteCarloResults')
    metab.(fn{i}) = diel_dtd_MC;

    dt2 = dateshift(metab.(fn{i}).daystart_dt,'start','day');
    metab.(fn{i}).daystart_dt = dt2;
    metab.(fn{i}) = removevars(metab.(fn{i}),{'daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});

    % Remove timepoints with anomalous GPP and ER values
    anomER = find(metab.(fn{i}).ER_avg > 0);
    anomGPP = find(metab.(fn{i}).GPP_avg < 0);
    metab.(fn{i})([anomER;anomGPP],:) = [];
end

clearvars diel_dtd diel_obs

%====Calculate day of year averages======================================
allDays = table;
allDays.day = (1:365)';
allDays.day = categorical(allDays.day);

% Find daily averages across years for each OWP
for i = 1:length(fn)
    temp = groupsummary(metab.(fn{i}),'date','dayofyear','mean');
    temp.Properties.VariableNames(1) = {'day'};

    dayofyear.(fn{i}) = outerjoin(allDays,temp);
    dayofyear.(fn{i}) = removevars(dayofyear.(fn{i}),{'day_temp'});
    dayofyear.(fn{i}).Properties.VariableNames(1) = {'day'};
end

% Find mean rates by day of year across the three OWPs
dayofyear_avg = table;
dayofyear_avg.day = (1:365)';
for i = 1:height(dayofyear_avg.day)
    dayofyear_avg.GPP(i) = mean([dayofyear.north.mean_GPP_avg(i) dayofyear.gull.mean_GPP_avg(i) dayofyear.south.mean_GPP_avg(i)],'omitmissing');
    dayofyear_avg.ER(i) = mean([dayofyear.north.mean_ER_avg(i) dayofyear.gull.mean_ER_avg(i) dayofyear.south.mean_ER_avg(i)],'omitmissing');
    dayofyear_avg.NEM(i) = mean([dayofyear.north.mean_NEM_avg(i) dayofyear.gull.mean_NEM_avg(i) dayofyear.south.mean_NEM_avg(i)],'omitmissing');

    dayofyear_avg.GPP_err(i) = mean([dayofyear.north.mean_GPP_sd(i) dayofyear.gull.mean_GPP_sd(i) dayofyear.south.mean_GPP_sd(i)],'omitmissing');
    dayofyear_avg.ER_err(i) = mean([dayofyear.north.mean_ER_sd(i) dayofyear.gull.mean_ER_sd(i) dayofyear.south.mean_ER_sd(i)],'omitmissing');
    dayofyear_avg.NEM_err(i) = mean([dayofyear.north.mean_NEM_sd(i) dayofyear.gull.mean_NEM_sd(i) dayofyear.south.mean_NEM_sd(i)],'omitmissing');

    % Check method
    % dayofyear_avg.GPP_lo(i) = mean([dayofyear.north.mean_GPP_avg(i)-dayofyear.north.mean_GPP_sd(i)...
    %     dayofyear.gull.mean_G PP_avg(i)-dayofyear.gull.mean_GPP_sd(i)...
    %     dayofyear.south.mean_GPP_avg(i)-dayofyear.south.mean_GPP_sd(i)],'omitmissing');
    % 
    % dayofyear_avg.ER_lo(i) = mean([dayofyear.north.mean_ER_avg(i)-dayofyear.north.mean_ER_sd(i)...
    %     dayofyear.gull.mean_ER_avg(i)-dayofyear.gull.mean_ER_sd(i)...
    %     dayofyear.south.mean_ER_avg(i)-dayofyear.south.mean_ER_sd(i)],'omitmissing');
    % 
    % dayofyear_avg.NEM_lo(i) = mean([dayofyear.north.mean_NEM_avg(i)-dayofyear.north.mean_NEM_sd(i)...
    %     dayofyear.gull.mean_NEM_avg(i)-dayofyear.gull.mean_NEM_sd(i)...
    %     dayofyear.south.mean_NEM_avg(i)-dayofyear.south.mean_NEM_sd(i)],'omitmissing');
end

%====Calculate month of year averages======================================
% Create monthly bins
monthly.gull = retime(metab.gull,'monthly','mean');
monthly.north = retime(metab.north,'monthly','mean');
monthly.south = retime(metab.south,'monthly','mean');

% Create monthly averages across years for each site/sonde
monthofyear.gull = groupsummary(metab.gull,'date','monthofyear','mean');
monthofyear.north = groupsummary(metab.north,'date','monthofyear','mean');
monthofyear.south = groupsummary(metab.south,'date','monthofyear','mean');

monthofyear.gull.Properties.VariableNames(1) = {'month'};
monthofyear.north.Properties.VariableNames(1) = {'month'};
monthofyear.south.Properties.VariableNames(1) = {'month'};

monthofyear_avg = table;
monthofyear_avg.month = (1:12)';
for i = 1:height(monthofyear_avg.month)
    monthofyear_avg.GPP(i) = mean([monthofyear.north.mean_GPP_avg(i) monthofyear.gull.mean_GPP_avg(i) monthofyear.south.mean_GPP_avg(i)]);
    monthofyear_avg.ER(i) = mean([monthofyear.north.mean_ER_avg(i) monthofyear.gull.mean_ER_avg(i) monthofyear.south.mean_ER_avg(i)]);
    monthofyear_avg.NEM(i) = mean([monthofyear.north.mean_NEM_avg(i) monthofyear.gull.mean_NEM_avg(i) monthofyear.south.mean_NEM_avg(i)]);

    monthofyear_avg.GPP_err(i) = mean([monthofyear.north.mean_GPP_sd(i) monthofyear.gull.mean_GPP_sd(i) monthofyear.south.mean_GPP_sd(i)]);
    monthofyear_avg.ER_err(i) = mean([monthofyear.north.mean_ER_sd(i) monthofyear.gull.mean_ER_sd(i) monthofyear.south.mean_ER_sd(i)]);
    monthofyear_avg.NEM_err(i) = mean([monthofyear.north.mean_NEM_sd(i) monthofyear.gull.mean_NEM_sd(i) monthofyear.south.mean_NEM_sd(i)]);

    % Check method
    % monthofyear_avg.GPP_lo(i) = mean([monthofyear.north.mean_GPP_avg(i)-monthofyear.north.mean_GPP_sd(i)...
    %     monthofyear.gull.mean_GPP_avg(i)-monthofyear.gull.mean_GPP_sd(i)...
    %     monthofyear.south.mean_GPP_avg(i)-monthofyear.south.mean_GPP_sd(i)],'omitmissing');
    % 
    % monthofyear_avg.ER_lo(i) = mean([monthofyear.north.mean_ER_avg(i)-monthofyear.north.mean_ER_sd(i)...
    %     monthofyear.gull.mean_ER_avg(i)-monthofyear.gull.mean_ER_sd(i)...
    %     monthofyear.south.mean_ER_avg(i)-monthofyear.south.mean_ER_sd(i)],'omitmissing');
    % 
    % monthofyear_avg.NEM_lo(i) = mean([monthofyear.north.mean_NEM_avg(i)-monthofyear.north.mean_NEM_sd(i)...
    %     monthofyear.gull.mean_NEM_avg(i)-monthofyear.gull.mean_NEM_sd(i)...
    %     monthofyear.south.mean_NEM_avg(i)-monthofyear.south.mean_NEM_sd(i)],'omitmissing');
end

%% ====Plots of metabolic rates on daily & monthly time scales=============
% Colors from https://personal.sron.nl/~pault/
gpp_clr = hex2rgb('#55AA22');
er_clr = hex2rgb('#BB0011');
nem_clr = hex2rgb('#88551a'); % Mix of gpp and er colors (https://www.w3schools.com/colors/colors_mixer.asp)

% Version with using shaded errorbar for Day of Year plots
% https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar?s_tid=srchtitle

fig(1) = figure(1);clf
t = tiledlayout(2,2,'TileSpacing','compact','padding','compact');

% Tile 1
nexttile
% GPP
x = dayofyear_avg.day;
y = dayofyear_avg.GPP;
yerr = dayofyear_avg.GPP_err;
s = shadedErrorBar(x,y,yerr,'lineprops','-','patchsaturation',0.33);
set(s.edge,'LineStyle','none','Color',gpp_clr)
s.mainLine.LineWidth = 1;
s.mainLine.Color = gpp_clr;
s.patch.FaceColor = gpp_clr;
hold on
% ER
y = dayofyear_avg.ER;
yerr = dayofyear_avg.ER_err;
s = shadedErrorBar(x,y,yerr,'lineprops','-','patchsaturation',0.33);
set(s.edge,'LineStyle','none','Color',er_clr)
s.mainLine.LineWidth = 1;
s.mainLine.Color = er_clr;
s.patch.FaceColor = er_clr;

yline(0,'-k','LineWidth',2)
text(0.03,0.9,'GPP','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',gpp_clr)
text(0.03,0.05,'ER','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',er_clr)
set(gca,'XTickLabel',[])
ylabel('mmol O_2 m^{-2} d^{-1}','fontsize',14)
title('Daily')
set(gca,'box','off','LineWidth',1.5)

% Tile 2
nexttile
b1 = bar(monthofyear_avg.month,monthofyear_avg.GPP,0.5,'FaceColor',gpp_clr,'EdgeColor',gpp_clr);
hold on
for k = 1:numel(b1)
    xtips = b1(k).XEndPoints;
    ytips = b1(k).YEndPoints;
    errorbar(xtips,ytips,monthofyear_avg.GPP_err,monthofyear_avg.GPP_err,'.k','MarkerSize',0.1)
end
b2 = bar(monthofyear_avg.month,monthofyear_avg.ER,0.5,'FaceColor',er_clr,'EdgeColor',er_clr);
for k = 1:numel(b2)
    xtips = b2(k).XEndPoints;
    ytips = b2(k).YEndPoints;
    errorbar(xtips,ytips,monthofyear_avg.ER_err,monthofyear_avg.ER_err,'.k','MarkerSize',0.1)
end
text(0.03,0.9,'GPP','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',gpp_clr)
text(0.03,0.05,'ER','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',er_clr)
set(gca,'XTickLabel',[])
title('Monthly')
set(gca,'box','off','LineWidth',1.5)
ylim([-400 300])

% Tile 3
nexttile
y = dayofyear_avg.NEM;
yerr = dayofyear_avg.NEM_err;
s = shadedErrorBar(x,y,yerr,'lineprops','-','patchsaturation',0.33);
set(s.edge,'LineStyle','none','Color',nem_clr)
s.mainLine.LineWidth = 1;
s.mainLine.Color = nem_clr;
s.patch.FaceColor = nem_clr;

yline(0,'-k','LineWidth',2)
text(0.03,0.05,'NEM','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',nem_clr)
set(gca,'box','off','LineWidth',1.5)
ylabel('mmol O_2 m^{-2} d^{-1}','fontsize',14)
xlabel('Day of Year','fontsize',14)

% Tile 4
nexttile
bar(monthofyear_avg.month,monthofyear_avg.NEM,0.5,'k')
hold on
b3 = bar(monthofyear_avg.month,monthofyear_avg.NEM,0.5,'FaceColor',nem_clr,'EdgeColor',nem_clr);
for k = 1:numel(b3)
    xtips = b3(k).XEndPoints;
    ytips = b3(k).YEndPoints;
    errorbar(xtips,ytips,monthofyear_avg.NEM_err,monthofyear_avg.NEM_err,'.k','MarkerSize',0.1)
end
text(0.03,0.05,'NEM','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',nem_clr)
set(gca,'box','off','LineWidth',1.5)
xlabel('Month of Year','fontsize',14)
% plot(monthofyear_avg.month,monthofyear_avg.NEM_lo,'--')   % Check method

% set(gcf,'Position',[-1200 -100 400 300])
fig.Units               = 'centimeters';
fig.Position(1) = 0;
fig.Position(2) = 0;
fig.Position(3)         = 25.2095 ;
fig.Position(4)         = 18;

%====Save the plot=========================================================
option = questdlg('Save plot?','Save plot','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\diel-analysis\matlab-results\all-sites'])
        saveas(fig(1),'metab_daily-monthly.fig')
        exportgraphics(fig(1),'metab_daily-monthly.jpg','Resolution',600)
        disp('Plot saved!')
    case 'No'
        disp('Plot not saved.')
end

%% Old version using bar for Day of Year plots
% fig = figure;clf
% t = tiledlayout(2,2,'TileSpacing','compact','padding','compact');
%
% % Tile 1
% nexttile
% bar(dayofyear_avg.day,dayofyear_avg.GPP,0.5,'FaceColor',gpp_clr,'EdgeColor',gpp_clr)
% hold on
% bar(dayofyear_avg.day,dayofyear_avg.ER,0.5,'FaceColor',er_clr,'EdgeColor',er_clr)
% text(0.03,0.9,'GPP','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',gpp_clr)
% text(0.03,0.05,'ER','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',er_clr)
% set(gca,'XTickLabel',[])
% ylabel('mmol O_2 m^{-2} d^{-1}','fontsize',14)
% title('Daily')
% set(gca,'box','off','LineWidth',1.5)
%
% % Tile 2
% nexttile
% b1 = bar(monthofyear_avg.month,monthofyear_avg.GPP,0.5,'FaceColor',gpp_clr,'EdgeColor',gpp_clr);
% hold on
% for k = 1:numel(b1)
%     xtips = b1(k).XEndPoints;
%     ytips = b1(k).YEndPoints;
%     errorbar(xtips,ytips,monthofyear_avg.GPP_err,monthofyear_avg.GPP_err,'.k','MarkerSize',0.1)
% end
% b2 = bar(monthofyear_avg.month,monthofyear_avg.ER,0.5,'FaceColor',er_clr,'EdgeColor',er_clr);
% for k = 1:numel(b2)
%     xtips = b2(k).XEndPoints;
%     ytips = b2(k).YEndPoints;
%     errorbar(xtips,ytips,monthofyear_avg.ER_err,monthofyear_avg.ER_err,'.k','MarkerSize',0.1)
% end
% text(0.03,0.9,'GPP','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',gpp_clr)
% text(0.03,0.05,'ER','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',er_clr)
% set(gca,'XTickLabel',[])
% title('Monthly')
% set(gca,'box','off','LineWidth',1.5)
%
% % Tile 3
% nexttile
% bar(dayofyear_avg.day,dayofyear_avg.NEM,0.5,'FaceColor',nem_clr,'EdgeColor',nem_clr)
% text(0.03,0.05,'NEM','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',nem_clr)
% set(gca,'box','off','LineWidth',1.5)
% ylabel('mmol O_2 m^{-2} d^{-1}','fontsize',14)
% xlabel('Day of Year','fontsize',14)
%
% % Tile 4
% nexttile
% bar(monthofyear_avg.month,monthofyear_avg.NEM,0.5,'k')
% hold on
% b3 = bar(monthofyear_avg.month,monthofyear_avg.NEM,0.5,'FaceColor',nem_clr,'EdgeColor',nem_clr);
% for k = 1:numel(b3)
%     xtips = b3(k).XEndPoints;
%     ytips = b3(k).YEndPoints;
%     errorbar(xtips,ytips,monthofyear_avg.NEM_err,monthofyear_avg.NEM_err,'.k','MarkerSize',0.1)
% end
% text(0.03,0.05,'NEM','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',nem_clr)
% set(gca,'box','off','LineWidth',1.5)
% xlabel('Month of Year','fontsize',14)
%
% set(gcf,'Position',[-1200 -100 400 300])
% fig.Units               = 'centimeters';
% fig.Position(3)         = 25.2095 ;
% fig.Position(4)         = 18;
