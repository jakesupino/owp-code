%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stormEvents.m
% This script picks out major storms during the 3-year OWP measurement period and 
% evaluates whether NEM values differ significantly during a storm event
% and the pre- and post-storm periods, using the Monte Carlo mean NEM
% calculated in "dielAnalysis_montecarlo.m".
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 8/21/2024
% Last updated: 12/17/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import windspeed and PAR data=====================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

cd([rootpath,'physical-data\final-dataset'])
load('par.mat')

%====Calculate daily means for each site===============================
wspd_dailyAvg = retime(era5Dat,'daily','mean');
par_dailyAvg = retime(parDat,'daily','mean');

%==========================================================================
% Run analyses on individual sites and save results in "NEM_tbl_all" for
% pooled t-test (after for loop)
%==========================================================================
site = {'north','gull','south'};
siteNames = {'North','Gull','South'};

pre_clr = [187 187 187]/255;
event_clr = [238 102 119]/255;
post_clr = [187 187 187]/255;

colors(1,:) = pre_clr;
colors(2,:) = event_clr;
colors(3,:) = post_clr;

north_clr = '#5D3A9B';
gull_clr = '#019E73';
south_clr = '#D55E00';

for kk = 1:length(site)
    clearvars start_pre start_stm start_post end_pre end_stm end_post pre_NEM event_NEM post_NEM

    %====Import final QC'd data & diel analysis results========================
    cd([rootpath,'open-water-platform-data\',site{kk},'\cleaned\final-qc'])
    load('finalQC.mat');
    params = finalQC;
    cd([rootpath,'diel-method\uncertainty-analysis\',site{kk}])
    load('MonteCarloResults')
    metab = diel_dtd_MC;

    clearvars finalQC diel_dtd diel_obs

    %====Calculate daily means for each site===============================
    params_dailyAvg = retime(params,'daily','mean');

    %====Create time table of all data=====================================
    % dt2 = dateshift(metab.daystart_dt,'start','day');
    % metab.daystart_dt = dt2;
    dailyAvg.(site{kk}) = synchronize(params_dailyAvg,wspd_dailyAvg,par_dailyAvg,metab);
    dailyAvg.(site{kk}) = removevars(dailyAvg.(site{kk}),{'deployment','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg',...
        'R_hourly_sd','P_hourly_sd','R_daily_sd','P_daily_sd','datetime_local','Tair','light_lux'});
    dailyAvg.(site{kk})(dailyAvg.(site{kk}).datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
    dailyAvg.(site{kk})(dailyAvg.(site{kk}).datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

    % Remove rows with anomalous GPP and ER values
    dailyAvg.(site{kk})((find(dailyAvg.(site{kk}).ER_avg > 0)),:) = [];
    dailyAvg.(site{kk})((find(dailyAvg.(site{kk}).GPP_avg < 0)),:) = [];

    % Remove rows with missing NEM values
    dailyAvg.(site{kk})(find(isnan(dailyAvg.(site{kk}).NEM_avg)),:) = [];

    %======================================================================
    % Find high wind speed events
    %======================================================================
    % Use daily wind speed
    dt_utc = dailyAvg.(site{kk}).datetime_utc;

    A = dailyAvg.(site{kk}).wspd;
    percentile = 99;
    th = prctile(A,percentile);   % Set threshold based on chosen percentile
    idx = find(diff(sign(dailyAvg.(site{kk}).wspd - th))); % Indices of all days where wspd crosses 'th'

    t_th = dt_utc(idx);
    t_thm = reshape(t_th, 2, []);       % Matrix Of Paired Events
    t_thm(2,:) = t_thm(2,:) + days(1);  % Add a day to end of each event

    % Remove second storm event for South, because event is cut off by missing data
    % if kk == 3
    %     t_thm(:,2) = [];
    % else
    %     % Do nothing
    % end

    eventDuration = diff(t_thm);

    % Find mean NEM during storm events
    start_stm = round(interp1(dailyAvg.(site{kk}).datetime_utc, 1:length(dailyAvg.(site{kk}).datetime_utc), t_thm(1,:)));
    end_stm = round(interp1(dailyAvg.(site{kk}).datetime_utc, 1:length(dailyAvg.(site{kk}).datetime_utc), t_thm(2,:)));
    for i = 1:length(t_thm)
        event_NEM(i) = mean(dailyAvg.(site{kk}).NEM_avg(start_stm(i):end_stm(i)));
    end

    % Find mean NEM during pre-storm period
    start_pre = start_stm - 6;
    end_pre = start_stm - 1;
    for i = 1:length(t_thm)
        pre_NEM(i) = mean(dailyAvg.(site{kk}).NEM_avg(start_pre(i):end_pre(i)));
    end

    % Find mean NEM during post-storm period
    start_post = end_stm + 1;
    end_post = end_stm + 6;
    for i = 1:length(t_thm)
        post_NEM(i) = mean(dailyAvg.(site{kk}).NEM_avg(start_post(i):end_post(i)));
    end

    pre = table(start_pre',end_pre','VariableNames',{'start','end'});
    event = table(start_stm',end_stm','VariableNames',{'start','end'});

    post = table(start_post',end_post','VariableNames',{'start','end'});
    events.(site{kk}) = struct('pre',pre,'event',event,'post',post,'eventDuration',days(eventDuration'));

    %======================================================================
    % Plots
    %======================================================================
    cd([rootpath,'figures\stats-analyses\storm events'])
    
    %----Daily wind, NEM, and >99%tiles marked-----------------------------
    fig(kk*2-1) = figure(kk*2-1);clf
    fig(kk*2-1).WindowState = 'maximized';

    yyaxis left
    plot(dailyAvg.(site{kk}).datetime_utc,dailyAvg.(site{kk}).wspd,'.-','DisplayName','Daily wind')
    hold on
    plot(dt_utc, ones(size(dt_utc))*th,':','Color',rgb('purple'),'DisplayName','99% tile (wind)')
    plot(t_thm, ones(size(t_thm))*th, '|-r','HandleVisibility','off')
    yline(th,'HandleVisibility','off')
    ylabel('Wind speed (m/s)')

    yyaxis right
    plot(dailyAvg.(site{kk}).datetime_utc,dailyAvg.(site{kk}).NEM_avg,'.-k','DisplayName',site{kk},'HandleVisibility','off')
    hold on
    for i = 1:width(t_thm)
        plot(dailyAvg.(site{kk}).datetime_utc(start_stm(i):end_stm(i)),dailyAvg.(site{kk}).NEM_avg(start_stm(i):end_stm(i)),'|-r','HandleVisibility','off')
        plot(dailyAvg.(site{kk}).datetime_utc(start_pre(i):end_pre(i)),dailyAvg.(site{kk}).NEM_avg(start_pre(i):end_pre(i)),'|-m','HandleVisibility','off')
        plot(dailyAvg.(site{kk}).datetime_utc(start_post(i):end_post(i)),dailyAvg.(site{kk}).NEM_avg(start_post(i):end_post(i)),'|-m','HandleVisibility','off')
    end
    legend('show','location','northwest')
    ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
    ax = gca;
    ax.YColor = 'k';

    xlim([min(dailyAvg.(site{kk}).datetime_utc) max(dailyAvg.(site{kk}).datetime_utc)])
    title(siteNames{kk})

    %======================================================================
    % Paired t-test
    %======================================================================
    x = pre_NEM;
    y = event_NEM;
    z = post_NEM;

    [hx,px] = ttest(x,y);
    [hz,pz] = ttest(z,y);

    % Save NEM means for pooled tests later
    pre_NEM_all.(site{kk}) = x;
    event_NEM_all.(site{kk}) = y;
    post_NEM_all.(site{kk}) = z;

    % Make table for box plot
    NEM_tbl.(site{kk}) = array2table(zeros(length(x)*3,2),'VariableNames',{'bin','NEM'});
    NEM_tbl.(site{kk}).bin = [repelem("Pre",length(x)),repelem("Event",length(x)),repelem("Post",length(x))]';
    NEM_tbl.(site{kk}).NEM = [x'; y'; z'];

    binOrder = ["Pre","Event","Post"];
    valueset = ["Pre","Event","Post"];
    binNames = categorical(NEM_tbl.(site{kk}).bin,valueset,binOrder);
    
    fig(kk*2) = figure(kk*2);clf    
    ax = axes;
    hold(ax)
    for i=1:3
        if isequal(kk,2)
            boxchart(binNames(i*6-5:i*6),NEM_tbl.(site{kk}).NEM(i*6-5:i*6),'BoxFaceColor',colors(i,:),...
                'BoxEdgeColor','k','BoxWidth',0.6,'MarkerColor',colors(i,:))
        elseif isequal(kk,1) || isequal(kk,3)
            boxchart(binNames(i*5-4:i*5),NEM_tbl.(site{kk}).NEM(i*5-4:i*5),'BoxFaceColor',colors(i,:),...
                'BoxEdgeColor','k','BoxWidth',0.6,'MarkerColor',colors(i,:))
        end
    end

    ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
    title(siteNames{kk})
    annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'String', {['p = ',num2str(px,1),' (pre)'],['p = ',num2str(pz,1),' (post)']},'LineWidth',1)

    disp('Press enter to continue to next site')
    pause
end
%%
%==========================================================================
% Paired t-test for all sites
%==========================================================================
x = [pre_NEM_all.north pre_NEM_all.gull pre_NEM_all.south];
y = [event_NEM_all.north event_NEM_all.gull event_NEM_all.south];
z = [post_NEM_all.north post_NEM_all.gull post_NEM_all.south];

[hx,px] = ttest(x,y);
[hz,pz] = ttest(z,y);

% Make table for box plot
NEM_tbl_all = array2table(zeros(length(x)*3,2),'VariableNames',{'bin','NEM'});
NEM_tbl_all.bin = [repelem("Pre",length(x)),repelem("Event",length(x)),repelem("Post",length(x))]';
NEM_tbl_all.NEM = [x'; y'; z'];

binOrder = ["Pre","Event","Post"];
valueset = ["Pre","Event","Post"];
binNames = categorical(NEM_tbl_all.bin,valueset,binOrder);

figure(7),clf;
ax = axes;
hold(ax)
for i=1:3
    boxchart(binNames(i*18-17:i*18),NEM_tbl_all.NEM(i*18-17:i*18),'BoxFaceColor',colors(i,:),...
        'BoxEdgeColor','k','BoxWidth',0.6,'MarkerColor',colors(i,:))
end
ylabel('NEM (mmol O_2 m^{-2} d^{-1})','FontSize',14)

set(gca,'box','on','LineWidth',1)

annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'String', {['p = ',num2str(px,1),' (pre)'],['p = ',num2str(pz,1),' (post)']},'LineWidth',1)

% title('All Sites')

%%
%==========================================================================
% Case study: TS Ian
%==========================================================================
figure(8),clf
yyaxis left
hold all
for kk = 1:length(site)
    h = plot(dailyAvg.(site{kk}).datetime_utc,dailyAvg.(site{kk}).NEM_avg,'k','linewidth',2,'DisplayName',siteNames{kk});
    if kk == 1
        h.LineStyle = '-';
        h.Color = north_clr;
    elseif kk == 2
        h.LineStyle = '-.';
        h.Color = gull_clr;
    elseif kk == 3
        h.LineStyle = ':';
        h.Color = south_clr;
    end
end
% Shade regions
xregion(dailyAvg.gull.datetime_utc(events.gull.pre.start),dailyAvg.gull.datetime_utc(events.gull.pre.end),'FaceColor',pre_clr,'HandleVisibility','off')
xregion(dailyAvg.gull.datetime_utc(events.gull.event.start),dailyAvg.gull.datetime_utc(events.gull.event.end),'FaceColor',event_clr,'HandleVisibility','off')
xregion(dailyAvg.gull.datetime_utc(events.gull.post.start),dailyAvg.gull.datetime_utc(events.gull.post.end),'FaceColor',post_clr,'HandleVisibility','off')
ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
legend('show','location','best')
ax = gca;
ax.YColor = 'k';

yyaxis right
plot(dailyAvg.north.datetime_utc,dailyAvg.north.wspd,'Color',rgb('DarkGray'),'linewidth',2,'HandleVisibility','off')
ylabel('Wind speed (m/s)')
ylim([1 15])
ax = gca;
ax.YColor = rgb('DarkGray');

set(gca,'box','on','LineWidth',1)

xlim([dailyAvg.gull.datetime_utc(events.gull.pre.start(2) - 2) dailyAvg.gull.datetime_utc(events.gull.post.end(2) + 2)])

%% Plot time series indicating storms
fig9=figure(9);clf
fig9.WindowState = 'maximized';
hold all
for kk = 1:length(site)
    h = plot(dailyAvg.(site{kk}).datetime_utc,dailyAvg.(site{kk}).NEM_avg,'k','linewidth',1.5,'DisplayName',siteNames{kk});
    if kk == 1
        h.LineStyle = '-';
        h.Color = north_clr;
    elseif kk == 2
        h.LineStyle = '-.';
        h.Color = gull_clr;
    elseif kk == 3
        h.LineStyle = ':';
        h.Color = south_clr;
    end
end
% Shade regions
xregion(dailyAvg.gull.datetime_utc(events.gull.event.start),dailyAvg.gull.datetime_utc(events.gull.event.end),'FaceColor',event_clr,'HandleVisibility','off')
ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
legend('show')
ylim([-500 400])

yyaxis right
plot(dailyAvg.north.datetime_utc,dailyAvg.north.wspd,'Color',rgb('DarkGray'),'linewidth',1.5,'HandleVisibility','off')
ylabel('Wind speed (m/s)')
ylim([-15 15])
ax = gca;
ax.YColor = rgb('DarkGray');

set(gca,'box','on','LineWidth',1)


%% Old: TS Ian
% figure(5),clf
% yyaxis left
% plot(dailyAvg.north.datetime_utc,dailyAvg.north.wspd,'linewidth',3,'HandleVisibility','off')
% ylabel('Wind speed (m/s)')
% ylim([1 15])
% 
% yyaxis right
% hold all
% for kk = 1:length(site)
%     h = plot(dailyAvg.(site{kk}).datetime_utc,dailyAvg.(site{kk}).NEM_avg,'k','linewidth',3,'DisplayName',site{kk});
%     if kk == 1
%         h.LineStyle = '-';
%     elseif kk == 2
%         h.LineStyle = ':';
%     elseif kk == 3
%         h.LineStyle = '-.';
%     end
% end
% xline(dailyAvg.gull.datetime_utc(events.gull.event.start),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
% xline(dailyAvg.gull.datetime_utc(events.gull.event.end),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
% xline(dailyAvg.gull.datetime_utc(events.gull.pre.start),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
% xline(dailyAvg.gull.datetime_utc(events.gull.pre.end),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
% xline(dailyAvg.gull.datetime_utc(events.gull.post.start),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
% xline(dailyAvg.gull.datetime_utc(events.gull.post.end),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
% ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
% legend('show')
% ax = gca;
% ax.YColor = 'k';
% xlim([dailyAvg.gull.datetime_utc(events.gull.pre.start(2) - 2) dailyAvg.gull.datetime_utc(events.gull.post.end(2) + 2)])
% % ylim([-80 40])
% annotation('textbox',[.3 .1 .1 .1],'String','PRE','FontSize',18,'LineStyle','none')
% annotation('arrow',[.30 .22],[.18 .18])
% annotation('arrow',[.34 .428],[.18 .18])
% annotation('textbox',[.49 .1 .1 .1],'String','EVENT','FontSize',18,'LineStyle','none')
% annotation('arrow',[.49 .479],[.18 .18])
% annotation('arrow',[.548 .559],[.18 .18])
% annotation('textbox',[.69 .1 .1 .1],'String','POST','FontSize',18,'LineStyle','none')
% annotation('arrow',[.69 .608],[.18 .18])
% annotation('arrow',[.74 .815],[.18 .18])

%% Extra plots
% %==========================================================================
% % Plot wind direction
% %==========================================================================
% x = datenum(wspd_dailyAvg.datetime);            % Need to convert datetime values to use quiver
% y = zeros(length(wspd_dailyAvg.wspd),1);        % Station remains stationary
% 
% x_storm = datenum(t_thm);
% 
% clear idx
% for kk = 1:length(x_storm)
%     idx(1,kk) = find(x == x_storm(1,kk));
%     idx(2,kk) = find(x == x_storm(2,kk));
% end
% 
% figure(5);clf
% quiver(x,y,wspd_dailyAvg.u10,wspd_dailyAvg.v10,'linewidth',2)
% hold on
% quiver(x(idx),y(idx),wspd_dailyAvg.u10(idx),wspd_dailyAvg.v10(idx),'r')
% datetick('x','dd-mmm-yyyy')
% xlim([min(x) max(x)])
% 
% %==========================================================================
% % Histograms of wind speed and precip data
% %==========================================================================
% figure,clf
% histogram(era5Dat.wspd,'DisplayName','Hourly')
% hold on
% histogram(wspd_6h,'DisplayName','6-h moving mean')
% histogram(dailyAvg.(site{kk}).wspd,'DisplayName','Daily')
% xlabel('Wind speed (m/s)')
% ylabel('Count')
% legend('show')
% 
% figure,clf
% histogram(precip1.PRCP,'DisplayName','Ventnor')
% hold on
% histogram(precip2.PRCP,'DisplayName','Dennis Twp')
% xlabel('Precip')
% ylabel('Count')
% legend('show')
% 
% %==========================================================================
% % Compare precip from two stations
% %==========================================================================
% precip1(find(isnan(precip1.PRCP)),:) = [];
% precip2(find(isnan(precip2.PRCP)),:) = [];
% precip1_TT = timetable(precip1.DATE,precip1.PRCP,'VariableNames',"PRCP");
% precip2_TT = timetable(precip2.DATE,precip2.PRCP,'VariableNames',"PRCP");
% precip12_TT = synchronize(precip1_TT,precip2_TT,'intersect');
% 
% % Create the model
% x = precip12_TT.PRCP_precip1_TT;
% y = precip12_TT.PRCP_precip2_TT;
% mdl = fitlm(x,y,'y~x1-1');  % Force intercept through zero; see Wilkinson notation
% 
% % Create equation string for plot
% eqn = ['y = ',num2str(mdl.Coefficients.Estimate,3),'x'];
% R2 = num2str(mdl.Rsquared.Ordinary,2);
% 
% figure,clf
% h = plot(mdl,'marker','.','markersize',20);
% hold on
% plot([min(min([x,y])) max(max([x,y]))], [min(min([x,y])) max(max([x,y]))],'--k')
% legend('',[eqn,newline,'R^2 = ',R2],'95% confidence interval','','1:1 line','Location','southeast')
% xlabel('Ventnor precip (mm)')
% ylabel('Dennis Twp precip (mm)')
% daspect([1 1 1])
% title('')
% 
% %==========================================================================
% %----Compare hourly, 6-h moving mean, and daily wind speed-----------------
% %==========================================================================
% wspd_6h = movmean(era5Dat.wspd,6);
% figure(3),clf
% yyaxis left
% plot(era5Dat.datetime,era5Dat.wspd,':','Color',rgb('darkblue'),'DisplayName','Hourly')
% hold on
% plot(era5Dat.datetime,wspd_6h,'-','Color',"#4DBEEE",'DisplayName','6-h moving mean')
% plot(dt_utc,dailyAvg.(site{kk}).wspd,'-','Color',"#0072BD",'DisplayName','Daily')
% plot(dt_utc, ones(size(dt_utc))*th,':','Color',rgb('purple'),'DisplayName','99% tile (wind)')
% plot(t_thm, ones(size(t_thm))*th, '|-r','HandleVisibility','off')
% plot(dt_utc, dailyAvg.(site{kk}).NEM,'-k','DisplayName','NEM')
% ylabel('Wind speed (m/s); NEM (mmol O_2 m^{-2} d^{-1})')
% ax = gca;
% ax.YColor = 'k';
% legend('show')
% 
% th_precip = prctile(precip2.PRCP,99);
% yyaxis right
% plot(precip1.DATE,precip1.PRCP,':','Color','#A2142F','DisplayName','Ventnor')
% hold on
% plot(precip2.DATE,precip2.PRCP,':','Color','#D95319','DisplayName','Dennis Twp')
% plot(precip1.DATE,ones(size(precip1.DATE))*th_precip, ':','Color','#EDB120','DisplayName','99% tile (precip)')
% ylabel('Precipitation (mm)')
% ax = gca;
% ax.YColor = 'k';
% 
% xlim([min(dailyAvg.(site{kk}).datetime_utc) max(dailyAvg.(site{kk}).datetime_utc)])
% 
% title(site{kk})