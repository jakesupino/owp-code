%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stormEvents.m
% This script picks out major storms during the 3-year OWP measurement period and 
% evaluates whether NEM values differ significantly during a storm event
% and the pre- and post-storm periods, using the "original" NEM calculated from "dielAnalysis_final.m".
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 8/21/2024
% Last updated: 8/28/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import windspeed and PAR data=====================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

cd([rootpath,'physical-data\final-dataset'])
load('par.mat')

%====Import precip data from NOAA======================================
cd([rootpath,'physical-data\NOAA'])
precip1 = readtable('ventnor.csv');
precip1.DATE = datetime(precip1.DATE,'TimeZone','America/New_York');

precip2 = readtable('dennis.csv');
precip2.DATE = datetime(precip2.DATE,'TimeZone','America/New_York');

%====Calculate daily means for each site===============================
wspd_dailyAvg = retime(era5Dat,'daily','mean');
par_dailyAvg = retime(parDat,'daily','mean');

%==========================================================================
% Run analyses on individual sites and save results in "NEM_tbl_all" for
% pooled t-test (after for for loop)
%==========================================================================
siteNames = {'North','Gull','South'};

for kk = 1:length(siteNames)
    clearvars start_pre start_stm start_post end_pre end_stm end_post pre_NEM event_NEM post_NEM

    %====Import final QC'd data & diel analysis results====================
    cd([rootpath,'open-water-platform-data\',siteNames{kk},'\cleaned\final-qc'])
    load([siteNames{kk},'-cleaned.mat']);
    cd([rootpath,'diel-method\matlab-results\final-qc\',siteNames{kk}])
    % cd([rootpath,'diel-method\matlab-results\final-qc\',siteNames{kk},'\wanninkhof'])   % HERE
    params = finalQC;
    load('diel_res.mat')
    metab = diel_dtd;

    clearvars finalQC diel_dtd diel_obs

    %====Calculate daily means for each site===============================
    params_dailyAvg = retime(params,'daily','mean');

    %====Create time table of all data=====================================
    % dt2 = dateshift(metab.daystart_dt,'start','day');
    % metab.daystart_dt = dt2;
    dailyAvg.(siteNames{kk}) = synchronize(params_dailyAvg,wspd_dailyAvg,par_dailyAvg,metab);
    dailyAvg.(siteNames{kk}) = removevars(dailyAvg.(siteNames{kk}),{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
    dailyAvg.(siteNames{kk})(dailyAvg.(siteNames{kk}).datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
    dailyAvg.(siteNames{kk})(dailyAvg.(siteNames{kk}).datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

    % Remove rows with anomalous GPP and ER values
    dailyAvg.(siteNames{kk})((find(dailyAvg.(siteNames{kk}).ER > 0)),:) = [];
    dailyAvg.(siteNames{kk})((find(dailyAvg.(siteNames{kk}).GPP < 0)),:) = [];

    % Remove rows with missing NEM values
    dailyAvg.(siteNames{kk})(find(isnan(dailyAvg.(siteNames{kk}).NEM)),:) = [];

    %======================================================================
    % Find high wind speed events
    %======================================================================
    % Use daily wind speed
    time = dailyAvg.(siteNames{kk}).datetime_utc;

    A = dailyAvg.(siteNames{kk}).wspd;
    th = prctile(A,99);   % Set threshold based on chosen percentile
    idx = find(diff(sign(dailyAvg.(siteNames{kk}).wspd - th))); % Indices of all days where wspd crosses 'th'

    t_th = time(idx);

    % figure,clf
    % plot(dailyAvg.(siteNames{kk}).datetime_utc,dailyAvg.(siteNames{kk}).wspd)
    % hold on
    % plot(dailyAvg.(siteNames{kk}).datetime_utc(idx),dailyAvg.(siteNames{kk}).wspd(idx),'o')
    % plot(dailyAvg.(siteNames{kk}).datetime_utc,dailyAvg.(siteNames{kk}).NEM,'k')
    % plot(dailyAvg.(siteNames{kk}).datetime_utc(idx),dailyAvg.(siteNames{kk}).NEM(idx),'o')
    % yline(th,'--')
    % ylabel('wspd')
    % title(siteNames{kk})

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
    start_stm = round(interp1(dailyAvg.(siteNames{kk}).datetime_utc, 1:length(dailyAvg.(siteNames{kk}).datetime_utc), t_thm(1,:)));
    end_stm = round(interp1(dailyAvg.(siteNames{kk}).datetime_utc, 1:length(dailyAvg.(siteNames{kk}).datetime_utc), t_thm(2,:)));
    for i = 1:length(t_thm)
        event_NEM(i) = mean(dailyAvg.(siteNames{kk}).NEM(start_stm(i):end_stm(i)));
    end

    % Find mean NEM during pre-storm period
    start_pre = start_stm - 6;
    end_pre = start_stm - 1;
    for i = 1:length(t_thm)
        pre_NEM(i) = mean(dailyAvg.(siteNames{kk}).NEM(start_pre(i):end_pre(i)));
    end

    % Find mean NEM during post-storm period
    start_post = end_stm + 1;
    end_post = end_stm + 6;
    for i = 1:length(t_thm)
        post_NEM(i) = mean(dailyAvg.(siteNames{kk}).NEM(start_post(i):end_post(i)));
    end

    pre = table(start_pre',end_pre','VariableNames',{'start','end'});
    event = table(start_stm',end_stm','VariableNames',{'start','end'});

    post = table(start_post',end_post','VariableNames',{'start','end'});
    events.(siteNames{kk}) = struct('pre',pre,'event',event,'post',post,'eventDuration',days(eventDuration'));

    %======================================================================
    % Plots
    %======================================================================
    cd([rootpath,'figures\stats-analyses\storm events'])
    %----Daily wind, NEM, and >99%tiles marked-----------------------------
    figure(1),clf
    yyaxis left
    plot(dailyAvg.(siteNames{kk}).datetime_utc,dailyAvg.(siteNames{kk}).wspd,'.-','DisplayName','Daily wind')
    hold on
    plot(time, ones(size(time))*th,':','Color',rgb('purple'),'DisplayName','99% tile (wind)')
    plot(t_thm, ones(size(t_thm))*th, '|-r','HandleVisibility','off')
    yline(th,'HandleVisibility','off')
    ylabel('Wind speed (m/s)')

    yyaxis right
    plot(dailyAvg.(siteNames{kk}).datetime_utc,dailyAvg.(siteNames{kk}).NEM,'.-k','DisplayName',siteNames{kk})
    hold on
    for i = 1:width(t_thm)
        plot(dailyAvg.(siteNames{kk}).datetime_utc(start_stm(i):end_stm(i)),dailyAvg.(siteNames{kk}).NEM(start_stm(i):end_stm(i)),'|-r','HandleVisibility','off')
        plot(dailyAvg.(siteNames{kk}).datetime_utc(start_pre(i):end_pre(i)),dailyAvg.(siteNames{kk}).NEM(start_pre(i):end_pre(i)),'|-m','HandleVisibility','off')
        plot(dailyAvg.(siteNames{kk}).datetime_utc(start_post(i):end_post(i)),dailyAvg.(siteNames{kk}).NEM(start_post(i):end_post(i)),'|-m','HandleVisibility','off')
    end
    legend('show','location','northwest')
    ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
    ax = gca;
    ax.YColor = 'k';

    xlim([min(dailyAvg.(siteNames{kk}).datetime_utc) max(dailyAvg.(siteNames{kk}).datetime_utc)])
    title(siteNames{kk})

    %----Daily wind, NEM, temperature, and >99%tiles marked----------------
    figure(2),clf
    yyaxis left
    plot(dailyAvg.(siteNames{kk}).datetime_utc,dailyAvg.(siteNames{kk}).wspd,'.-','DisplayName','Daily wind')
    hold on
    plot(time, ones(size(time))*th,':','Color',rgb('purple'),'DisplayName','99% tile (wind)')
    plot(t_thm, ones(size(t_thm))*th, '|-r','HandleVisibility','off')
    yline(th,'HandleVisibility','off')
    plot(dailyAvg.(siteNames{kk}).datetime_utc,dailyAvg.(siteNames{kk}).NEM,'.-k','DisplayName',siteNames{kk},'DisplayName','NEM')
    for i = 1:width(t_thm)
        plot(dailyAvg.(siteNames{kk}).datetime_utc(start_stm(i):end_stm(i)),dailyAvg.(siteNames{kk}).NEM(start_stm(i):end_stm(i)),'|-r','HandleVisibility','off')
        plot(dailyAvg.(siteNames{kk}).datetime_utc(start_pre(i):end_pre(i)),dailyAvg.(siteNames{kk}).NEM(start_pre(i):end_pre(i)),'|-m','HandleVisibility','off')
        plot(dailyAvg.(siteNames{kk}).datetime_utc(start_post(i):end_post(i)),dailyAvg.(siteNames{kk}).NEM(start_post(i):end_post(i)),'|-m','HandleVisibility','off')
    end
    ylabel('Wind speed (m/s) and NEM (mmol O_2 m^{-2} d^{-1})')

    yyaxis right
    plot(dailyAvg.(siteNames{kk}).datetime_utc,dailyAvg.(siteNames{kk}).temperature,':','HandleVisibility','off')
    ylabel('Temperature (^oC)')

    legend('show','location','southwest')

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
    pre_NEM_all.(siteNames{kk}) = x;
    event_NEM_all.(siteNames{kk}) = y;
    post_NEM_all.(siteNames{kk}) = z;

    % Make table for box plot
    NEM_tbl.(siteNames{kk}) = array2table(zeros(length(x)*3,2),'VariableNames',{'bin','NEM'});
    NEM_tbl.(siteNames{kk}).bin = [repelem("Pre",length(x)),repelem("Event",length(x)),repelem("Post",length(x))]';
    NEM_tbl.(siteNames{kk}).NEM = [x'; y'; z'];

    figure(3),clf
    boxplot(NEM_tbl.(siteNames{kk}).NEM,NEM_tbl.(siteNames{kk}).bin,'Colors','k')
    ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
    title(siteNames{kk})
    annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'String', {['p = ',num2str(px,1),' (pre)'],['p = ',num2str(pz,1),' (post)']})

    disp('Press enter to continue to next site')
    pause

end

%==========================================================================
% Paired t-test for all sites
%==========================================================================
x = [pre_NEM_all.North pre_NEM_all.Gull pre_NEM_all.South];
y = [event_NEM_all.North event_NEM_all.Gull event_NEM_all.South];
z = [post_NEM_all.North post_NEM_all.Gull post_NEM_all.South];

[hx,px] = ttest(x,y);
[hz,pz] = ttest(z,y);

% Make table for box plot
NEM_tbl_all = array2table(zeros(length(x)*3,2),'VariableNames',{'bin','NEM'});
NEM_tbl_all.bin = [repelem("Pre",length(x)),repelem("Event",length(x)),repelem("Post",length(x))]';
NEM_tbl_all.NEM = [x'; y'; z'];

figure(4),clf
boxplot(NEM_tbl_all.NEM,NEM_tbl_all.bin,'Colors','k')
ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
title('All Sites')
annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'String', {['p = ',num2str(px,1),' (pre)'],['p = ',num2str(pz,1),' (post)']})

%%
%==========================================================================
% Case study: TS Ian
%==========================================================================
figure(5),clf
yyaxis left
plot(dailyAvg.North.datetime_utc,dailyAvg.North.wspd,'linewidth',3,'HandleVisibility','off')
ylabel('Wind speed (m/s)')
ylim([1 15])

yyaxis right
hold all
for kk = 1:length(siteNames)
    h = plot(dailyAvg.(siteNames{kk}).datetime_utc,dailyAvg.(siteNames{kk}).NEM,'k','linewidth',3,'DisplayName',siteNames{kk});
    if kk == 1
        h.LineStyle = '-';
    elseif kk == 2
        h.LineStyle = ':';
    elseif kk == 3
        h.LineStyle = '-.';
    end
end
xline(dailyAvg.Gull.datetime_utc(events.Gull.event.start),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
xline(dailyAvg.Gull.datetime_utc(events.Gull.event.end),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
xline(dailyAvg.Gull.datetime_utc(events.Gull.pre.start),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
xline(dailyAvg.Gull.datetime_utc(events.Gull.pre.end),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
xline(dailyAvg.Gull.datetime_utc(events.Gull.post.start),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
xline(dailyAvg.Gull.datetime_utc(events.Gull.post.end),'--','linewidth',1,'Color',rgb('grey'),'HandleVisibility','off')
ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
legend('show')
ax = gca;
ax.YColor = 'k';
xlim([dailyAvg.Gull.datetime_utc(events.Gull.pre.start(2) - 2) dailyAvg.Gull.datetime_utc(events.Gull.post.end(2) + 2)])
% ylim([-80 40])
annotation('textbox',[.3 .1 .1 .1],'String','PRE','FontSize',18,'LineStyle','none')
annotation('arrow',[.30 .22],[.18 .18])
annotation('arrow',[.34 .428],[.18 .18])
annotation('textbox',[.49 .1 .1 .1],'String','EVENT','FontSize',18,'LineStyle','none')
annotation('arrow',[.49 .479],[.18 .18])
annotation('arrow',[.548 .559],[.18 .18])
annotation('textbox',[.69 .1 .1 .1],'String','POST','FontSize',18,'LineStyle','none')
annotation('arrow',[.69 .608],[.18 .18])
annotation('arrow',[.74 .815],[.18 .18])

%% Extra plots
%==========================================================================
% Plot wind direction
%==========================================================================
x = datenum(wspd_dailyAvg.datetime);            % Need to convert datetime values to use quiver
y = zeros(length(wspd_dailyAvg.wspd),1);        % Station remains stationary

x_storm = datenum(t_thm);

clear idx
for i = 1:length(x_storm)
    idx(1,i) = find(x == x_storm(1,i));
    idx(2,i) = find(x == x_storm(2,i));
end

figure(5);clf
quiver(x,y,wspd_dailyAvg.u10,wspd_dailyAvg.v10,'linewidth',2)
hold on
quiver(x(idx),y(idx),wspd_dailyAvg.u10(idx),wspd_dailyAvg.v10(idx),'r')
datetick('x','dd-mmm-yyyy')
xlim([min(x) max(x)])

%==========================================================================
% Histograms of wind speed and precip data
%==========================================================================
figure,clf
histogram(era5Dat.wspd,'DisplayName','Hourly')
hold on
histogram(wspd_6h,'DisplayName','6-h moving mean')
histogram(dailyAvg.(siteNames{kk}).wspd,'DisplayName','Daily')
xlabel('Wind speed (m/s)')
ylabel('Count')
legend('show')

figure,clf
histogram(precip1.PRCP,'DisplayName','Ventnor')
hold on
histogram(precip2.PRCP,'DisplayName','Dennis Twp')
xlabel('Precip')
ylabel('Count')
legend('show')

%==========================================================================
% Compare precip from two stations
%==========================================================================
precip1(find(isnan(precip1.PRCP)),:) = [];
precip2(find(isnan(precip2.PRCP)),:) = [];
precip1_TT = timetable(precip1.DATE,precip1.PRCP,'VariableNames',"PRCP");
precip2_TT = timetable(precip2.DATE,precip2.PRCP,'VariableNames',"PRCP");
precip12_TT = synchronize(precip1_TT,precip2_TT,'intersect');

% Create the model
x = precip12_TT.PRCP_precip1_TT;
y = precip12_TT.PRCP_precip2_TT;
mdl = fitlm(x,y,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl.Rsquared.Ordinary,2);

figure,clf
h = plot(mdl,'marker','.','markersize',20);
hold on
plot([min(min([x,y])) max(max([x,y]))], [min(min([x,y])) max(max([x,y]))],'--k')
legend('',[eqn,newline,'R^2 = ',R2],'95% confidence interval','','1:1 line','Location','southeast')
xlabel('Ventnor precip (mm)')
ylabel('Dennis Twp precip (mm)')
daspect([1 1 1])
title('')

%==========================================================================
%----Compare hourly, 6-h moving mean, and daily wind speed-----------------
%==========================================================================
wspd_6h = movmean(era5Dat.wspd,6);
figure(3),clf
yyaxis left
plot(era5Dat.datetime,era5Dat.wspd,':','Color',rgb('darkblue'),'DisplayName','Hourly')
hold on
plot(era5Dat.datetime,wspd_6h,'-','Color',"#4DBEEE",'DisplayName','6-h moving mean')
plot(time,dailyAvg.(siteNames{kk}).wspd,'-','Color',"#0072BD",'DisplayName','Daily')
plot(time, ones(size(time))*th,':','Color',rgb('purple'),'DisplayName','99% tile (wind)')
plot(t_thm, ones(size(t_thm))*th, '|-r','HandleVisibility','off')
plot(time, dailyAvg.(siteNames{kk}).NEM,'-k','DisplayName','NEM')
ylabel('Wind speed (m/s); NEM (mmol O_2 m^{-2} d^{-1})')
ax = gca;
ax.YColor = 'k';
legend('show')

th_precip = prctile(precip2.PRCP,99);
yyaxis right
plot(precip1.DATE,precip1.PRCP,':','Color','#A2142F','DisplayName','Ventnor')
hold on
plot(precip2.DATE,precip2.PRCP,':','Color','#D95319','DisplayName','Dennis Twp')
plot(precip1.DATE,ones(size(precip1.DATE))*th_precip, ':','Color','#EDB120','DisplayName','99% tile (precip)')
ylabel('Precipitation (mm)')
ax = gca;
ax.YColor = 'k';

xlim([min(dailyAvg.(siteNames{kk}).datetime_utc) max(dailyAvg.(siteNames{kk}).datetime_utc)])

title(siteNames{kk})

%%
%==========================================================================
% Find high wind speed events -- METHOD 1
%==========================================================================

% % Find wind speeds above chosen percentile
% A = daily_wspd_mph;
% P90 = prctile(A,90);   % 90th percentile
% P99 = prctile(A,99);   % 99th percentile
% id90 = find(A > P90);
% id99 = find(A > P99);
%
% daily_high_wspd = table(dat_TT.datetime_utc(id99),daily_wspd_mph(id99),'VariableNames',{'date','wspd'});
%
% % Find mean wind speed of events that last more than one day
% id_consec = find(diff(id99) == 1);
% for i = 1:length(id_consec)
%     mean_consec_wspd(i) = mean([daily_high_wspd.wspd(id_consec(i)) daily_high_wspd.wspd(id_consec(i)+1)]);
% end
%
% % Create new high wspd table containing means of consecutive days
% high_wspd = array2table(zeros(height(daily_high_wspd),4),'VariableNames',{'start_date','end_date','wspd','num_days'});
% high_wspd.start_date = datetime(daily_high_wspd.date,'TimeZone','UTC');
% high_wspd.end_date = datetime(daily_high_wspd.date,'TimeZone','UTC');
%
% % Find the end date of events that are >1 day long
% is_consec = diff(id99) == 1;
% for i = 1:length(is_consec)
%     if is_consec(i) == 1
%         high_wspd.end_date(i) = high_wspd.start_date(i) + days(1);
%     end
% end
%
% % Add the wind speed values to the table
% high_wspd.wspd = daily_high_wspd.wspd;
% high_wspd.wspd(id_consec) = mean_consec_wspd;
% high_wspd(id_consec+1,:) = [];
%
% % Find the length (number of days) of each event
% high_wspd.num_days = days(high_wspd.end_date - high_wspd.start_date + 1);
%
% %==========================================================================
% % Find mean NEM during, pre, and post event
% %==========================================================================
% ind_match = interp1(dat_TT.datetime_utc,1:length(dat_TT.datetime_utc),daily_high_wspd.date);
%
% daily_event_NEM = table(dat_TT.datetime_utc(ind_match),dat_TT.NEM(ind_match),'VariableNames',{'date','NEM'});
%
% id_consec = find(diff(ind_match) == 1);
% for i = 1:length(id_consec)
%     mean_consec_NEM(i) = mean([daily_event_NEM.NEM(id_consec(i)) daily_event_NEM.NEM(id_consec(i)+1)]);
% end
%
% % New NEM table containing means of consecutive days
% event_NEM = array2table(zeros(height(daily_event_NEM),4),'VariableNames',{'start_date','end_date','NEM','num_days'});
% event_NEM.start_date = datetime(daily_event_NEM.date,'TimeZone','UTC');
% event_NEM.end_date = datetime(daily_event_NEM.date,'TimeZone','UTC');
%
% is_consec = diff(ind_match) == 1;
% for i = 1:length(is_consec)
%     if is_consec(i) == 1
%         event_NEM.end_date(i) = event_NEM.start_date(i) + days(1);
%     end
% end
%
% event_NEM.NEM = daily_event_NEM.NEM;
% event_NEM.NEM(id_consec) = mean_consec_NEM;
% event_NEM(id_consec+1,:) = [];
%
% event_NEM.num_days = days(event_NEM.end_date - event_NEM.start_date + 1);
%
% % Pre-event NEM means
% pre_NEM = array2table(zeros(height(event_NEM),4),'VariableNames',{'start_date','end_date','NEM','num_days'});
% pre_NEM.start_date = event_NEM.start_date - days(6);
% pre_NEM.end_date = event_NEM.start_date - days(1);
%
% ind_start = interp1(dat_TT.datetime_utc,1:length(dat_TT.datetime_utc),pre_NEM.start_date);
% ind_end = interp1(dat_TT.datetime_utc,1:length(dat_TT.datetime_utc),pre_NEM.end_date);
%
% for i = 1:length(ind_start)
%     pre_NEM.NEM(i) = mean([dat_TT.NEM(ind_start(i)) dat_TT.NEM(ind_end(i))]);
% end
%
% % Post-event NEM means
% post_NEM = array2table(zeros(height(event_NEM),4),'VariableNames',{'start_date','end_date','NEM','num_days'});
% post_NEM.start_date = event_NEM.start_date + days(1);
% post_NEM.end_date = event_NEM.start_date + days(6);
%
% ind_start = interp1(dat_TT.datetime_utc,1:length(dat_TT.datetime_utc),post_NEM.start_date);
% ind_end = interp1(dat_TT.datetime_utc,1:length(dat_TT.datetime_utc),post_NEM.end_date);
%
% for i = 1:length(ind_start)
%     post_NEM.NEM(i) = mean([dat_TT.NEM(ind_start(i)) dat_TT.NEM(ind_end(i))]);
% end
%
% % Plot the results
% figure(1),clf
% yyaxis left
% plot(dat_TT.datetime_utc,daily_wspd_mph,'DisplayName','Daily wind')
% hold on
% plot(high_wspd.start_date,high_wspd.wspd,'xm','DisplayName','>99 percentile')
% yline(P90)
% yline(P99)
% ylabel('Daily wind speed (mph)')
%
% yyaxis right
% plot(dat_TT.datetime_utc,dat_TT.NEM,'.-k','DisplayName',site)
% hold on
% plot(event_NEM.start_date,event_NEM.NEM,'xm','DisplayName','Event NEM')
% plot(pre_NEM.start_date,pre_NEM.NEM,'og','DisplayName','Pre NEM')
% plot(post_NEM.end_date,post_NEM.NEM,'or','DisplayName','Post NEM')
% legend('show')
% ylabel('NEM')
%
% % Delete any rows with a missing NEM value
% % pre_NEM(isnan(pre_NEM.NEM),:) = [];
% % event_NEM(isnan(event_NEM.NEM),:) = [];
% % post_NEM(isnan(post_NEM.NEM),:) = [];
%
% % event_NEM(4,:) = [];
% % pre_NEM(4,:) = [];
% % post_NEM(4,:) = [];
%
% % Make table for box plot
% NEM_tbl = array2table(zeros(length(pre_NEM.NEM)*3,2),'VariableNames',{'bin','NEM'});
% NEM_tbl.bin = [repelem("Pre",length(pre_NEM.NEM)),repelem("Event",length(pre_NEM.NEM)),repelem("Post",length(pre_NEM.NEM))]';
% NEM_tbl.NEM = [pre_NEM.NEM; event_NEM.NEM; post_NEM.NEM];
%
% figure(2),clf
% boxplot(NEM_tbl.NEM,NEM_tbl.bin,'Colors','k')
% ylabel('NEM (mmol O_2 m^{-2} d^{-1})')
% title(site)
%
% %==========================================================================
% % Paired t-test
% %==========================================================================
% x = pre_NEM.NEM;
% y = event_NEM.NEM;
% z = post_NEM.NEM;
%
% [h,p] = ttest(x,y)
%
% [h,p] = ttest(z,y)