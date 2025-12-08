%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dupCheck_d_S_T.m
% This script compares the QC'd results from the BC and ERDC sondes for
% salinity, depth, and temperature and produces the "best-guess" time
% series data for each of these parameters.
%
% Code that requires manual input is commented with "INPUTS".
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 4/10/2024
% Last updated: 12/05/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import cleaned BC and ERDC data=======================================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\movmed'])
load([site,'-bc-cleaned.mat']);
load([site,'-erdc-cleaned.mat'])
dat1 = sonde1_cleaned;
dat2 = sonde2_cleaned;
% Synchronize the BC and ERDC data to a common datetime vector
dat_syn = synchronize(dat1,dat2);
dt_utc = dat_syn.datetime_utc;

clearvars dat1 dat2 sonde1_cleaned sonde2_cleaned

%====Import lab S from discrete samples====================================
cd('G:\Shared drives\SMIIL\Shared Data')
discreteS = readtable('SMIIL DIC_TA Sample Salinities - OWP.xlsx');

varUnits = ["","","","","","","","","","psu","psu","degC","",""];
discreteS.Properties.VariableUnits = varUnits;

discreteS.Sampling_Time = datetime(discreteS.Sampling_Time,'ConvertFrom','datenum');
myDatetime = discreteS.Date_Sampled + timeofday(discreteS.Sampling_Time);
myDatetime.TimeZone = 'America/New_York';
datetime_utc = myDatetime;
datetime_utc.TimeZone = 'UTC';
discreteS.datetime_utc = datetime_utc;
discreteS = movevars(discreteS,'datetime_utc','Before','trip_ID');

% Find indices of samples taken from this site
ind_site = find(ismember(discreteS.Location_ID,site));

%====Find the indices of deployment changes================================
switch site
    case 'South'    % South BC Dep. 15 is missing, so need to manually add in dep number
        ind_d15 = find(dat_syn.deployment_dat2 == 15,1);
        ind_d16 = find(dat_syn.deployment_dat1 == 16,1);
        dat_syn.deployment_dat1(ind_d15:ind_d16) = 15;
end

dep1 = fillmissing(dat_syn.deployment_dat1,'previous'); % Replace NaNs in BC deployment column with previous deployment #
dat_syn.deployment_dat1 = dep1;
ind_dep = find(diff(dat_syn.deployment_dat1) > 0);

switch site
    case 'Gull'
        dep = table([1;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
    case 'North'
        dep = table([2;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
    case 'South'
        dep = table([1;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
end
dep.Properties.VariableNames = {'depNum','ind'};

switch site
    case 'Gull'
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
            'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
            'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
            'Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
    case 'North'
        label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8',...
            'Deployment 9','Deployment 10','Deployment 11','Deployment 12',...
            'Deployment 13','Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
    case 'South'
        label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
            'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
            'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
            'Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
end

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = dat_syn.datetime_utc(end);
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 6;
circlesize = 4;

%====SALINITY==============================================================
cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\salinity'])
% Check agreement between BC and ERDC sondes
S_bc = dat_syn.salinity_dat1;
S_erdc = dat_syn.salinity_dat2;
threshold = 4; % psu
ind_diverge = find(abs(S_bc - S_erdc) > threshold);

pc_exceed = height(ind_diverge)/height(dat_syn)*100; % Percentage of synchronized data that exceed threshold

txt = ['Differs by >',num2str(threshold),' psu (',num2str(pc_exceed,'%.1f'),'%)'];

window = 144;
mmed_bc = movmedian(S_bc,window,'omitmissing');
mmed_erdc = movmedian(S_erdc,window,'omitmissing');
mmed_all = [mmed_bc,mmed_erdc];
mmed_mean = mean(mmed_all,2);

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(dt_utc,S_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC');
hold on
plot(dt_utc,S_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
% plot(dt_utc(ind_diverge),depth_bc(ind_diverge),'o','MarkerSize',circlesize,'Color',rgb('chartreuse'),'DisplayName',txt)
plot(dt_utc,mmed_bc,'-','Color',rgb('darkred'),'DisplayName','BC Moving Median')
plot(dt_utc,mmed_erdc,'-','Color',rgb('darkblue'),'DisplayName','ERDC Moving Median')
plot(dt_utc,mmed_mean,'-','Color',rgb('magenta'),'DisplayName','Mean Moving Median')
plot(discreteS.datetime_utc(ind_site),discreteS.Lab_Salinity(ind_site),'o','Color',rgb('gold'),'MarkerSize',6,'LineWidth',2,'DisplayName','Lab S')
xline(dt_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Salinity (psu)')
legend('show','location','best')
title([site,' - After Initial QC and Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

switch site
    case 'Gull'
        % Funky deployments
        % Dep 6 -- define points where sensors do not track together
        kk = find(ind_diverge > 29223 & ind_diverge < dep.ind(5));
        ind_d6 = ind_diverge(kk);

        % Choose "better" deployment
        d1 = mmed_erdc(dep.ind(1):dep.ind(2));        % Dep 1
        d2 = mmed_mean(dep.ind(2)+1:dep.ind(3));      % Dep 2
        d5 = mmed_mean(dep.ind(3)+1:dep.ind(4));      % Dep 5
        d6 = [mmed_mean(dep.ind(4)+1:ind_d6(1)); mmed_erdc(ind_d6(1)+1:dep.ind(5))];  % Dep 6
        d7 = mmed_mean(dep.ind(5)+1:dep.ind(6));      % Dep 7
        d8 = [mmed_mean(dep.ind(6)+1:41331); mmed_erdc(41332:dep.ind(7))];        % Dep 8
        d9 = mmed_erdc(dep.ind(7)+1:dep.ind(8));      % Dep 9
        d10 = mmed_bc(dep.ind(8)+1:dep.ind(9));       % Dep 10
        d11 = mmed_mean(dep.ind(9)+1:dep.ind(10));    % Dep 11
        d12 = mmed_mean(dep.ind(10)+1:dep.ind(11));   % Dep 12
        d13 = [mmed_erdc(dep.ind(11)+1:81445); mmed_bc(81446:dep.ind(12))];   % Dep 13
        d14 = mmed_bc(dep.ind(12)+1:dep.ind(13));     % Dep 14
        d15 = mmed_bc(dep.ind(13)+1:dep.ind(14));     % Dep 15
        d16 = [mmed_erdc(dep.ind(14)+1:111691); mmed_bc(111692:dep.ind(15))];     % Dep 16
        d17 = mmed_mean(dep.ind(15)+1:dep.ind(16));   % Dep 17
        d18 = mmed_bc(dep.ind(16)+1:end);           % Dep 18

        S_bestguess = [d1;d2;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];

        % Calculate mean difference for Deployments 11 & 12
        diff_S = diff(mmed_all(dep.ind(9):dep.ind(11),:),1,2);
        two_sigma.S = mean(abs(diff_S),'omitmissing');

    case 'North'
        % Funky deployments
        % Dep 7 -- offset adjustment
        delta_d7 = mmed_bc(20745) - mmed_erdc(20786);

        % Dep 15 -- offset adjustment
        delta_d15 = mmed_erdc(76160) - mmed_erdc(76159);

        % Choose "better" deployment
        d2 = [mmed_erdc(dep.ind(1):dep.ind(2))];          % Dep 2
        d6 = mmed_bc(dep.ind(2)+1:dep.ind(3));          % Dep 6
        d7 = mmed_bc(dep.ind(3)+1:dep.ind(4)) - delta_d7; % Dep 7
        d8 = mmed_erdc(dep.ind(4)+1:dep.ind(5));        % Dep 8
        d9 = mmed_erdc(dep.ind(5)+1:dep.ind(6));        % Dep 9
        d10 = mmed_erdc(dep.ind(6)+1:dep.ind(7));       % Dep 10
        d11 = mmed_mean(dep.ind(7)+1:dep.ind(8));       % Dep 11
        d12 = mmed_erdc(dep.ind(8)+1:dep.ind(9));       % Dep 12
        d13 = mmed_erdc(dep.ind(9)+1:dep.ind(10));      % Dep 13
        d14 = mmed_erdc(dep.ind(10)+1:dep.ind(11));     % Dep 14
        d15 = mmed_erdc(dep.ind(11)+1:dep.ind(12)) - delta_d15; % Dep 15
        d16 = mmed_mean(dep.ind(12)+1:dep.ind(13));     % Dep 16
        d17 = mmed_mean(dep.ind(13)+1:dep.ind(14));     % Dep 17
        d18 = mmed_erdc(dep.ind(14)+1:end);             % Dep 18

        S_bestguess = [d2;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];

        % Remove "obvious" outliers
        S_bestguess(4351:dep.ind(2)) = NaN;
        S_bestguess(11090:15492) = NaN;
        S_bestguess(20785) = NaN;
        S_bestguess(80182:dep.ind(12)) = NaN;
        S_bestguess(84381:84385) = NaN;

        % Calculate mean difference for Deployments 16 & 17
        diff_S = diff(mmed_all(dep.ind(12):dep.ind(14),:),1,2);
        two_sigma.S = mean(abs(diff_S),'omitmissing');

    case 'South'
        % Funky deployments
        % Dep 6 -- offset adjustment
        ind = dep.ind(5);
        delta_d6 = mmed_bc(ind+1) - mmed_bc(ind-1);

        % Dep 7 -- offset adjustment
        ind = dep.ind(7);
        delta_d7 = mmed_erdc(ind-11) - mmed_erdc(ind);

        % Dep 15 -- offset adjustment
        delta_d15 = S_erdc(98253) - S_erdc(98251);

        % Choose "better" deployment
        d1 = mmed_erdc(dep.ind(1):dep.ind(2));        % Dep 1
        d2 = mmed_mean(dep.ind(2)+1:dep.ind(3));      % Dep 2
        d4 = mmed_bc(dep.ind(3)+1:dep.ind(4));        % Dep 4
        d5 = mmed_bc(dep.ind(4)+1:dep.ind(5));        % Dep 5
        d6 = mmed_bc(dep.ind(5)+1:dep.ind(6)) - delta_d6;        % Dep 6
        d7 = mmed_erdc(dep.ind(6)+1:dep.ind(7)) - delta_d7;      % Dep 7
        d8 = mmed_erdc(dep.ind(7)+1:dep.ind(8));      % Dep 8
        d9 = mmed_mean(dep.ind(8)+1:dep.ind(9));      % Dep 9
        d10 = mmed_bc(dep.ind(9)+1:dep.ind(10));      % Dep 10
        d11 = mmed_erdc(dep.ind(10)+1:dep.ind(11));   % Dep 11
        d12 = mmed_mean(dep.ind(11)+1:dep.ind(12));   % Dep 12
        d13 = mmed_erdc(dep.ind(12)+1:dep.ind(13));   % Dep 13
        d14 = mmed_mean(dep.ind(13)+1:dep.ind(14));   % Dep 14
        d15 = mmed_erdc(dep.ind(14)+1:dep.ind(15)) + delta_d15;   % Dep 15
        d16 = mmed_mean(dep.ind(15)+1:dep.ind(16));   % Dep 16
        d17 = mmed_erdc(dep.ind(16)+1:dep.ind(17));   % Dep 17
        d18 = mmed_mean(dep.ind(17)+1:end);           % Dep 18
        
        S_bestguess = [d1;d2;d4;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];
        
        % Remove "obvious" outliers
        S_bestguess(5896:5898) = NaN;
        S_bestguess(11043:dep.ind(3)) = NaN;
        S_bestguess(11928) = NaN;
        S_bestguess(29482) = NaN;
        S_bestguess(33012) = NaN;
        S_bestguess(38194:38204) = NaN;
        S_bestguess(49118:49149) = NaN;
        S_bestguess(77744:77797) = NaN;
        S_bestguess(80846:80852) = NaN;

        % Calculate mean absolute difference for Deployments 9 & 12
        diff_S = diff(mmed_all([dep.ind(8):dep.ind(9),dep.ind(11):dep.ind(12)],:),1,2);
        two_sigma.S = mean(abs(diff_S),'omitmissing');
end

clearvars d1 d2 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 d17 d18

fig2 = figure(2);clf
fig2.WindowState = 'maximized';
plot(dt_utc,S_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC');
hold on
plot(dt_utc,S_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
plot(dt_utc,mmed_mean,'.','Color',rgb('magenta'),'DisplayName','Mean Moving Median')
plot(dt_utc(1:length(S_bestguess)),S_bestguess,'.k','DisplayName','"Best guess"')
plot(discreteS.datetime_utc(ind_site),discreteS.Lab_Salinity(ind_site),'o','Color',rgb('gold'),'MarkerSize',6,'LineWidth',2,'DisplayName','Lab S')
xline(dt_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Salinity (psu)')
title([site,' - "Best Guess" Salinity'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)
legend('show')

% % Not using moving median
% S_all = [S_bc,S_erdc];
% S_mean = mean(S_all,2,'omitmissing');
% 
% fig1 = figure(1);clf
% fig1.WindowState = 'maximized';
% plot(dt_utc,S_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC');
% hold on
% plot(dt_utc,S_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
% plot(dt_utc,S_mean,'-','Color',rgb('magenta'),'DisplayName','Mean Value')
% plot(discreteS.datetime_utc(ind_site),discreteS.Lab_Salinity(ind_site),'o','Color',rgb('gold'),'MarkerSize',6,'LineWidth',2,'DisplayName','Lab S')
% xline(dt_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('Salinity (psu)')
% legend('show','location','best')
% title([site,' - After Initial QC and Moving Median Test'])
% xlim([dt1 dt2])                 % Use same x limits for comparing sites
% set(gca,'FontSize',fontsize)
% 
% switch site
%     case 'Gull'
%         % Funky deployments
%         % Dep 6 -- define points where sensors do not track together
%         kk = find(ind_diverge > 29223 & ind_diverge < dep.ind(5));
%         ind_d6 = ind_diverge(kk);
% 
%         % Choose "better" deployment
%         d1 = S_erdc(dep.ind(1):dep.ind(2));        % Dep 1
%         d2 = S_mean(dep.ind(2)+1:dep.ind(3));      % Dep 2
%         d5 = S_mean(dep.ind(3)+1:dep.ind(4));      % Dep 5
%         d6 = [S_mean(dep.ind(4)+1:ind_d6(1)); S_erdc(ind_d6(1)+1:dep.ind(5))];  % Dep 6
%         d7 = S_mean(dep.ind(5)+1:dep.ind(6));      % Dep 7
%         d8 = [S_mean(dep.ind(6)+1:41331); S_erdc(41332:dep.ind(7))];        % Dep 8
%         d9 = S_erdc(dep.ind(7)+1:dep.ind(8));      % Dep 9
%         d10 = S_erdc(dep.ind(8)+1:dep.ind(9));       % Dep 10
%         d11 = S_mean(dep.ind(9)+1:dep.ind(10));    % Dep 11
%         d12 = S_mean(dep.ind(10)+1:dep.ind(11));   % Dep 12
%         d13 = [S_erdc(dep.ind(11)+1:81445); S_erdc(81446:dep.ind(12))];   % Dep 13
%         d14 = S_erdc(dep.ind(12)+1:dep.ind(13));     % Dep 14
%         d15 = S_bc(dep.ind(13)+1:dep.ind(14));     % Dep 15
%         d16 = [S_erdc(dep.ind(14)+1:111691); S_bc(111692:dep.ind(15))];     % Dep 16
%         d17 = S_mean(dep.ind(15)+1:dep.ind(16));   % Dep 17
%         d18 = S_bc(dep.ind(16)+1:end);           % Dep 18
% 
%         S_bestguess = [d1;d2;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];
% 
%         % Calculate mean difference for Deployments 11 & 12
%         diff_S = diff(S_all(dep.ind(9):dep.ind(11),:),1,2);
%         two_sigma.S = mean(abs(diff_S),'omitmissing');
% 
%     case 'North'
%         % Funky deployments
%         % Dep 7 -- offset adjustment
%         delta_d7 = S_bc(20745) - S_erdc(20786);
% 
%         % Dep 15 -- offset adjustment
%         delta_d15 = S_erdc(76160) - S_erdc(76159);
% 
%         % Choose "better" deployment
%         d2 = [S_erdc(dep.ind(1):dep.ind(2))];          % Dep 2
%         d6 = S_bc(dep.ind(2)+1:dep.ind(3));          % Dep 6
%         d7 = S_bc(dep.ind(3)+1:dep.ind(4)) - delta_d7; % Dep 7
%         d8 = S_erdc(dep.ind(4)+1:dep.ind(5));        % Dep 8
%         d9 = S_erdc(dep.ind(5)+1:dep.ind(6));        % Dep 9
%         d10 = S_erdc(dep.ind(6)+1:dep.ind(7));       % Dep 10
%         d11 = S_mean(dep.ind(7)+1:dep.ind(8));       % Dep 11
%         d12 = S_erdc(dep.ind(8)+1:dep.ind(9));       % Dep 12
%         d13 = S_erdc(dep.ind(9)+1:dep.ind(10));      % Dep 13
%         d14 = S_erdc(dep.ind(10)+1:dep.ind(11));     % Dep 14
%         d15 = S_erdc(dep.ind(11)+1:dep.ind(12)) - delta_d15; % Dep 15
%         d16 = S_mean(dep.ind(12)+1:dep.ind(13));     % Dep 16
%         d17 = S_mean(dep.ind(13)+1:dep.ind(14));     % Dep 17
%         d18 = S_erdc(dep.ind(14)+1:end);             % Dep 18
% 
%         S_bestguess = [d2;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];
% 
%         % Remove "obvious" outliers
%         S_bestguess(4351:dep.ind(2)) = NaN;
%         S_bestguess(11090:15492) = NaN;
%         S_bestguess(20785) = NaN;
%         S_bestguess(80182:dep.ind(12)) = NaN;
%         S_bestguess(84381:84385) = NaN;
% 
%         % Calculate mean difference for Deployments 16 & 17
%         diff_S = diff(S_all(dep.ind(12):dep.ind(14),:),1,2);
%         two_sigma.S = mean(abs(diff_S),'omitmissing');
% 
%     case 'South'
%         % Funky deployments
%         % Dep 6 -- offset adjustment
%         ind = dep.ind(5);
%         delta_d6 = S_bc(ind+1) - S_bc(ind-1);
% 
%         % Dep 7 -- offset adjustment
%         ind = dep.ind(7);
%         delta_d7 = S_erdc(ind-11) - S_erdc(ind);
% 
%         % Dep 15 -- offset adjustment
%         delta_d15 = S_erdc(98253) - S_erdc(98251);
% 
%         % Choose "better" deployment
%         d1 = S_erdc(dep.ind(1):dep.ind(2));        % Dep 1
%         d2 = S_mean(dep.ind(2)+1:dep.ind(3));      % Dep 2
%         d4 = S_bc(dep.ind(3)+1:dep.ind(4));        % Dep 4
%         d5 = S_bc(dep.ind(4)+1:dep.ind(5));        % Dep 5
%         d6 = S_bc(dep.ind(5)+1:dep.ind(6)) - delta_d6;        % Dep 6
%         d7 = S_erdc(dep.ind(6)+1:dep.ind(7)) - delta_d7;      % Dep 7
%         d8 = S_erdc(dep.ind(7)+1:dep.ind(8));      % Dep 8
%         d9 = S_mean(dep.ind(8)+1:dep.ind(9));      % Dep 9
%         d10 = S_bc(dep.ind(9)+1:dep.ind(10));      % Dep 10
%         d11 = S_erdc(dep.ind(10)+1:dep.ind(11));   % Dep 11
%         d12 = S_mean(dep.ind(11)+1:dep.ind(12));   % Dep 12
%         d13 = S_erdc(dep.ind(12)+1:dep.ind(13));   % Dep 13
%         d14 = S_mean(dep.ind(13)+1:dep.ind(14));   % Dep 14
%         d15 = S_erdc(dep.ind(14)+1:dep.ind(15)) + delta_d15;   % Dep 15
%         d16 = S_mean(dep.ind(15)+1:dep.ind(16));   % Dep 16
%         d17 = S_erdc(dep.ind(16)+1:dep.ind(17));   % Dep 17
%         d18 = S_mean(dep.ind(17)+1:end);           % Dep 18
% 
%         S_bestguess = [d1;d2;d4;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17;d18];
% 
%         % Remove "obvious" outliers
%         S_bestguess(5896:5898) = NaN;
%         S_bestguess(11043:dep.ind(3)) = NaN;
%         S_bestguess(11928) = NaN;
%         S_bestguess(29482) = NaN;
%         S_bestguess(33012) = NaN;
%         S_bestguess(38194:38204) = NaN;
%         S_bestguess(49118:49149) = NaN;
%         S_bestguess(77744:77797) = NaN;
%         S_bestguess(80846:80852) = NaN;
% 
%         % Calculate mean absolute difference for Deployments 9 & 12
%         diff_S = diff(S_all([dep.ind(8):dep.ind(9),dep.ind(11):dep.ind(12)],:),1,2);
%         two_sigma.S = mean(abs(diff_S),'omitmissing');
% end
% 
% clearvars d1 d2 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 d17 d18
% 
% fig2 = figure(2);clf
% fig2.WindowState = 'maximized';
% plot(dt_utc,S_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC');
% hold on
% plot(dt_utc,S_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
% plot(dt_utc,S_mean,'.','Color',rgb('magenta'),'DisplayName','Mean Value')
% plot(dt_utc(1:length(S_bestguess)),S_bestguess,'.k','DisplayName','"Best guess"')
% plot(discreteS.datetime_utc(ind_site),discreteS.Lab_Salinity(ind_site),'o','Color',rgb('gold'),'MarkerSize',6,'LineWidth',2,'DisplayName','Lab S')
% xline(dt_utc(dep.ind),'--',label,'HandleVisibility','off')
% ylabel('Salinity (psu)')
% title([site,' - "Best Guess" Salinity'])
% xlim([dt1 dt2])                 % Use same x limits for comparing sites
% set(gca,'FontSize',fontsize)
% legend('show')

disp('Press enter to continue on to Depth')
pause

%====DEPTH=================================================================
cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\depth'])
% Check agreement between BC and ERDC sondes
depth_bc = dat_syn.depth_dat1;
depth_erdc = dat_syn.depth_dat2;
threshold = 0.25; 
ind_diverge = find(abs(depth_bc - depth_erdc) > threshold);

pc_exceed = height(ind_diverge)/height(dat_syn)*100; % Percentage of synchronized data that exceed threshold

txt = ['Differs by >',num2str(threshold),' m (',num2str(pc_exceed,'%.1f'),'%)'];

window = 144;

% Find mean depth
depth_all = [depth_bc,depth_erdc];
depth_mean = mean(depth_all,2,'omitmissing');
depth_bestguess = depth_mean;

% Find mean offset b/w BC & ERDC sondes
bc_mean = mean(depth_bc,'omitmissing');     % [m]
erdc_mean = mean(depth_erdc,'omitmissing'); % [m]
delta = abs(bc_mean - erdc_mean);           % [m]

fig3 = figure(3);clf
fig3.WindowState = 'maximized';
plot(dt_utc,depth_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC');
hold on
plot(dt_utc,depth_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
plot(dt_utc(ind_diverge),depth_bc(ind_diverge),'o','MarkerSize',circlesize,'Color',rgb('chartreuse'),'DisplayName',txt)
xline(dt_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Depth (m)')
legend('show','location','best')
title([site,' - After Initial QC and Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits

switch site
    case 'Gull'
        % Funky deployments
        % Dep 7 -- replace the mean depth at certain points
        kk = find(ind_diverge > 35479 & ind_diverge < dep.ind(6));
        ind_d7 = ind_diverge(kk);     
        depth_bestguess(ind_d7,1) = depth_bc(ind_d7) + delta/2;

    case 'North'
        % Funky deployments
        % Dep 13 -- replace the mean depth for entire deployment
        ind_d13 = (dep.ind(9):dep.ind(10))';
        depth_bestguess(ind_d13) = depth_bc(ind_d13) + delta/2;
        
        % Remove "obvious" outliers
        depth_bestguess(30420:30429) = NaN;

    case 'South'
        % Funky deployments       
        % Replace the mean depth for these deployments
        % Dep 4
        depth_bestguess(dep.ind(3):dep.ind(4)) = depth_bc(dep.ind(3):dep.ind(4)) + delta/2;
        % Dep 6
        depth_bestguess(dep.ind(5):dep.ind(6)) = depth_bc(dep.ind(5):dep.ind(6)) + delta/2;
        % Dep 10
        depth_bestguess(54311:dep.ind(10)) = depth_bc(54311:dep.ind(10)) + delta/2;
        % Dep 11
        depth_bestguess(dep.ind(10):dep.ind(11)) = depth_erdc(dep.ind(10):dep.ind(11)) - delta/2;
        % Dep 14
        depth_bestguess(dep.ind(13):dep.ind(14)) = depth_erdc(dep.ind(13):dep.ind(14)) - delta/2;
        % Dep 15
        depth_bestguess(dep.ind(14):dep.ind(15)) = depth_erdc(dep.ind(14):dep.ind(15)) - delta/2;
        % Dep 16
        depth_bestguess(110977:dep.ind(16)) = depth_bc(110977:dep.ind(16)) + delta/2;

        % Remove "obvious" outliers
        depth_bestguess(47969:47971) = NaN;

end

fig4 = figure(4);clf
fig4.WindowState = 'maximized';
plot(dt_utc,depth_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC');
hold on
plot(dt_utc,depth_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
plot(dt_utc,depth_mean,'-','Color',rgb('magenta'),'DisplayName','Mean Value')
plot(dt_utc(1:length(depth_bestguess)),depth_bestguess,'-k','DisplayName','"Best guess"')
xline(dt_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Depth (m)')
title([site,' - "Best Guess" Depth'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)
legend('show')
%%
disp('Press enter to continue on to Temperature')
pause

%====TEMPERATURE===========================================================
cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\temperature'])
% Check agreement between BC and ERDC sondes
T_bc = dat_syn.temperature_dat1;
T_erdc = dat_syn.temperature_dat2;
threshold = 0.5; 
ind_diverge = find(abs(T_bc - T_erdc) > threshold);

pc_exceed = height(ind_diverge)/height(dat_syn)*100; % Percentage of synchronized data that exceed threshold

txt = ['Differs by >',num2str(threshold),'^oC (',num2str(pc_exceed,'%.1f'),'%)'];

T_all = [T_bc,T_erdc];
T_mean = mean(T_all,2,'omitmissing');
T_bestguess = T_mean;

switch site
    case 'Gull'
        % Calculate mean difference for Deployments 11 & 12
        diff_T = diff(T_all(dep.ind(9):dep.ind(11),:),1,2);
        two_sigma.T = mean(abs(diff_T),'omitmissing');

    case 'North'
        T_bestguess(30420:30430) = NaN;

        % Calculate mean difference for Deployments 11 & 12
        diff_T = diff(T_all(dep.ind(7):dep.ind(9),:),1,2);
        two_sigma.T = mean(abs(diff_T),'omitmissing');

    case 'South'
        % Remove "obvious" outliers
        T_bestguess(47968:47970) = NaN;

        % Calculate mean difference for Deployments 11 & 12
        diff_T = diff(T_all(dep.ind(10):dep.ind(11),:),1,2);
        two_sigma.T = mean(abs(diff_T),'omitmissing');
end

fig5 = figure(5);clf
fig5.WindowState = 'maximized';
plot(dt_utc,T_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC');
hold on
plot(dt_utc,T_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
xline(dt_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Temperature (^oC)')
legend('show','location','best')
title([site,' - After Initial QC and Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

fig6 = figure(6);clf
fig6.WindowState = 'maximized';
plot(dt_utc,T_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC');
hold on
plot(dt_utc,T_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
plot(dt_utc,T_mean,'-','Color',rgb('magenta'),'DisplayName','Mean Value')
plot(dt_utc(1:length(T_bestguess)),T_bestguess,'-k','DisplayName','"Best guess"')
xline(dt_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Temperature (^oC)')
title([site,' - "Best Guess" Temperature'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)
legend('show')

%==========================================================================
% Save "best-guess" data into a structure
%==========================================================================
bestguess.deployment = table(dat_syn.datetime_utc,dat_syn.deployment_dat1,'VariableNames',{'datetime_utc','deployment'});
bestguess.depth = table(dat_syn.datetime_utc,depth_bestguess,'VariableNames',{'datetime_utc','depth'});
bestguess.salinity = table(dat_syn.datetime_utc,S_bestguess,'VariableNames',{'datetime_utc','salinity'});
bestguess.temperature = table(dat_syn.datetime_utc,T_bestguess,'VariableNames',{'datetime_utc','temperature'});

bestguess.deployment = table2timetable(bestguess.deployment);
bestguess.depth = table2timetable(bestguess.depth);
bestguess.salinity = table2timetable(bestguess.salinity);
bestguess.temperature = table2timetable(bestguess.temperature);

%====Save "best guess" data================================================
option = questdlg('Save best guess data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\',site,'\cleaned\dupcheck'])
        save([site,'-bestGuess.mat'],'bestguess','two_sigma')
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end

%====Save the plots=======================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\salinity'])
        saveas(fig1,'duplicateComparison.png')
        saveas(fig1,'duplicateComparison.fig')
        saveas(fig2,'bestGuessComparison.png')
        saveas(fig2,'bestGuessComparison.fig')
        cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\depth'])
        saveas(fig3,'duplicateComparison.png')
        saveas(fig3,'duplicateComparison.fig')
        saveas(fig4,'bestGuessComparison.png')
        saveas(fig4,'bestGuessComparison.fig')
        cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\temperature'])
        saveas(fig5,'duplicateComparison.png')
        saveas(fig5,'duplicateComparison.fig')
        saveas(fig6,'bestGuessComparison.png')
        saveas(fig6,'bestGuessComparison.fig')
        disp('Plots saved!')
    
    case 'No'
        disp('Plots not saved.')
end
