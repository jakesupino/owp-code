clear all;close all;clc

site = 'gull';

cd(['H:\My Drive\Postdoc\SMIIL\raw-data\open-water-platform-data\',site])

ds = fileDatastore(['H:\My Drive\Postdoc\SMIIL\raw-data\open-water-platform-data\',site],"ReadFcn",@load,"FileExtensions",".mat");

dat = readall(ds);

sonde1_all = table();
sonde2_all = table();
for i = 1:length(dat)
    sonde1_all = [sonde1_all;dat{i}.sonde1];
    sonde2_all = [sonde2_all;dat{i}.sonde2];
end

%%
cd(['H:\My Drive\Postdoc\SMIIL\figures\open-water-platform-figures\',site])

% Global plotting settings
FontSize = 12;
NumTicks = 13;
LineWidth = 1;
XTick = linspace(sonde1_all.datetime_utc(1),sonde1_all.datetime_utc(end),NumTicks);
XTickFormat = "M/yy";
Legend = {'BC','ERDC'};
XLabel = 'Month/Year';
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde

figure,clf
plot(sonde1_all.datetime_utc,sonde1_all.depth,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.depth,'Color',blue);
hold off
ax1 = gcf().CurrentAxes;
ylabel('Depth (m)')
xlim('tight');ylim('tight')
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure,clf
plot(sonde1_all.datetime_utc,sonde1_all.salinity,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.salinity,'Color',blue)
hold off
ax2 = gcf().CurrentAxes;
ylabel('Salinity (PSU)')
xlim('tight');ylim('tight')
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure,clf
plot(sonde1_all.datetime_utc,sonde1_all.temperature,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.temperature,'Color',blue)
hold off
ax3 = gcf().CurrentAxes;
ylabel('Temperature (^oC)')
xlim('tight');ylim('tight')
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure,clf
plot(sonde2_all.datetime_utc,sonde2_all.turbidity,'Color',blue)
ax4 = gcf().CurrentAxes;
ylabel('Turbidity (NTU)')
xlim('tight');ylim('tight')
xtickformat(XTickFormat)
legend(Legend{2})
xlabel(XLabel)

figure,clf
plot(sonde1_all.datetime_utc,sonde1_all.pH,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.pH,'Color',blue)
hold off
ax5 = gcf().CurrentAxes;
ylabel('pH')
xlim('tight');ylim('tight')
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure,clf
plot(sonde1_all.datetime_utc,sonde1_all.DO_conc,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.DO_conc,'Color',blue)
hold off
ax6 = gcf().CurrentAxes;
ylabel('DO concentration (\mumol/L)')
xlim('tight');ylim('tight')
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure,clf
plot(sonde1_all.datetime_utc,sonde1_all.ORP,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.ORP,'Color',blue)
hold off
ax7 = gcf().CurrentAxes;
ylabel('ORP (mV)')
xlim('tight');ylim('tight')
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure,clf
plot(sonde1_all.datetime_utc,sonde1_all.chla,'Color',red)
ax8 = gcf().CurrentAxes;
ylabel('Chl a (RFU)')
xlim('tight');ylim('tight')
xtickformat(XTickFormat)
legend(Legend{1})
xlabel(XLabel)

% Set these properties for all figures
set([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'FontSize',FontSize,'LineWidth',LineWidth,'XTick',XTick)