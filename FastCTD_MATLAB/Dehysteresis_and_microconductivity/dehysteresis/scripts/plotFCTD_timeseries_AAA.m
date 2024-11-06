function plotFCTD_timeseries_AAA(FCTD,figpath)
% This function plots pressure, conductivity, and temperature vs time for
% the full fastCTD dataset as an initial check to make sure the right data
% is loaded.
%
% Alex Andriatis
% 2021-08-19

fpath = fullfile(figpath,'raw');
fname = 'FCTDfull';
mkdir(fpath);

figure;
    ax(1)=subplot(3,1,1);
        plot(FCTD.time,FCTD.pressure);
        ylabel('Pressure');
        set(ax(1),'YDir','reverse');
        
     ax(2)=subplot(3,1,2);
        plot(FCTD.time,FCTD.temperature);
        ylabel('Temperature');
        
     ax(3)=subplot(3,1,3);
        plot(FCTD.time,FCTD.conductivity);
        ylabel('Conductivity');
     
     dynamicDateTicks(ax,'linked');
     linkaxes(ax,'x');
     

saveas(gcf,fullfile(fpath,fname),'png');
saveas(gcf,fullfile(fpath,fname),'fig');
end
     