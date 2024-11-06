function plotFCTD_profiles_AAA(FCTD,figpath)
% This function plots pressure, conductivity, and temperature vs time for
% the full fastCTD dataset as an initial check to make sure the right data
% is loaded.
%
% Alex Andriatis
% 2021-08-19

fpath = fullfile(figpath,'profiles');
mkdir(fpath);


% The whole timeseries proved to be really big so lets break it up into
% 100-profile segments

seglength=100;
tmp=FCTD.up;
nprofs = length(tmp.time);
nfigs = ceil(nprofs/seglength);

sz = 25;
figure
set(gcf,'Position',[382,143,1489,1446]);
for n=1:nfigs
    clf;
    profs=[(n-1)*seglength+1:1:min(n*seglength,nprofs)];
    
    fname = ['FCTDprofiles_' num2str(profs(1))];
    
    ax(1)=subplot(2,2,1);
        hold on;
    tmp = FCTD.up;
        for t=profs(1):profs(end)
            y=tmp.pressure{t};
            % x=mean(tmp.time{t})*ones(size(y)); % a time x-axis 
            x=t*ones(size(y)); % x-axis of profile index
            z=tmp.temperature{t};
            scatter(x,y,sz,z,'filled')
        end
        c=colorbar;
        cmocean('thermal');
        set(gca,'ydir','reverse');
        xlabel('Profile number');
        ylabel('Pressure');
        title(['Upcast Temperature for profiles ' num2str(profs(1)) ' to ' num2str(profs(end))]);
        
    ax(2)=subplot(2,2,3);
        hold on;
    tmp = FCTD.down;
        for t=profs(1):profs(end)
            y=tmp.pressure{t};
            % x=mean(tmp.time{t})*ones(size(y)); % a time x-axis 
            x=t*ones(size(y)); % x-axis of profile index
            z=tmp.temperature{t};
            scatter(x,y,sz,z,'filled')
        end
        c=colorbar;
        cmocean('thermal');
        set(gca,'ydir','reverse');
        xlabel('Profile number');
        ylabel('Pressure');
        title(['Downcast Temperature for profiles ' num2str(profs(1)) ' to ' num2str(profs(end))]);
        
    ax(3)=subplot(2,2,2);
        hold on;
    tmp = FCTD.up;
        for t=profs(1):profs(end)
            y=tmp.pressure{t};
            % x=mean(tmp.time{t})*ones(size(y)); % a time x-axis 
            x=t*ones(size(y)); % x-axis of profile index
            z=tmp.conductivity{t};
            scatter(x,y,sz,z,'filled')
        end
        c=colorbar;
        cmocean('haline');
        set(gca,'ydir','reverse');
        xlabel('Profile number');
        ylabel('Pressure');
        title(['Upcast Conductivity for profiles ' num2str(profs(1)) ' to ' num2str(profs(end))]);
        
    ax(4)=subplot(2,2,4);
        hold on;
    tmp = FCTD.down;
        for t=profs(1):profs(end)
            y=tmp.pressure{t};
            % x=mean(tmp.time{t})*ones(size(y)); % a time x-axis 
            x=t*ones(size(y)); % x-axis of profile index
            z=tmp.conductivity{t};
            scatter(x,y,sz,z,'filled')
        end
        c=colorbar;
        cmocean('haline');
        set(gca,'ydir','reverse');
        xlabel('Profile number');
        ylabel('Pressure');
        title(['Downcast Conductivity for profiles ' num2str(profs(1)) ' to ' num2str(profs(end))]);
        
     linkaxes;
     
    if exist(fullfile(fpath,[fname '.png']),'file')
        warning('Figure already exists, not saved');
    else
        saveas(gcf,fullfile(fpath,fname),'png');
        saveas(gcf,fullfile(fpath,fname),'fig');
    end
     
end
     disp('Finished plotting profiles');
end
     