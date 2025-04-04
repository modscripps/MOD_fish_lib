disp('Start plotting')


clims.temperature = [9 13];
clims.salinity = [33.5 34.5];
clims.chi = [-10 -6];
clims.epsilon = [-10 -6];
clims.n2 = [-6 -3];

%% Plot temperature, epsilon, and chi with density contours
data=load(fullfile(ec.Meta_Data.paths.profiles,'griddedProfiles'));

%% Plot profiles if you haven't yet
profile_plot_list = dir(fullfile(ec.Meta_Data.paths.figures,'Profile*.png'));

% Find the number of the last profile with a profile figure
profnum_cells = cellfun(@(x) x(8:11), {profile_plot_list(:).name}, 'UniformOutput',false);
if ~isempty(profnum_cells)
    profnum_list = str2num(cell2mat(profnum_cells.'));
else
    profnum_list = 1;
end

% Make a figure for all the profiles that haven't been plotted yet, and the
% most recent one, in case it wasn't complete when the figure was made.
for profID = max(profnum_list):length(data.GRID.dnum)


    curr_profile_name=sprintf('Profile%04i.mat',profID);
    lastProfile=load(fullfile(ec.Meta_Data.paths.profiles,curr_profile_name));

    data.GRID.bottom_depth=filloutliers(data.GRID.bottom_depth,'linear');
    close all

    fig1=figure('units','inches','position',[10         0   18 13]);
    QuickEpsiProfilePlotMHA_accel(lastProfile.Profile);
    fig1.PaperPosition=[0 0 18 13];
    print('-dpng2',fullfile(ec.Meta_Data.paths.figures,[curr_profile_name(1:end-4) '.png']))
    eval(['savefig ' fullfile(ec.Meta_Data.paths.figures,curr_profile_name(1:end-4))])
    
end

%% Plot TS of last few profiles

cols = cmocean('thermal',5);
profList = length(data.GRID.dnum)-4:length(data.GRID.dnum);
figure
for pp=1:length(profList)
    p = profList(pp);
    plot(data.GRID.s(:,p),data.GRID.t(:,p),'.','Color',cols(pp,:),'displayname',num2str(p));
    hold on
end

legend

%% Plot section
nProf = length(data.GRID.dnum);
if 5+nProf*0.1<18
    figWidth = 5+nProf*0.1;
else
    figWidth = 18;
end
figure('units','inches','position',[10         0   figWidth 13])

dnummask=find(~isnan(data.GRID.dnum));
[~,iun]=unique(data.GRID.dnum(dnummask)); dnummask=dnummask(iun); %size(dnummask)

%grid for plot
nr = 7; nc = 1;

%let's just plot the most recent 24 hours
dnumm=data.GRID.dnum(dnummask); 
irecent=find(dnumm>(dnumm(end)-1));
dnummask=dnummask(irecent);

% depth limits for plotting
zlim=[depth_array(1),depth_array(end)];
sgf=19.9; % front we are tracking 27 may

clear ax

ax(1) =subtightplot(nr,nc,1);
pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,data.GRID.s(:,dnummask));
hold on
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
caxis(clims.salinity) %TODO make it a param in RUN_Auto_process
ylim(ax(1),zlim)
cax1=colorbar;
grid(ax(1),'on');
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
title('Salinity')
ylabel('Depth','fontname','Times New Roman','fontsize',20)
ylabel(cax1,'psu','fontname','Times New Roman','fontsize',20)
set(ax(1),'ydir','reverse');
datetick(ax(1),'x','keeplimits')
colormap(gca,cmocean('delta'))

ax(2) =subtightplot(nr,nc,2);
pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,data.GRID.t(:,dnummask));
hold on
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
caxis(clims.temperature) 
cax2=colorbar;
grid(ax(2),'on');
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
title('Temperature')
ylabel('Depth','fontname','Times New Roman','fontsize',20)
ylabel(cax2,'\circC','fontname','Times New Roman','fontsize',20)
set(ax(2),'ydir','reverse');
datetick(ax(2),'x','keeplimits')
ylim(ax(2),zlim)
colormap(gca,lansey)

ax(3) = subtightplot(nr,nc,3);
pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,log10(data.GRID.epsilon_co1(:,dnummask)));
cax3=colorbar;
title('Epsilon 1')
ylabel('Depth','fontname','Times New Roman','fontsize',20);
hold on
caxis(clims.epsilon)
ylim(ax(3),zlim)
grid(ax(3),'on');

%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
ylabel(cax3,'W kg^{-1}','fontname','Times New Roman','fontsize',20)
set(ax(3),'ydir','reverse');
datetick(ax(3),'x','keeplimits')
colormap(gca,cmocean(['curl']))

ax(4) = subtightplot(nr,nc,4);
pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,log10(data.GRID.epsilon_co2(:,dnummask)));
cax4=colorbar;
title('Epsilon 2')
ylabel('Depth','fontname','Times New Roman','fontsize',20);
hold on
caxis(clims.epsilon)
ylim(ax(4),zlim)
grid(ax(4),'on');

%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
ylabel(cax4,'W kg^{-1}','fontname','Times New Roman','fontsize',20)
set(ax(4),'ydir','reverse');
datetick(ax(4),'x','keeplimits')
colormap(gca,cmocean(['curl']))

ax(5) = subtightplot(nr,nc,5);
pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,log10(data.GRID.chi1(:,dnummask)));
hold on
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(:,dnummask))
title('Chi 1')
cax5=colorbar;
grid(ax(5),'on');
ylim(ax(5),zlim)

caxis(clims.chi)
ylabel('Depth','fontname','Times New Roman','fontsize',20)
% xlabel(datestr(nanmin(data.GRID.dnum(:)),"dd-mm"),'fontname','Times New Roman','fontsize',20)

ylabel(cax5,'K^2 s^{-1}','fontname','Times New Roman','fontsize',20)
set(ax(5),'ydir','reverse');
% lp = linkprop([ax(:)],'xlim');
% linkaxes(ax,'xy');
datetick(ax(5),'x','keeplimits')
colormap(gca,cmocean('amp'))

ax(6) = subtightplot(nr,nc,6);
pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,log10(data.GRID.chi2(:,dnummask)));
hold on
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(:,dnummask))
title('Chi 2')
cax6=colorbar;
grid(ax(6),'on');
ylim(ax(6),zlim)
caxis(clims.chi)
ylabel('Depth','fontname','Times New Roman','fontsize',20)
ylabel(cax6,'K^2 s^{-1}','fontname','Times New Roman','fontsize',20)
set(ax(6),'ydir','reverse');
colormap(gca,cmocean('amp'))

% add n2
ax(7) = subtightplot(nr,nc,7);
[bfrq,vort,p_ave] = sw_bfrq(data.GRID.s,data.GRID.t,data.GRID.pr*ones(1,length(data.GRID.dnum)),12);
pcolorjw(data.GRID.dnum(dnummask),p_ave(:,1),real(log10(bfrq(:,dnummask))));
hold on
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(:,dnummask))
title('N^2')
cax7=colorbar;
grid(ax(7),'on');
ylim(ax(7),zlim)
caxis(clims.n2)
ylabel('Depth','fontname','Times New Roman','fontsize',20)
%xlabel(datestr(nanmin(data.GRID.dnum(:)),"dd-mmm"),'fontname','Times New Roman','fontsize',20)
[ax(:).YDir] = deal('reverse');
ylabel(cax7,'N^2 s^{-2}','fontname','Times New Roman','fontsize',20)
set(ax(7),'ydir','reverse');
lp = linkprop([ax(:)],'xlim');
linkaxes(ax,'xy');
colormap(gca,cmocean('speed'))

% Add density contours
for iAx=1:length(ax)
    axes(ax(iAx))
    [c,ch]=contour(data.GRID.dnum(dnummask),data.GRID.z,real(data.GRID.sgth(:,dnummask))-1e3,[19:.1:30],'color',.5*[1 1 1]);
end

for a=1:length(ax)
    ax(a).XTick=data.GRID.dnum(1):2/24:data.GRID.dnum(end);
    ax(a).XTickLabel='';
end
try
    ax(end).XTickLabel=datestr(data.GRID.dnum(1):2/24:data.GRID.dnum(end),'DD - HH:MM');
    ax(end).XTickLabelRotation=45;
catch
    datetick(ax(end),'x','keeplimits')
end

datetick(ax(end),'x','keeplimits')
[ax(:).YLim] = deal([0 160]);

% Save section plot
save_name = fullfile(ec.Meta_Data.paths.figures,'deployment_sections');
% eval(['export_fig ' save_name ' -png -r150 -nocrop']);
print('-dpng',save_name)
eval(['savefig ' save_name])

