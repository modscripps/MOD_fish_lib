% plot_fctd_sections

%% Decide what to plot in Fig 1 - sections
fluor_cols = get_fluorometer_colors;

axdata{1} = 'temperature';
axdata{2} = 'salinity';
axdata{3} = 'chi';
axdata{4} = 'N2';
axdata{5} = 'fluorometer';

clims.temperature = [9 18];
clims.salinity = [33.5 34.5];
clims.chi = [-9 -6];
clims.n2 = [-6 -3];
clims.fluoromter = fluor_cols.cblimits;

ylims = [0 200];

%%
fctd_mat_dir = fullfile(ec.Meta_Data.paths.data,'fctd_mat');

%% First, concatenate the individual fctd files in the deployment directory
%[FCTDall,FCTDgrid] = concatenate_and_grid_fctd(fctd_mat_dir,vars2grid_list);
if ~exist('FCTDgrid')
    disp('No FCTDgrid to plot')
else
    if ~isempty(FCTDgrid)
        %% Plot some stuff
        %clf
        %%ax = [];
        fig = figure(501);
        clf(fig)
        fig.Units = 'normalized';
        fig.Position = [0.2848    0.0370    0.4840    0.8667];
        zlim = [input_struct.depth_array(1),input_struct.depth_array(end)];
        % clim_temp = [18 27];
        % clim_sal = [34.5 35];
        % % clim_chi = [0.2 1];
        % clim_chla = [0 100];
        % clim_chi = [-10 -6];
        levels_dens = [19:1:25.5 26:0.2:27.7 27.71:0.01:27.8];

        % which data to plot? How about the most recent 1 day
        iplot=find(FCTDgrid.time>FCTDgrid.time-1);
        iplot=iplot(2:2:end); % just the up-casts till we correct the hysteresis later

        if length(iplot)>1
            % Temperature
            ax(1) = subtightplot(5,1,1);
            pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.temperature(:,iplot));
            %imagesc(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.temperature(:,iplot));
            ax(1).CLim = clims.temperature;
            cb(1) = colorbar;
            colormap(ax(1),lansey)
            cb(1).Label.String = 'Temperature';

            % Salinity
            ax(2) = subtightplot(5,1,2);
            pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,real(FCTDgrid.salinity(:,iplot)));
            ax(2).CLim = clims.salinity;
            cb(2) = colorbar;
            colormap(ax(2),cmocean('delta'));
            cb(2).Label.String = 'Salinity';
            %set(cb(2),'ydir','reverse');

            % chi
            ax(3) = subtightplot(5,1,3);
            pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi(:,iplot)));
            %imagesc(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi(:,iplot)));
            % pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi2(:,iplot)));
            ax(3).CLim = clims.chi;
            cb(3) = colorbar;
            colormap(ax(3),parula)
            cb(3).Label.String = '\chi';

            switch axdata{3}
                case 'chla'
                    if isfield(FCTDgrid,'chla')
                        pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,((FCTDgrid.chla(:,iplot))/2^16-0.5)*500.0);
                        %pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.chla(:,iplot));
                        ax(3).CLim = [0 3];
                        cb(3) = colorbar;
                        cb(3).Label.String = 'Chla';
                    end
                case 'fluor'
                    if isfield(FCTDgrid,'fluor')
                        pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.fluor(:,iplot));
                        ax(3).CLim = [3.25e4 3.33e4];
                        cb(3) = colorbar;
                        cb(3).Label.String = 'fluor';
                    end
            end

            ax(4) = subtightplot(5,1,4);
            switch axdata{4}
                case 'chi'
                    pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi(:,iplot)));
                    colormap(ax(4),cmocean('matter'))
                    cb(4) = colorbar;
                    ax(4).CLim = clims.chi;
                    cb(4).Label.String = 'chi';
                case 'chi2'
                    pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi2(:,iplot)));
                    colormap(ax(4),cmocean('matter'))
                    cb(4) = colorbar;
                    ax(4).CLim = clims.chi;
                    cb(4).Label.String = 'chi2';
                case 'N2'
                    % N^2
                    %pcolorjw(FCTDgrid.time,FCTDgrid.depth,((FCTDgrid.chla)/2^16-0.5)*500.0);
                    %ALB adiabtic sorting to plot N2 adia
                    sort_T    = FCTDgrid.temperature.*nan;
                    sort_S    = FCTDgrid.salinity.*nan;
                    OT        = FCTDgrid.salinity.*nan;
                    for p=1:length(FCTDgrid.time)
                        local_dens=FCTDgrid.density(:,p);
                        sort_dens=FCTDgrid.density(:,p).*nan;
                        dens_Inan=find(~isnan(local_dens));
                        [~,IA]=sort(local_dens(dens_Inan),'ascend');
                        sort_dens(dens_Inan)=local_dens(dens_Inan(IA));
                        delta_dens=local_dens-sort_dens(:);
                        delta_dens(delta_dens==0)=nan;

                        sort_T(dens_Inan,p)=FCTDgrid.temperature(dens_Inan(IA),p);
                        sort_S(dens_Inan,p)=FCTDgrid.salinity(dens_Inan(IA),p);
                        OT(:,p)=abs(delta_dens)>.002;
                    end
                    [bfrq,vort,p_ave] = sw_bfrq(FCTDgrid.salinity,FCTDgrid.temperature,FCTDgrid.pressure,mean(FCTDgrid.latitude,'omitmissing'));
                    Ndepth=FCTDgrid.depth(1:end-1)+diff(FCTDgrid.depth);
                    % [bfrq,vort,p_ave] = sw_bfrq(sort_S,sort_T,FCTDgrid.pressure,mean(FCTDgrid.latitude,'omitmissing'));
                    pcolorjw(FCTDgrid.time(iplot),Ndepth,real(log10(bfrq(:,iplot))));
                    colormap(ax(4),cmocean('speed'))
                    ax(4).CLim = (clims.n2);
                    cb(4) = colorbar;
                    cb(4).Label.String = 'N^2';


            end % End ax4

            % Axes 5 - fluorometer
            ax(5) = subtightplot(5,1,5);

            % We're going to plot fluorometer data by "painting by numbers". Get
            % the color the corresponds to each fluorometer bin.
            fluorometer = FCTDall.fluorometer(:,1);
            fluorBins = nan(size(fluorometer));
            for iRange=1:length(fluor_cols.levels)-1
                fluorBins(fluorometer > fluor_cols.levels_actual(iRange) &...
                    fluorometer <= fluor_cols.levels_actual(iRange+1)) = iRange;
            end
            scatter(FCTDall.time,FCTDall.depth,20,fluorometer,'filled');
            

            % Apply colormap
            % cb(5) = colorbar;
            % colormap(ax(5),fluor_cols.cmap);
            % set(gca,'clim',[fluor_cols.levels(1),fluor_cols.levels(end)]);
            % cb(5).Ticks = fluor_cols.ticks;
            % cb(5).TickLabels = fluor_cols.ticklabels;
            % cb(5).Label.String = fluor_cols.units;
            cb(5) = colorbar;
            set(gca,'clim',clims.fluorometer)
            colormap(gca,lansey)

            % Add density contours and datetick
            for iAx=1:5
                axes(ax(iAx))
                hold(ax(iAx),'on')
                [c,ch] = contour(FCTDgrid.time,FCTDgrid.depth,real(FCTDgrid.density-1000),['k'],'levellist',levels_dens);
                %contour(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature,'m','levellist',13);
                %  contour(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature,'c','levellist',15);
                clabel(c,ch);

                datetick(ax(iAx),'x','HH:MM','keeplimits')
            end
            [ax(1:4).XTickLabel] = deal('');
            [ax(:).XGrid] = deal('on');

            % Depth axes
            [ax(1:5).YLim] = deal(ylims);
            [ax(1:5).YDir] = deal('reverse');

            % Link axes
            lp = linkprop([ax(:)],{'xlim','ylim'});
            drawnow
            pos1 = get(ax(1),'position');
            ax(5).Position(3)= pos1(3);

            %MHA hack
            %[ax(1:2).YLim] = deal([zlim(1) 500]);
        end %end if iplot>1

        %% Plot TS of last few profiles
        openFigs = findobj('type','figure');
        iTS = find(contains({openFigs(:).Tag},'ts_plot'));
        if ~isempty(iTS)
            close(openFigs(iTS))
        end

        firstProf = max([length(FCTDgrid.time)-4,1]);
        profList = firstProf:length(FCTDgrid.time);
        cols = cmocean('thermal',length(profList));
        figTS = figure(2);
        figTS.Tag = 'ts_plot';
        for pp=1:length(profList)
            p = profList(pp);
            if mod(p,2)==0
                plot(FCTDgrid.salinity(:,p),FCTDgrid.temperature(:,p),'.','MarkerSize',12,'Color',cols(pp,:),'displayname',sprintf('Down %i',p/2));
            else
                plot(FCTDgrid.salinity(:,p),FCTDgrid.temperature(:,p),'.','MarkerSize',12,'Color',cols(pp,:),'displayname',sprintf('Up %i',floor(p/2)));
            end
            hold on
        end
        grid on
        xlim([33.2 35]);
        ylim([8 18])
        xlabel('S [psu]','FontName','times new roman','FontSize',13);
        ylabel('T [˚C]','FontName','times new roman','FontSize',13);
        legend

    else
        disp('No FCTDgrid to plot')
    end %end of ~isempty(FCTDgrid)
end %end if exist(FCTDgrid)


function [fluor_cols] = get_fluorometer_colors

blue = [80 117 186]./255;
%clevels = unique([0,0.2:0.02:0.4, 0.4:0.05:1.2, 1:1:10, 30:30:150]);
clevels = unique([0,linspace(0,0.2,50),linspace(0.2,1,20),linspace(1,2.5,15)]);
matter = cmocean('matter',length(clevels)+4);

clevel_ticks = [0.2,1,2.5];
for t=1:length(clevel_ticks)
    cb_ticks(t) = find(clevels==clevel_ticks(t));
end

fluor_cols.cmap = [blue; matter(1:end-5,:)];
fluor_cols.levels_actual = clevels;
fluor_cols.levels = 1:length(clevels);
fluor_cols.ticks = cb_ticks;
fluor_cols.ticklabels = fluor_cols.levels_actual(cb_ticks);
fluor_cols.cblimits = [0,length(clevels)];
fluor_cols.cblimits_actual = [0,2.5];
fluor_cols.units = 'fluor (V)';

end %end get_fluorometer_colors