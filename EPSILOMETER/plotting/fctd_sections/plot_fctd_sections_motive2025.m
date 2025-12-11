% plot_fctd_sections

%% Decide what to plot in ax(4)
%ax4data = 'chla';
% ax4data = 'chi2';
ax4data = 'N2';
ax3data = 'chi';
%ax3data = 'chla';

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
        fig = figure(1);
        clf
        % fig.Units = 'normalized';
        % fig.Position = [0.1    0.0370    0.3    0.8667];
        
        % set depth to deepest observationtemp = FCTDgrid.temperature

        ttemp = FCTDgrid.temperature(~all(isnan(FCTDgrid.temperature), 2), :);
        [n, ~] = size(ttemp);
        deep_lim = 2500;%FCTDgrid.depth(n+10);


        % zlim = [input_struct.depth_array(1),input_struct.depth_array(end)];
        zlim = [0 deep_lim];
        % clim_temp = [18 27];
        clims.temperature = [22.5 25.2];
        % clim_sal = [34.5 35];
        clims.salinity = [34.7 35.4];
        % % clim_chi = [0.2 1];
        clim_chla = [0.2e-5 7e-5];
        % clim_chi = [-10 -6];
        clims.chi = [-10 -4];
        levels_dens = [19:0.25:25.5 26:0.2:27.7 27.71:0.01:27.8];
        clims.dens = [19:0.25:25.5 26:0.2:27.7 27.71:0.01:27.8];

        clims.n2 = [-6 -2.5];

        % which data to plot? How about the most recent 1 day
        iplot=find(FCTDgrid.time>FCTDgrid.time-1);
        iplot=iplot(2:2:end); % just the up-casts till we correct the hysteresis later

        if length(iplot)>1
        % Temperature
        ax(1) = subtightplot(4,2,1);
        % pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.temperature(:,iplot));
        pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.temperature(:,iplot));
        ax(1).CLim = clims.temperature;
        colormap(ax(1),lansey)

        ax(2) = subtightplot(4,2,2);
        % pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.temperature(:,iplot));
        pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.temperature(:,iplot));
        ax(2).CLim = clims.temperature;
        cb(2) = colorbar;
        colormap(ax(2),lansey)
        cb(2).Label.String = 'Temperature';

        % Salinity
        ax(3) = subtightplot(4,2,3);
        pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,real(FCTDgrid.salinity(:,iplot)));
        ax(3).CLim = clims.salinity;
        colormap(ax(3),cmocean('delta'));
        %set(cb(2),'ydir','reverse');

        ax(4) = subtightplot(4,2,4);
        pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,real(FCTDgrid.salinity(:,iplot)));
        ax(4).CLim = clims.salinity;
        cb(4) = colorbar;
        colormap(ax(4),cmocean('delta'));
        cb(4).Label.String = 'Salinity';
        %set(cb(2),'ydir','reverse');


        switch ax3data
            case 'chla'
                if isfield(FCTDgrid,'chla')
                    pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.chla(:,iplot))%/2^16-0.5)*500.0);
                    %pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.chla(:,iplot));
                    ax(3).CLim = clim_chla;
                    %ax(3).YLim([-200,0])
                    ylim([-200,0])
                    cb(3).Label.String = 'Chla';
                end
            case 'fluor'
                if isfield(FCTDgrid,'fluor')
                    pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.fluor(:,iplot));
                    ax(3).CLim = [3.25e4 3.33e4];
                end

            case 'chi'
                ax(5) = subtightplot(4,2,5);
                % pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi(:,iplot)));
                imagesc(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi(:,iplot)));
                % pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi2(:,iplot)));
                ax(5).CLim = clims.chi;
                colormap(ax(5),cmocean('matter'))

                ax(6) = subtightplot(4,2,6);
                % pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi(:,iplot)));
                imagesc(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi(:,iplot)));
                % pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi2(:,iplot)));
                ax(6).CLim = clims.chi;
                cb(6) = colorbar;
                colormap(ax(6),cmocean('matter'))
                cb(6).Label.String = '\chi';
        end

        switch ax4data
            case 'chi'
                ax(7) = subtightplot(4,2,7);
                pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi(:,iplot)));
                ax(8) = subtightplot(4,2,8);
                pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi(:,iplot)));
                colormap(ax(7),cmocean('matter'))
                colormap(ax(8),cmocean('matter'))
                ax(7).CLim = clims.chi;
                cb(8) = colorbar;
                ax(8).CLim = clims.chi;
                cb(8).Label.String = 'chi';
            case 'chi2'
                ax(7) = subtightplot(4,2,7);
                pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi2(:,iplot)));
                ax(8) = subtightplot(4,2,8);
                pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,log10(FCTDgrid.chi2(:,iplot)));
                colormap(ax(7),cmocean('matter'))
                ax(7).CLim = clims.chi;
                colormap(ax(8),cmocean('matter'))
                cb(8) = colorbar;
                ax(8).CLim = clims.chi;
                cb(8).Label.String = 'chi2';
            case 'chla'
                if isfield(FCTDgrid,'chla')
                    ax(7) = subtightplot(4,2,7);
                    pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.chla(:,iplot))%/2^16-0.5)*500.0);
                    ax(8) = subtightplot(4,2,8);
                    pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.chla(:,iplot))
                    %pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.chla(:,iplot));
                    ax(7).CLim = clim_chla;
                    ax(8).CLim = clim_chla;
                    cb(8) = colorbar;
                    cb(8).Label.String = 'Chla';
                end
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
                ax(7) = subtightplot(4,2,7);
                pcolorjw(FCTDgrid.time(iplot),Ndepth,real(log10(bfrq(:,iplot))));
                ax(8) = subtightplot(4,2,8);
                pcolorjw(FCTDgrid.time(iplot),Ndepth,real(log10(bfrq(:,iplot))));
                colormap(ax(7),cmocean('speed'))
                ax(7).CLim = (clims.n2);
                colormap(ax(8),cmocean('speed'))
                ax(8).CLim = (clims.n2);
                cb(8) = colorbar;
                cb(8).Label.String = 'N^2';


        end % End ax4

        % Add density contours and datetick
        for iAx=1:2:8
            axes(ax(iAx))
            ylabel("depth [m]")
            hold(ax(iAx),'on')
            [c,ch] = contour(FCTDgrid.time(iplot),FCTDgrid.depth,real(FCTDgrid.density(:,iplot)-1000),['k'],'levellist',levels_dens);
            %contour(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature,'m','levellist',13);
            %  contour(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature,'c','levellist',15);
            clabel(c,ch);

            datetick(ax(iAx),'x','HH:MM','keeplimits')
        end
        for iAx=2:2:8
            axes(ax(iAx))
            hold(ax(iAx),'on')
            [c,ch] = contour(FCTDgrid.time(iplot),FCTDgrid.depth,real(FCTDgrid.density(:,iplot)-1000),['k'],'levellist',levels_dens);
            %contour(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature,'m','levellist',13);
            %  contour(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature,'c','levellist',15);
            clabel(c,ch);

            datetick(ax(iAx),'x','HH:MM','keeplimits')
        end

        % Depth axes
        [ax(1:8).YLim] = deal([zlim(1) zlim(2)]);
        [ax(1:8).YDir] = deal('reverse');
        [ax(1).YLim] = deal([0,1300]);
        [ax(3).YLim] = deal([0,1300]);
        [ax(5).YLim] = deal([0,1300]);
        [ax(7).YLim] = deal([0,1300]);

        [ax(2).YLim] = deal([0,200]);
        [ax(4).YLim] = deal([0,200]);
        [ax(6).YLim] = deal([0,200]);
        [ax(8).YLim] = deal([0,200]);

        title(ax(1),'Full Depth');
        title(ax(2),'Upper 200m');

        xticklabels(ax(1), {});
        xticklabels(ax(2), {});
        xticklabels(ax(3), {});
        xticklabels(ax(4), {});
        xticklabels(ax(5), {});
        xticklabels(ax(6), {});

        yticklabels(ax(2), {});
        yticklabels(ax(4), {});
        yticklabels(ax(6), {});
        yticklabels(ax(8), {});

        % Link axes
        % lp = linkprop([ax(:)],{'xlim','ylim'});
        lp = linkprop([ax(:)],{'xlim'});

        % Save figure
        fig_name = fullfile(fig_path,['sections_',datestr(now,'yy_mmdd_HHMMSS')]);
        % n_savepng(fig_name);

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
        figTS.Units = 'normalized';
        figTS.Position = [0.02    0.50    0.35    0.5];
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
        xlim([34.5 35.0]);
        ylim([2 20]);
        xlabel('S [psu]','FontSize',13);
        ylabel('T [˚C]','FontSize',13);
        legend

        fig_name = fullfile(fig_path,['TS_',datestr(now,'yy_mmdd_HHMMSS')]);
        % n_savepng(fig_name);
    else
        disp('No FCTDgrid to plot')
    end %end of ~isempty(FCTDgrid)
end