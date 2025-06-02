% Open basemap
basemap_dir = '/Volumes/Software_current_cruise/MOD_fish_lib/EPSILOMETER/plotting/basemaps/';
openfig(fullfile(basemap_dir,'map_50m_contours.fig'))

% We're going to plot fluorometer data by "painting by numbers". Get
% the color the corresponds to each fluorometer bin.
fluorBins = nan(size(FCTDgrid.fluorometer));
for iRange=1:length(fluor_cols.levels)-1
    fluorBins(FCTDgrid.fluorometer > fluor_cols.levels_actual(iRange) &...
        FCTDgrid.fluorometer <= fluor_cols.levels_actual(iRange+1)) = iRange;
end


% Add FCTDgrid fluor data - max value below 50 m
hold on
idxDeep = FCTDgrid.depth>50;
scatter(nanmean(FCTDgrid.longitude),nanmean(FCTDgrid.latitude),30,nanmean(fluorBins(idxDeep,:)),'filled');

% Apply colormap
cb = colorbar;
colormap(fluor_cols.cmap);
set(gca,'clim',[fluor_cols.levels(1),fluor_cols.levels(end)]);
cb.Ticks = fluor_cols.ticks;
cb.TickLabels = fluor_cols.ticklabels;
cb.Label.String = fluor_cols.units;
