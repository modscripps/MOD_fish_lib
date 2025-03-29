            
FCTDall.fluorometer_volts = FCTDall.fluorometer(:,1);

ax(5) = axes;

% We're going to plot fluorometer data by "painting by numbers". Get
% the color the corresponds to each fluorometer bin.
fluorBins = nan(size(FCTDall.fluorometer_volts));
for iRange=1:length(fluor_cols.levels)-1
    fluorBins(FCTDall.fluorometer_volts > fluor_cols.levels_actual(iRange) &...
        FCTDall.fluorometer_volts <= fluor_cols.levels_actual(iRange+1)) = iRange;
end
scatter(FCTDall.time,FCTDall.depth,30,fluorBins);

% Apply colormap
cb(5) = colorbar;
colormap(ax(5),fluor_cols.cmap);
set(gca,'clim',[fluor_cols.levels(1),fluor_cols.levels(end)]);
cb(5).Ticks = fluor_cols.ticks;
cb(5).TickLabels = fluor_cols.ticklabels;
cb(5).Label.String = fluor_cols.units;