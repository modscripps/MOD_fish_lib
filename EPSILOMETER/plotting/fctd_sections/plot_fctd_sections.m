% The section plots that you choose to plot for each cruise will vary. The
% colorbar limits will change based on where and when you are sampling and
% you might add extra plots like TS or maps. 
%
% Keep all the versions in
% /MOD_fish_lib/EPSILOMETER/plotting/plot_fctd_sections
%
% Here, just uncomment the one you want to use. (For now, that's easier
% than making a function like plot_fctd_sections('motive2024') because
% you'd need a lot of inputs.

fig_path = ec.Meta_Data.paths.figures;

%% TLC
if 0
    plot_fctd_sections_tlc2025
    save_sections_fig(fig)
    fprintf('plot_fctd_sections.m - Plotting FCTD sections for TLC 2025 \n')

    plot_fctd_map_tlc2025
    save_map_fig(fig)
    fprintf('plot_fctd_sections.m - Plotting map for TLC 2025 \n')
end

%% MOTIVE
if 1
    plot_fctd_sections_motive2024
    % fprintf('plot_fctd_sections.m - Plotting FCTD sections for MOTIVE 2024 \n')
end

%% SUBFUNCTIONS
function [] = save_sections_fig()

end %end save_sections_fig

function [] = save_map_fig()

end %end save_map_fig