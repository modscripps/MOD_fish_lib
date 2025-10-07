function [] = n_savepng(fig_name)

eval(['export_fig ' fig_name ' -png -r150 -nocrop'])