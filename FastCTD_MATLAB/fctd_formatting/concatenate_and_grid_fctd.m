function [FCTDall,FCTDgrid] = concatenate_and_grid_fctd(fctd_mat_dir)
% [FCTDall,FCTDgrid] = concatenate_and_grid_fctd(fctd_mat_dir)
%
% Concatenates all .modraw FCTD files in the directory 'fctd_mat_dir' into FCTDall.
% Grids FCTDall into FCTDgrid.

% To do: First check that any of these exist. Don't remake them all the
% time!
FCTDall = make_FCTDall_L0(fctd_mat_dir);
FCTDall = make_FCTDall_L1(FCTDall,fctd_mat_dir);
FCTDgrid = make_FCTDgrid(FCTDall,fctd_mat_dir);

end %end main function




% -------------------------------------------------------------------------
