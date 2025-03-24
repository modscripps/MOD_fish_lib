function [FCTDall,FCTDgrid] = concatenate_and_grid_fctd(fctd_mat_dir,vars2grid_list)
% [FCTDall,FCTDgrid] = concatenate_and_grid_fctd(fctd_mat_dir,vars2grid_list)
%
% Concatenates all .modraw FCTD files in the directory 'fctd_mat_dir' into FCTDall.
% Grids FCTDall into FCTDgrid.


FCTDall = make_FCTDall_L0(fctd_mat_dir);
FCTDall = make_FCTDall_L1(FCTDall,fctd_mat_dir);
FCTDgrid = grid_FCTDall(FCTDall);

end %end main function




% -------------------------------------------------------------------------
