function [vars2grid_list] = get_FCTD_fields(FCTDall)
% Get the field names that match the number of elements in FCTDall.time. This
% could include uConductivity which is a Nx20 array with N =
% size(FCTDall.time,1)

size_time = size(FCTDall.time,1);
field_list = fields(FCTDall);
is_same_size = structfun(@(x) isequal(size(x,1), size_time), FCTDall);
vars2grid_list = field_list(is_same_size);
% Remove time from list
vars2grid_list = setdiff(vars2grid_list, 'time');