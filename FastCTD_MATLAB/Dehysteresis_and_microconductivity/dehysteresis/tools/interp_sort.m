function [x_sort,y_sort]=interp_sort(x,y)
% interp_sort helps prepare an array for 1-D interpolation
%
% This function takes the x-axis and values used in the interpolation, sorts it
% in ascending order, and keeps only unique x-values to prevent interpolation
% problems. For repeated x-values the function returns the mean y-value
% 
% Alex Andriatis
% 2020-12-07

[x_sort,~,ic] = unique(x);
y_sort = accumarray(ic,y,[],@nanmean);
end

