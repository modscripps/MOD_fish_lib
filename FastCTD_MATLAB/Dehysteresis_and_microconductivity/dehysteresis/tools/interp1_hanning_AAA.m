function data_grid = interp1_hanning_AAA(x,y,xg)
% This function is designed to grid unevenly-spaced data onto a
% uniformly-spaced grid. It addresses the problem of improper weighting when 
% binning data by applying a running hanning window mean to the data before
% interpolating onto a uniform grid.
%
% Inputs:
%   x: The data position vector
%   y: The data at the positions x
%   xg: A uniformly-spaced target position vector
% 
% Outputs:
%   data_grid: Gridded data at the positions xg
%
% Alex Andriatis
% 2021-03-11

% The width of the gridded bins
dx = abs(mean(diff(xg)));

% Make sure data has no NaNs
I = ~isnan(y);
x = x(I);
y = y(I);

B = NaN(1,length(x));
% Loop through every datapoint
for i=1:length(x)
    dist = x-x(i); % At each datapoint calculate the distance
    A = cos(pi*dist/dx).^2; % Make a hanning window
    A(abs(dist)>(dx/2))=0; % Window is limited to the size of a grid cell
    A = A./sum(A); % Scale the window so that the sum of the weights is 1
    B(i) = sum(A.*y); % Hanning-average the data around each point
end

% Make sure the data points are unique before interpolation. This is done
% after averaging so that all the data are included.
[x,I]=unique(x);
B = B(I);

% Grid the data onto the uniform grid
data_grid = interp1(x,B,xg);
end