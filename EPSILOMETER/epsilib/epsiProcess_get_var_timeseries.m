function [var_timeseries] = epsiProcess_get_var_timeseries(obj,var_name)

% [var_timeseries] = epsiProcess_get_var_timeseries(obj,var_name)
% Get the full timeseries of a single variable from an epsi_class object
%
% INPUTS
%   obj - epsi_class object
%   var_name - name of variable (e.g. t1_volt,s2_volt,T,dPdt)
%
% OUTPUTS
%   var_timeseries.dnum - datenum array
%   var_timeseries.data - data array
%
% APPROACH
%   First, load just the first file and look for the variable name.
%   The input variable name is two strucutre levels deep. Loop through the
%   first level to find it. Once you do, grab the data for the whole
%   deployment and grab the dnum from the same structure.

% Load TimeIndex
load(fullfile(obj.Meta_Data.paths.mat_data,'TimeIndex'),'TimeIndex');

% NC 11/19/24 -  Epsi Minnow puts out files called modsom_1, modsom_2,
% modsom_10, etc so they look out of order when reading the list of files
% in a directory alphabetically (1, 10, 100, 101, 102, 103). TimeIndex gets out of order. Before
% cropping the time series, sort TimeIndex by dnumStart.
[~,iS] = sort(TimeIndex.dnumStart);
field_list = fields(TimeIndex);
for iF=1:length(field_list)
    TimeIndex.(field_list{iF}) = TimeIndex.(field_list{iF})(iS); 
end

% Load the first file
file_data = load(fullfile(obj.Meta_Data.paths.mat_data,TimeIndex.filenames{1}));

% Find the requested variable
field_list = fields(file_data);

for iF=1:length(field_list)
if ~isempty(file_data.(field_list{iF}))
    sub_list = fields(file_data.(field_list{iF}));
    iS = find(strcmp(sub_list,var_name));
    if ~isempty(iS)
        data = file_data.(field_list{iF}).(sub_list{iS});
        dnum = file_data.(field_list{iF}).dnum;
        break
    end
end
end

% Now concatenate data from the rest of the files
if ~isempty(data)
    for d=2:length(TimeIndex.filenames)
        file_data = load(fullfile(obj.Meta_Data.paths.mat_data,TimeIndex.filenames{d}));
        data = [data; file_data.(field_list{iF}).(sub_list{iS})];
        dnum = [dnum; file_data.(field_list{iF}).dnum];
    end
else
    data = [];
    dnum = [];
end

var_timeseries.dnum = dnum;
var_timeseries.data = data;
