function [FCTDall] = make_FCTDall_L0(fctd_mat_dir)
% This function concatenates all the individual FCTD .mat files into a one
% big FCTDall structure. 
disp(fprintf('Creating FCTDall_L0 from .mat files in %s',fctd_mat_dir))

% Get FCTD files in current directory
f = dir(fullfile(fctd_mat_dir,'*.mat'));
file_list = {f.name};
% Get rid of any concatenated or gridded files - they have 'all' or 'grid' in the name
is_all = contains(file_list,'all');
file_list(is_all) = [];
is_grid = contains(file_list,'grid');
file_list(is_grid) = [];

% Make a list of datenums for each file
for iF=1:length(file_list)
    yyyy = str2num(['20' file_list{iF}(5:6)]);
    mm = str2num(file_list{iF}(8:9));
    dd = str2num(file_list{iF}(11:12));
    HH = str2num(file_list{iF}(14:15));
    MM = str2num(file_list{iF}(16:17));
    SS = str2num(file_list{iF}(18:19));
    file_dnums(iF) = datenum(yyyy,mm,dd,HH,MM,SS);
end

% Look for FCTDall so you can just add the new ones
if exist(fullfile(fctd_mat_dir,'FCTDall.mat'),'file')==2 %If FCTDall already exists
    % Load FCTDall
    load(fullfile(fctd_mat_dir,'FCTDall'))

    % What's the last datenum?
    last_dnum = FCTDall.time(end);

    % Keep any file newer than 10 minutes before last_dnum
    idx = file_dnums > last_dnum - days(minutes(10));
    file_list = file_list(idx);

else % If FCTDall doesn't already exist
    % Load the first file and add microconductivity data
    load(fullfile(fctd_mat_dir,file_list{1}),'FCTD');
    disp(file_list{1})


    FCTDall = FCTD;
end

% Load the rest of the files and merge into FCTD all
for iF=2:length(file_list)
    fprintf(sprintf('Adding file %3.0f of %3.0f to FCTDall \n',iF,length(file_list)));
    disp(file_list{iF})

    % Load file
    load(fullfile(fctd_mat_dir,file_list{iF}),'FCTD');

    % Concatenate new file into FCTDall
    FCTDall = FastCTD_MergeFCTD(FCTDall,FCTD);
end

% You might have reprocessed some of the same data. Find unique dnums and
% only keep those data
nTime = max(size(FCTDall.time));
[~,iU] = unique(FCTDall.time);

% Get field list and keep only unique dnum data for all fields
field_names = fields(FCTDall);
for f=1:length(field_names)
    if any(size(FCTDall.(field_names{f}))==nTime)
        FCTDall.(field_names{f}) = FCTDall.(field_names{f})(iU,:);
    end
end

% Add depth
if isfield(FCTDall,'latitude')
    LAT = nanmean(FCTDall.latitude);
else
    LAT = 30;
end
if isnan(LAT)
    LAT = 30;
end
FCTDall.depth = sw_dpth(FCTDall.pressure,LAT);

% Output data and save concatenated file
save(fullfile(fctd_mat_dir,'FCTDall_L0'),'FCTDall','-v7.3')

end %end make_FCTDall