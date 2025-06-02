%% Change these parameters for post-processing data
epsi_depth_array = 0:200;
fctd_depth_array = 0:300;

survey_name = '25_0413_d_fctd1_dye_survey';
MDP_file = 'MD_2025_04_12_032652.txt';

% Colorbar limits for section plots
clims.temperature = [5 25];
clims.salinity = [34.3 35.5];
clims.epsilon = [-10 -7.5];
clims.chi = [-10 -6];
clims.chi = [-10 -6];
clims.n2 = [-6 -2.7];
clims.fluor = [];

%% ----------------------------------------------------------------------
% Define directory containing raw files
raw_dir = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/ReProcessed/';

% Path to setup file
root_software='/Volumes/DEV1_HD/Users/Shared/Software_current_cruise/MOD_fish_lib/';
Meta_Data_process_file = fullfile(raw_dir,survey_name,MDP_file);

% Define parent directory into which you will create this deployment directory with raw, mat, profiles directories inside
process_dir_root = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/ReProcessed';

% Define the name of the process directory based on survey_name
process_dir = fullfile( ...
    process_dir_root, ...
    strrep(survey_name,'''','')); %This will create a directory with this name

%% Begin processing
fprintf('Processing survey %s\n', survey_name)
fprintf('-----------------------------------------------------\n')

%% Make directories
% Make the process directory
if ~exist(fullfile(process_dir,'raw'),'dir')
    mkdir(process_dir);
end

% Define the directories for epsiProcess_convert_new_raw_to_mat
dirs.raw_incoming = fullfile( ...
                    raw_dir, ...
                    strrep(survey_name,'''',''),...
                    'raw'); %This will create a directory with this name

dirs.raw_copy  = fullfile(process_dir,'raw');
dirs.mat       = fullfile(process_dir,'mat');
dirs.fctd_deployment = fullfile(process_dir,'fctd_mat');
dirs.fctd_cruise = fullfile(process_dir,'fctd_cruise'); %Redundant, just to get processing to run

% Create directories if they don't exist
if ~exist(process_dir,'dir')
    eval([ '!mkdir ' strrep(process_dir,' ','\ ')]);
end
if ~exist(dirs.raw_copy,'dir')
    eval([ '!mkdir ' strrep(dirs.raw_copy,' ','\ ')]);
end
if ~exist(dirs.mat,'dir')
    eval([ '!mkdir ' strrep(dirs.mat,' ','\ ')]);
end
if ~exist(dirs.fctd_deployment,'dir')
    eval([ '!mkdir ' strrep(dirs.fctd_deployment,' ','\ ')]);
end
if ~exist(dirs.fctd_cruise,'dir')
    eval([ '!mkdir ' strrep(dirs.fctd_cruise,' ','\ ')]);
end

%% Copy files
% Copy the first file that matches str_to_match from raw_incoming into
% raw_copy - you need to have one file there for epsi_class to read the
% configuration information
file_list_all = dir(fullfile(dirs.raw_incoming,'*.modraw'));

% Loop through files and find the ones with survey_name
idx_in_survey = false(length(file_list_all),1);
for i=1:length(file_list_all)

    % Open file
    fid = fopen(fullfile(file_list_all(i).folder,file_list_all(i).name));
    fseek(fid,0,1);
    frewind(fid);
    str = fread(fid,'*char')';
    fclose(fid);

    % Find line that has survey name
    survey_flag=contains(str,'CTD.survey');

    if survey_flag
        surveyflag_str      = str(strfind(str,'CTD.survey')+(0:100));
        surveyflag_str      = surveyflag_str(1:find(uint8(surveyflag_str)==10,1,'first'));
        surveyflag_name     = strsplit(surveyflag_str,'=');
        survey_name_in_file = surveyflag_name{2}(1:end-1);

        % Does survey name in file match the survey name we're looking for?
        if contains(survey_name_in_file,survey_name)
            idx_in_survey(i) = true;
        end
    end
end %End loop through all files

% Keep only files in survey
file_list = file_list_all(idx_in_survey);
already_copied_list = dir(fullfile(dirs.raw_copy,'*.modraw'));
already_copied_list = {already_copied_list(:).name};

for i=1:length(file_list)
    % Copy each file unless it has already been copied
    if ~any(contains(already_copied_list,file_list(i).name))
        eval(['!cp ' fullfile(file_list(i).folder,file_list(i).name) ' ' dirs.raw_copy]);
    end
end

%% Get the fish name from one of the copied files
% Open file
fid = fopen(fullfile(dirs.raw_copy,file_list(1).name));
fseek(fid,0,1);
frewind(fid);
str = fread(fid,'*char')';
fclose(fid);

fish_flag=contains(str,'CTD.fishflag=');
if fish_flag
    fishflag_str      = str(strfind(str,'CTD.fishflag=')+(0:100));
    fishflag_str      = fishflag_str(1:find(uint8(fishflag_str)==10,1,'first'));
    fishflag_name      = strsplit(fishflag_str,'=');
    fishflag_name      = fishflag_name{2}(2:end-2);
    instrument = fishflag_name;

else
    instrument = input('What fish are we using? [''epsi'',''fctd'']');

end


%%  Get the depth array based on instrument choice
switch lower(instrument)
    case 'epsi'
        depth_array = epsi_depth_array;
    case 'fctd'
        depth_array = fctd_depth_array;
end
input_struct.depth_array=depth_array;

%% Process and plot the data
ec = epsi_class(process_dir,Meta_Data_process_file);

switch lower(instrument)
    case 'epsi'
        fprintf('Reading data...\n')
        ec.f_readData;
        fprintf('Making new profiles...\n')
        ec.f_makeNewProfiles;
        fprintf('Computing turbulence...\n')
        ec.f_computeTurbulence
        fprintf('Gridding profiles...\n')
        ec.f_gridProfiles(depth_array);
        plot_epsi_sections
    case 'fctd'
        fprintf('Reading data...\n')
        epsiProcess_convert_new_raw_to_mat(dirs,ec.Meta_Data,"doFCTD");
        [FCTDall,FCTDgrid] = concatenate_and_grid_fctd(dirs.fctd_deployment);
        plot_fctd_sections
end


