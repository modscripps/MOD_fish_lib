% RUN_epsiAuto_process_data.m
%
% This is a wrapper script to set up input values and run
% epsiAuto_process_data.m
%
% CREATE:
% input_struct, a structure containing the following fields
%
% Required:
%  .raw_dir = path to directory where the data are streaming in
%  .process_dir = path to directory where data will be copied and where
%              subdirectories raw, mat, and FCTDmat will be created
%  .Meta_Data_process_file = path to Meta_Data_Process text file
%
% Optional:
%  .refresh_time_sec = (default 5*60), refresh period in seconds
%  .version       = (default 4), version of mod_som_read_epsi_files.m to use
% -------------------------------------------------------------------------
%
% USER INPUTS
% These will probably be the same for the whole cruise
clear input_struct
input_struct.Meta_Data_process_file = '/Volumes/DEV1_HD/Users/Shared/Software_current_cruise/MOD_fish_lib/EPSILOMETER/Meta_Data_Process/MDP_motive_2024.txt';
input_struct.refresh_time_sec =  30*60;
% input_struct.cruise_specifics = 'tfo_2024';
epsi_depth_array = 0:1000;
fctd_depth_array = 0:1100;

% Set clims
clims.temperature = [20 27];
clims.salinity = [34.3 35.5];
clims.epsilon = [-10 -7.5];
clims.chi = [-9 -6];
clims.n2 = [-6 -2.7];

% Set xlims

% List FCTD variables to grid
vars2grid_list = {'pressure','temperature','conductivity','longitude','latitude','bb','chla','fDOM','chi','chi2'};

% Realtime or Simulator mode
data_mode = 'realtime'; %'realtime' or 'simulator'

% -------------------------------------------------------------------------

%% Find the fish name and the path to store processed data
% 1. If you're running this in realtime mode, the path to processed data will
% be in the most recent Setup file.
%
% 2. If you're running this is simulated mode, the most recent Setup file is
% probably not what you want or the path to it might not exist (because
% you're running this on a personal computer. In that case, ask for the
% directory where you want to store processed data.

% First, look for the instrument name and the survey name from the most
% recent .modraw file. Create the output directory based on the survey
% name. If it doesn't exist in the file, ask the user to define it. If the
% instrument name doesn't exist in the file, ask the user to define it.

% Then look for any (up to 10) .modraw files you might have missed that
% have the same



% Get the
switch data_mode

    case 'simulator'

        % Define raw directory for simulated data
        input_struct.raw_dir = '/Volumes/DEV3_HD/Users/Shared/EPSI_PROCESSING/Simulated_Data/Realtime_RAW/raw/';

        % Path to setup
        path2setup = '/Volumes/DEV1_HD/Users/Shared/FCTD_EPSI_DATA/Simulated_Data/Setup';

        process_dir_root = '/Volumes/DEV3_HD/Users/Shared/EPSI_PROCESSING/Simulated_Data/Processed';

    case 'realtime'

        % Define raw directory for realtime data
        input_struct.raw_dir = '/Volumes/DEV3_HD/Users/Shared/EPSI_PROCESSING/Current_Cruise/Realtime_RAW/raw/';

        % Path to setup file
        root_software='/Volumes/DEV1_HD/Users/Shared/Software_current_cruise/MOD_fish_lib/';
        path2setup=fullfile(root_software,'Acquisition/fctd_epsi_acq/build/fctd_epsi/Build/Products/Debug/Setup');

        process_dir_root = '/Volumes/DEV3_HD/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed';
end

%% Read the Setup file
% Look for fish flag in Setup file
fid=fopen(path2setup,'r');
fseek(fid,0,1);
frewind(fid);
str = fread(fid,'*char')';
fclose(fid);
newSetup_flag=contains(str,'CTD.fishflag=');
if newSetup_flag
    fishflag_str      = str(strfind(str,'CTD.fishflag=')+(0:100));
    fishflag_str      = fishflag_str(1:find(uint8(fishflag_str)==10,1,'first'));
    fishflag_name      = strsplit(fishflag_str,'=');
    fishflag_name      = fishflag_name{2}(2:end-2);
    instrument_setupfile = fishflag_name;

else
    instrument_setupfile = input('What fish are we using? [''epsi'',''fctd'']');

end

% Look for survey name in Setup file
newSurvey_flag=contains(str,'CTD.survey');
if newSetup_flag
    surveyflag_str      = str(strfind(str,'CTD.survey')+(0:100));
    surveyflag_str      = surveyflag_str(1:find(uint8(surveyflag_str)==10,1,'first'));
    surveyflag_name     = strsplit(surveyflag_str,'=');
    surveylag_name      = surveyflag_name{2}(1:end-1);
    survey_name_setupfile = surveylag_name;

else
    survey_name_setupfile = input('What is the survey name? [''yyyymmdd_d##_NAME'']');

end

% Define the name of the process directory based on either what you
% found in the Setup file or what you input
input_struct.process_dir = fullfile( ...
    process_dir_root, ...
    strrep(survey_name_setupfile,'''','')); %This will create a directory with this name

%% Look for the files that match the survey name and copy them to the processing directory
% Instead of looking through the last 10 files, just look at ALL of the
% files and find the ones that match the survey name.
%
%

% Make the process directory
if ~exist(fullfile(input_struct.process_dir,'raw'),'dir')
    mkdir(input_struct.process_dir);
    mkdir(fullfile(input_struct.process_dir),'raw/')
end

% % Loop through files in incoming raw directory. Find the ones with
% % survey_name and copy those to the processing raw directory
% file_list_all = dir(fullfile(input_struct.raw_dir,'*.modraw'));
% 
% idx_in_survey = false(length(file_list_all),1);
% for i=1:length(file_list_all)
% 
%     % Open file
%     fid = fopen(fullfile(file_list_all(i).folder,file_list_all(i).name));
%     fseek(fid,0,1);
%     frewind(fid);
%     str = fread(fid,'*char')';
%     fclose(fid);
% 
%     % Find line that has survey name
%     survey_flag=contains(str,'CTD.survey');
% 
%     if survey_flag
%         surveyflag_str      = str(strfind(str,'CTD.survey')+(0:100));
%         surveyflag_str      = surveyflag_str(1:find(uint8(surveyflag_str)==10,1,'first'));
%         surveyflag_name     = strsplit(surveyflag_str,'=');
%         survey_name_in_file = surveyflag_name{2}(1:end-1);
% 
%         % Does survey name in file match the survey name we're looking for?
%         if contains(survey_name_in_file,survey_name_setupfile)
%             idx_in_survey(i) = true;
%         end
%     end
% end %End loop through all files
% 
% % Keep only files in survey
% file_list_struct = file_list_all(idx_in_survey);
% file_list = {file_list_struct(:).name};
% % Rsync the files in file_list
% raw_files_to_copy = strjoin(fullfile(input_struct.raw_dir, file_list), ' '); % Create a space-separated list of full paths
% com = sprintf('/usr/bin/rsync -av %s %s', raw_files_to_copy, fullfile(input_struct.process_dir,'raw/'));
% unix(com);


%% All options have been determined. Now get depth array, create output directory, and get ready to run the processing script
% Get the depth array based on instrument choice
switch instrument_setupfile
    case {'epsi','EPSI'}
        input_struct.depth_array = epsi_depth_array;
    case {'fctd','FCTD'}
        input_struct.depth_array = fctd_depth_array;
end

% Set command window color
set_window_color('cyan')

%% Copy Setup file onto the end of Meta_Data_Process file to make a new file with datetime stamped
% Make a new file in the process directory called 'MD' with the
% current datetime
newfile_name = fullfile(input_struct.process_dir,...
    sprintf('MD_%04.0f_%02.0f_%02.0f_%02.0f%02.0f%02.0f.txt',...
    year(now),month(now),day(now),hour(now),minute(now),second(now)));

% Open Meta_Data_Process file for reading and get content
fid1 = fopen(input_struct.Meta_Data_process_file,'r');
content1 = fread(fid1,'*char')';
fclose(fid1);

% Open Setup file for reading and get content
fid2 = fopen(path2setup,'r');
content2 = fread(fid2,'*char')';
fclose(fid2);

% Open the new file for writing amd write content of both files
fidNew = fopen(newfile_name,'w');
fwrite(fidNew,content1,'char');
fwrite(fidNew,content2,'char');
fclose(fidNew);

input_struct.Meta_Data_process_file = newfile_name;


%% Run the processing script on a timer
switch instrument_setupfile
    case {'epsi','EPSI'}
        epsiAuto_process_data
    case {'fctd','FCTD'}
        fctdAuto_process_data
end

