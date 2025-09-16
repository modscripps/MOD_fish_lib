% Before running this script...
% 
% (1) Open your yaml_file and edit

process_dir = '/Volumes/Epsidrive/COMCEPTS/DATA/25_0911_minnow3_005';
yml_file = '25_0911_minnow3_005.yml';

depth_array = 0:1:300;

%% Process and plot the data
addpath(genpath('/Volumes/Epsidrive/COMCEPTS/MOD_fish_lib'))

Meta_Data_yaml = fullfile(process_dir,yml_file);
ec = epsi_class_yaml(process_dir,Meta_Data_yaml);
%%
fprintf('Reading data...\n')
ec.f_readData;
%%
fprintf('Making new profiles...\n')
ec.f_makeNewProfiles;

fprintf('Computing turbulence...\n')
ec.f_computeTurbulence

fprintf('Gridding profiles...\n')
ec.f_gridProfiles(depth_array);
%%
plot_epsi_sections