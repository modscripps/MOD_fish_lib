function [metadata] = setup_make_metadata_from_yaml(yaml_file)
% [metadata] = setup_make_metadata_from_yaml(yaml_file)
% 
% Makes metadata structure from .yml file input 
% yaml_file: has metadata for MOD instruments

%% Read yaml file
yaml = ReadYaml(yaml_file);

%% Organize fields in metadata...


% Add mission and instrument fields
metadata.cruise_name = yaml.cruise_name;
metadata.experiment_name = yaml.experiment_name;
metadata.vehicle_name = yaml.vehicle_name;
metadata.telemetry_case_name = yaml.telemetry_case;
metadata.pressure_case_name = yaml.pressure_case;
metadata.fishflag_name = yaml.fish_flag;
metadata.deployment_name = lower(yaml.deployment_name);
metadata.start_dnum = 0; %Still need this in mod_som_read_epsi_files_v4

% Add paths fields
paths_fields = fields(yaml.paths_from_processing_machine);
for ii=1:length(paths_fields)
    metadata.paths.(paths_fields{ii}) = yaml.paths_from_processing_machine.(paths_fields{ii});
end
% Reorganize and rename some of the paths names for MOD_fish_lib
metadata.paths.data = fullfile(metadata.paths.processed_data,metadata.deployment_name);
metadata.paths.raw_incoming = metadata.paths.raw_incoming_process;
metadata.paths = rmfield(metadata.paths,'processed_data');
metadata.paths = rmfield(metadata.paths,'raw_incoming_process');
% Add path to yaml_file
metadata.paths.setup_file = yaml_file;
% Add paths to raw, mat, profiles, and figs
metadata = epsiSetup_set_epsi_paths(metadata);


%metadata.dnum_range = [datenum(yaml.dnum_min),datenum(yaml.dnum_max)];

% Add AFE fields
AFE_fields = fields(yaml.afe.channels);
for ii=1:length(AFE_fields)
    metadata.AFE.(AFE_fields{ii}) = yaml.afe.channels.(AFE_fields{ii});
end

% Add CTD fields
metadata.CTD.name = yaml.ctd.type;
metadata.CTD.SN = yaml.sn.ctd;
metadata.CTD.sample_per_record = yaml.ctd.sample_per_record;

channel_list = fields(yaml.sn);
for iF=1:numel(channel_list)
    chan = channel_list{iF};
    metadata.AFE.(chan).SN = yaml.sn.(chan);
end

% Add PROCESS fields
metadata.PROCESS.use_file_headers = yaml.use_file_headers;
metadata.PROCESS.raw_file_version = yaml.raw_files.version;
metadata.PROCESS.raw_file_suffix = yaml.raw_files.suffix;
metadata.PROCESS.channels = fields(yaml.afe.channels);
metadata.PROCESS.nb_channels = numel(metadata.PROCESS.channels);
metadata.PROCESS.latitude = yaml.latitude;
metadata.PROCESS.profile_dir = yaml.profiling_direction;

% Add GEOMETRY fields
fish = lower(metadata.fishflag_name);
metadata.GEOMETRY.alt_angle_deg = yaml.altimeter.(fish).angle_deg;
metadata.GEOMETRY.alt_dist_from_crashguard_ft = yaml.altimeter.(fish).dist_from_crashguard_ft;
metadata.GEOMETRY.alt_probe_dist_from_crashguard_in = yaml.altimeter.(fish).probe_dist_from_crashguard_in;

% Add plot_properties fields
plot_properties_fields = fields(yaml.plot_properties);
for ii=1:length(plot_properties_fields)
    metadata.plot_properties.(plot_properties_fields{ii}) = ...
        yaml.plot_properties.(plot_properties_fields{ii});
end

%% Check that CTD and AFE serial numbers have calibration files.
% Ask for new serial numbers if they don't
metadata.CTD.cal = get_CalSBE(fullfile(metadata.paths.calibrations.ctd,[metadata.CTD.SN,'.cal']));
