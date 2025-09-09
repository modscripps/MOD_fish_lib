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

% If raw_incoming path includes '---deployment_name---', replace it with
% what's in meta_data.deployment_name
idx = strfind(metadata.paths.raw_incoming,'---');
if numel(idx)==2
    new_path = fullfile(metadata.paths.raw_incoming(1:idx(1)-1),...
                        yaml.deployment_name,...
                        metadata.paths.raw_incoming(idx+4:end));
    metadata.paths.raw_incoming = new_path;
end


% Add path to yaml_file
metadata.paths.setup_file = yaml_file;
% Add paths to raw, mat, profiles, and figs
metadata = epsiSetup_set_epsi_paths(metadata);

%metadata.dnum_range = [datenum(yaml.dnum_min),datenum(yaml.dnum_max)];

% Add AFE fields
AFE_fields = fields(yaml.afe);
for ii=1:length(AFE_fields)
    if ~strcmp(AFE_fields{ii},'channels')
        metadata.AFE.(AFE_fields{ii}) = yaml.afe.(AFE_fields{ii});
    end
end

AFE_fields_2 = fields(yaml.afe.channels);
for ii=1:length(AFE_fields_2)
    metadata.AFE.(AFE_fields_2{ii}) = yaml.afe.channels.(AFE_fields_2{ii});
end

% Add CTD fields
metadata.CTD.name = yaml.ctd.type;
metadata.CTD.SN = yaml.sn.ctd;
metadata.CTD.sample_per_record = yaml.ctd.sample_per_record;

% Add probe serial numbers
channel_list = fields(yaml.sn);
for iF=1:numel(channel_list)
    chan = channel_list{iF};
    metadata.AFE.(chan).SN = yaml.sn.(chan);
end

% Get calibrations for shear probe channels
metadata = mod_som_get_shear_probe_calibration_v2(metadata);

% Add more PROCESS fields
metadata.PROCESS.use_file_headers = yaml.use_file_headers;
metadata.PROCESS.ctd_in_modraw = yaml.ctd_in_modraw;
metadata.PROCESS.raw_file_version = yaml.raw_files.version;
metadata.PROCESS.raw_file_suffix = yaml.raw_files.suffix;
metadata.PROCESS.channels = fields(yaml.afe.channels);
metadata.PROCESS.nb_channels = numel(metadata.PROCESS.channels);
metadata.PROCESS.latitude = yaml.latitude;
metadata.PROCESS.profile_dir = yaml.profiling_direction;
metadata.PROCESS.Fs_epsi = yaml.afe.sample_rate;
metadata.PROCESS.Fs_ctd = yaml.ctd.sample_rate;
metadata.PROCESS.nfft = yaml.spectral.nfft;
metadata.PROCESS.nfftc = yaml.spectral.nfftc;
metadata.PROCESS.dof = yaml.spectral.dof;
metadata.PROCESS.dof_coherence = yaml.spectral.dof_coherence;
metadata.PROCESS.dz = yaml.spectral.dz;
metadata.PROCESS.fc1 = yaml.spectral.frequency_cutoff_1;
metadata.PROCESS.fc2 = yaml.spectral.frequency_cutoff_2;
metadata.PROCESS.ctd_fc = yaml.spectral.ctd_frequency_cutoff;
metadata.PROCESS.movmean_window_time = yaml.spectral.movmean_window_time;
metadata.PROCESS.adjustTemp = yaml.spectral.adjust_temp;

metadata.PROCESS.Prmin = @(x)min(x)+0.2*range(x);
metadata.PROCESS.Prmax = @(x)min(x)+0.8*range(x);

% Make 'timeseries' from 'channels'. (Add _g to acceleration channels and
for n=1:numel(metadata.PROCESS.channels)
    if contains(metadata.PROCESS.channels(n),'a')
        metadata.PROCESS.timeseries{n} = sprintf('%s_g',metadata.PROCESS.channels{n});
    elseif contains(metadata.PROCESS.channels(n),{'s','t'})
        metadata.PROCESS.timeseries{n} = sprintf('%s_volt',metadata.PROCESS.channels{n});
    elseif contains(metadata.PROCESS.channels(n),'c')
        metadata.PROCESS.timeseries{n} = sprintf('%s_count',metadata.PROCESS.channels{n});
    elseif contains(metadata.PROCESS.channels(n),'f')
        metadata.PROCESS.timeseries{n} = sprintf('%s_count',metadata.PROCESS.channels{n});
    end
end

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

%% Get CTD calibration data
if metadata.PROCESS.ctd_in_modraw
    metadata.CTD.cal = get_CalSBE(fullfile(metadata.paths.calibrations.ctd,[metadata.CTD.SN,'.cal']));
elseif ~metadata.PROCESS.ctd_in_modraw
    metadata.CTD.cal = get_CalSBE_nan;
end

%% Save Meta_Data
Meta_Data = metadata;
save(fullfile(Meta_Data.paths.data,'Meta_Data'),'Meta_Data');
