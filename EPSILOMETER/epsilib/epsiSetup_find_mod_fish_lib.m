function [Meta_Data] = epsiSetup_find_mod_fish_lib(Meta_Data)


remake_paths = 0;


% If paths are listed as cells, something is wrong, remake paths.
all_dirs_single_path = ~iscell(Meta_Data.paths.process_library) && ...
                       ~iscell(Meta_Data.paths.calibrations.ctd) && ...
                       ~iscell(Meta_Data.paths.calibrations.shear) && ...
                       ~iscell(Meta_Data.paths.calibrations.fpo7);

if ~all_dirs_single_path

    remake_paths = 1;

else
    
    % If paths are single directories, not cells, but still not found on
    % current computer, remake paths for this computer.
    all_dirs_found = isdir(Meta_Data.paths.process_library) && ...
        isdir(Meta_Data.paths.calibrations.ctd) && ...
        isdir(Meta_Data.paths.calibrations.shear) && ...
        isdir(Meta_Data.paths.calibrations.fpo7);

    if ~all_dirs_found

        remake_paths = 1;
    end
end

if remake_paths
        % Find path to MOD_fish_lib
        % Split MATLAB path into entries
        p = strsplit(path, pathsep);

        % Helper: get the last folder name in each path
        lastseg = @(d) regexp(regexprep(d, [filesep '+$'], ''), '[^\\/]+$', 'match','once');

        % Find all path entries whose last folder name is MOD_fish_lib
        mod_dirs = p(cellfun(@(d) strcmpi(lastseg(d), 'MOD_fish_lib'), p));

        % If multiple matches are found, choose the shortest (closest to root)
        if numel(mod_dirs) > 1
            [~, idx] = min(cellfun(@(d) numel(strfind(d, filesep)), mod_dirs));
            mod_fish_lib = mod_dirs{idx};
            warning('Multiple MOD_fish_lib paths found. Using: %s', mod_fish_lib);
        elseif numel(mod_dirs) == 1
            mod_fish_lib = mod_dirs{1};
        else
            error('No MOD_fish_lib directory found on path.');
        end

        Meta_Data.paths.process_library = mod_fish_lib;
        Meta_Data.paths.calibrations.ctd = fullfile(mod_fish_lib,'Acquisition','SBECAL');
        Meta_Data.paths.calibrations.shear = fullfile(mod_fish_lib,'EPSILOMETER','CALIBRATION','SHEAR_PROBES');
        Meta_Data.paths.calibrations.fpo7 = fullfile(mod_fish_lib,'EPSILOMETER','CALIBRATION','FPO7');
    end

end