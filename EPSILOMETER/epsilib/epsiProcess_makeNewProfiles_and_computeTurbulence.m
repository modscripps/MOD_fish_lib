function obj = epsiProcess_makeNewProfiles_and_computeTurbulence(obj,varargin)
% obj = f_makeNewProfiles_and_computeTurbulence(obj,varargin)
%
% - Divides data into profiles
% - Computes turbulence variables
% - Option to grid data
%
% INPUTS:
%   obj = epsi_class object
%
% OPTIONAL INPUTS:
%   'grid',z - 'grid' flag to grid some of the profile variables onto a
%               standard pressure grid
%            - z = the depth array to use
%
% If there are optional arguments they are 'grid' and P, the pressure array
% to grid onto
disp('--- epsiProcess_makeNewProfiles_and_computeTurbulence.m ---')


if nargin>1
    makeGrid = true;
    z = varargin{1}{2};
else
    makeGrid = false;
end

% % Pick out profile drops from CTD pressure timeseries
% obj.f_getProfileIndices;
% if Profile.Meta_Data.AFE.t1.cal==0 && Profile.Meta_Data.AFE.t2.cal==0
%     warning(sprintf(['\n !!!!!!!!! \n'...
%                     'The calibration value (dTdV) for both temperature probes is 0.\n',...
%                     'Run ec.f_calibrateTemperature or process_calibrate_dTdV.m before continuing.\n'...
%                     ' !!!!!!!!! \n',...
%                     'Hit any key to continue or ctrl-C to stop processing.']))
%     pause
% end

% Load PressureTimeseries
load(fullfile(obj.Meta_Data.paths.mat_data,'PressureTimeseries.mat'),'PressureTimeseries')

%ALB VErsion control issue. The ESPIlib on  DEV3 has a slightly
%different EPSIlib version
if~isfield(PressureTimeseries,'startprof')
if(~isfield(PressureTimeseries,'startprof') && isfield(PressureTimeseries,'startdown'))
    PressureTimeseries.startprof=PressureTimeseries.startdown;
    PressureTimeseries.endprof=PressureTimeseries.startup;
else

    error('PressureTimeseries has the wrong format')
end
end
% Look for the current list of profiles
profList = dir(fullfile(obj.Meta_Data.paths.profiles,'Profile*.mat'));
profNumChar = cell2mat(cellfun(@(C) C(8:11),{profList(:).name},'uniformoutput',0).');
if ~isempty(profNumChar)
    lastProfNum = str2double(profNumChar(end,:));
    % Load the last profile
    lastProf = load(fullfile(obj.Meta_Data.paths.profiles,sprintf('Profile%04.f',lastProfNum)));

    % Does the last profile have all its data? Or was more collected in the
    % last batch of files?
    lastProf_maxTime = nanmax(lastProf.Profile.ctd.dnum);

    % If there is more data for the last profile in the PressureTimeseries,
    % rerun the last profile. In the next step, we'll run everything after the
    % lastProfNum, so subtract one from that value.
    if PressureTimeseries.dnum(PressureTimeseries.endprof(lastProfNum))>lastProf_maxTime
        lastProfNum = lastProfNum-1;
    end
elseif isempty(profNumChar)
    lastProfNum = 0;
end


% Loop through the profile indices in PressureTimeseries. Process the new
% ones
for iProf=1:length(PressureTimeseries.startprof)
    if iProf>lastProfNum
        profIdx = PressureTimeseries.startprof(iProf):PressureTimeseries.endprof(iProf);
        tMin = PressureTimeseries.dnum(PressureTimeseries.startprof(iProf));
        tMax = PressureTimeseries.dnum(profIdx(end));
        fprintf('Building Profile%03.0f of %03.0f\n',iProf,length(PressureTimeseries.startprof))

        Profile = obj.f_cropTimeseries(tMin,tMax);
        Profile.profNum = iProf;

        %ALB I do not know yet how to deal with dTdV with the new design. 
        if Profile.Meta_Data.AFE.t1.cal==0 || Profile.Meta_Data.AFE.t2.cal==0
            % warning(sprintf(['\n !!!!!!!!! \n'...
            %     'The calibration value (dTdV) for both temperature probes is 0.\n',...
            %     'Run ec.f_calibrateTemperature or process_calibrate_dTdV.m before continuing.\n'...
            %     ' !!!!!!!!! \n',...
            %     'Hit any key to continue or ctrl-C to stop processing.']))
            % warning(sprintf(['\n !!!!!!!!! \n'...
            %     'The calibration value (dTdV) for one or both temperature probes is 0.\n',...
            %     'Running ec.f_calibrateTemperature before continuing.\n'...
            %     ' !!!!!!!!! \n']))


            % Profile.Meta_Data=mod_epsi_temperature_spectra_v4(Profile.Meta_Data,Profile,1,1);
            Profile.Meta_Data=mod_epsi_linear_calibration_FP07(Profile,1);

            % pause
        end





        Profile = obj.f_computeTurbulence(Profile);
        % Sort Profile by standard field order
        if ~isempty(Profile)
            Profile = sort_profile(Profile);
            
            % % Save new profile
            % saveName = fullfile(obj.Meta_Data.paths.profiles,sprintf('Profile%04i',iProf));
            % save(saveName,'Profile');
        end
        clear Profile
    end
end
end %end f_makeNewProfiles_and_computeTurbulence
