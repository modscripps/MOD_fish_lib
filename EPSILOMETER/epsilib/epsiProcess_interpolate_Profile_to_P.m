function [GRID] = epsiProcess_interpolate_Profile_to_P(Profile,P)

% GRID.experiment = Profile.Meta_Data.experiment;
GRID.mission = Profile.Meta_Data.mission;
GRID.vehicle_name = Profile.Meta_Data.vehicle_name;
% GRID.station = Profile.Meta_Data.station;
GRID.deployment = Profile.Meta_Data.deployment;
GRID.filenames = Profile.filenames;

varList0 = {'profNum','dnum','latitude','longitude'};
varList1 = {'w','t','s','th','sgth','epsilon_final','epsi_qc'}; %epsilon_final becomes epsilon
varList2 = {'epsilon_co','chi','fom','calib_volt','calib_vel'}; 
varList3 = {'a1','a2','a3'}; 
varListqc = {'a1','a2','a3','pitch','roll'}; 

% Add varList0
for iVar=1:length(varList0)
    if isfield(Profile,varList0{iVar})
        eval(['GRID.',varList0{iVar},' =  nanmean(Profile.(varList0{iVar}));']);
    end
end

% Add depth field and interpolate
GRID.pr = P(:);
try
    GRID.z = sw_dpth(GRID.pr,mean(GRID.latitude,'omitnan'));
catch
    GRID.z = sw_dpth(GRID.pr,Profile.Meta_Data.PROCESS.latitude);
end

if ~isempty(Profile.w)

    % Add varList1
    notNan = ~isnan(Profile.pr);
    for iVar=1:length(varList1)
        if isfield(Profile,varList1{iVar})
            GRID.(varList1{iVar}) = interp1(Profile.pr(notNan),Profile.(varList1{iVar})(notNan),GRID.pr);
        end
    end

    for iVar=1:length(varList3)
        if isfield(Profile.sumPa,varList3{iVar})
            GRID.(varList3{iVar}) = interp1(Profile.pr(notNan),Profile.sumPa.(varList3{iVar})(notNan),GRID.pr);
        end
    end

    % Change epsilon to epsilon
    GRID.epsilon = GRID.epsilon_final;
    GRID = rmfield(GRID,'epsilon');

    % Add varList2
    for iVar=1:length(varList2)
        eval(['GRID.',varList2{iVar},'1 = interp1(Profile.pr(notNan),Profile.(varList2{iVar})(notNan,1),GRID.pr);']);
        eval(['GRID.',varList2{iVar},'2 = interp1(Profile.pr(notNan),Profile.(varList2{iVar})(notNan,2),GRID.pr);']);
    end

    % Acceleration - find average power spectral density of verical
    % acceleration within a frequency range.
    f = Profile.f;
    f1 = 10;
    f2 = 45;
    for iScan=1:size(Profile.Pa_g_f.a1,1)
        [PSDavg(iScan,1)] = avg_psd_in_frange(Profile.Pa_g_f.a1(iScan,:),f,f1,f2);
    end
    GRID.a1_avg = interp1(Profile.pr(notNan),PSDavg(notNan),GRID.pr);

    % Use altimeter and pressure to get bottom depth
    % Interpolate pressure to altimeter time, and use all times when altimeter < 35
    if isfield(Profile,'alt')
        if length(Profile.alt.dnum)>=20 %Only calculate bottom depth for profiles with at least 20 points
            [~,iU] = unique(Profile.ctd.time_s);
            ctdZ = interp1(Profile.ctd.time_s(iU),Profile.ctd.z(iU),Profile.alt.time_s);
            hab = Profile.alt.hab;
            hab(hab>35) = nan;
            bottom_depth = hab+ctdZ;
            % Take the average of the deepest 20 measurements (20 seconds). This step
            % is an attempt to get rid of any spurious readings that might have
            % occurred further up in the profile
            [~,idxDeep] = sort(ctdZ);
            bottom_depth_median = nanmedian(bottom_depth(idxDeep(end-19:end)));
            GRID.bottom_depth = bottom_depth_median;
        else
            GRID.bottom_depth = nan;
        end %End if Profile is at least 20 points long
    end

end %end if ~isempty(Profile.w)

