% ALB matDataDir is the FCTDmat folder of the current deployment/section
% Be carefull to update last_file_idx and rot_count_offset
%
% rot_count_offset is the last rot_count_offset obtained at the previous
% deployment
%
% last_file_idx is the
%
%%

matDataDir = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/MAT_full_cruise_twist_counter/'; %Directory where FCTD*.mat are stored
rotDataDir = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/ROT_full_cruise_twist_counter/'; %Directory where rotation data for each FCTD*.mat file will be stored
matFiles = dir([matDataDir 'EPSI*.mat']);
rotFiles = dir([rotDataDir 'EPSI*.mat']);

% -----------------------------------------------------------------------------

%Redo the most recent rotFile in case it wasn't complete, and continue through the end of the matFiles
for iFile = numel(rotFiles)-1:numel(matFiles)
    if iFile<1
        iFile=1;
    end

    % Load vnav and ctd data
    load([matDataDir '/' matFiles(iFile).name],'vnav','ctd');

    % Convert mat data file to rotation accumulation file
    if exist('vnav','var') && exist('ctd','var') && ...
            isstruct(vnav) && ~isempty(vnav.time_s)

        diff_not_neg = [0;diff(vnav.dnum)]>0;
        keep = ~isnan(vnav.dnum) & ~isinf(vnav.dnum) & diff_not_neg;
        rot.compass      = vnav.compass(keep,:);
        rot.gyro         = vnav.gyro(keep,:);
        rot.acceleration = vnav.acceleration(keep,:)./9.81;
        rot.time_s       = vnav.time_s(keep);
        rot.dnum         = vnav.dnum(keep);
        if ~isempty(ctd)
            rot.pressure     = interp1(ctd.dnum,ctd.P,rot.dnum);
        else
            rot.pressure     = rot.dnum.*nan;
        end

        % Add a variable to keep track of the file number
        rot.file_num = ones(size(rot.pressure))*iFile;

        % Define rotation matrix
        Rot_Mat = @(p,t,s)[ cos(t)*cos(s), -cos(p)*sin(s) + sin(p)*sin(t)*cos(s),  sin(p)*sin(s) + cos(p)*sin(t)*cos(s);
            cos(t)*sin(s),  cos(p)*cos(s) + sin(p)*sin(t)*sin(s), -sin(p)*cos(s) + cos(p)*sin(t)*sin(s);
            -sin(t),         sin(p)*cos(t),                         cos(p)*cos(t)];

        % Prepare empty arrays for rotated data
        rotated_comp = nan(size(rot.compass));
        rotated_gyro = nan(size(rot.gyro));
        rotated_acce = nan(size(rot.acceleration));

        % Find the number of data points
        time = rot.dnum;
        pts = numel(time);

        % For each point, rotate acceleration, compass,and gyro data
        for k = 1:pts
            [Phi, Theta, Psi, Rot_mat] = SN_RotateToZAxis(rot.acceleration(k,:));
            rotated_comp(k,:) = Rot_mat*(rot.compass(k,:).');
            rotated_gyro(k,:) = Rot_mat*(rot.gyro(k,:).');
            rotated_acce(k,:) = Rot_mat*(rot.acceleration(k,:).');
        end

        % Collect variables to save
        pressure     = rot.pressure;
        acceleration = rot.acceleration;
        gyro         = rotated_gyro;
        acce         = rotated_acce;
        file_num     = rot.file_num;

        % Normalize compass data and find the magnitude
        comp_mag = repmat(sqrt(sum(rotated_comp(:,1:2).*rotated_comp(:,1:2),2)),[1 3]);
        compass_norm = rotated_comp./comp_mag;
        compass = compass_norm(:,1)+1i*compass_norm(:,2);

        % The total rotation counts from the accelerometer is the phase of
        % the normalized compass vector.
        tot_rot_acc = phase(compass);

        % Find time interval, dt
        dt = diff(time)*24*3600;
        dt = [dt; nanmedian(dt)];

        % The total rotation counts from the gyro is its cummulative sum
        tot_rot_gyro(:,1) = cumsum(gyro(:,1).*dt);
        tot_rot_gyro(:,2) = cumsum(gyro(:,2).*dt);
        tot_rot_gyro(:,3) = cumsum(gyro(:,3).*dt);

        % Save rotation data
        save(fullfile(rotDataDir,matFiles(iFile).name),'compass','tot_rot_gyro','tot_rot_acc','time','gyro','acce','pressure','file_num');

        % Clear processed data
        clear compass tot_rot_gyro tot_rot_acc time gyro acce pressure file_num
    end
    clear FCTD;
end

%% Load rotation data that you saved in the cell above

rotFiles = dir([rotDataDir 'EPSI*.mat']);

if exist([rotDataDir 'Latest_rot_acc_count.mat'],'file')
    load([rotDataDir 'Latest_rot_acc_count.mat'], ...
        'tot_rot_gyro', ...
        'tot_rot_acc',...
        'tot_time',...
        'tot_pressure',...
        'last_file_idx',...
        'tot_file_num');
    nonan_tot_rot_gyro=tot_rot_gyro(~isnan(tot_rot_gyro));
    rot_count_offset=-nonan_tot_rot_gyro(end)/pi/2;

else
    tot_rot_gyro     = [];
    tot_rot_acc      = [];
    last_file_idx    = 1;
    tot_file_num     = [];
    tot_time         = [];
    tot_pressure     = [];
    rot_count_offset = 0;
end


for iFile = last_file_idx:numel(rotFiles)

    % Load rotation file
    rot = load([rotDataDir '/' rotFiles(iFile).name]);

    if iFile==last_file_idx

        if ~isempty(tot_file_num)
            % For the first file in the list find the indices of file_num < last_file_idx
            before_iFile = find(tot_file_num<last_file_idx);

            %If this is empty, you're at the beginning of the file
            if isempty(before_iFile)
                tot_rot_gyro = rot.tot_rot_gyro;
                tot_rot_acc  = rot.tot_rot_acc;
                tot_pressure = rot.pressure;
                tot_time = rot.time;
                tot_file_num = rot.file_num;
            elseif ~isempty(before_iFile)
                % Concatenate gyro and acc data. The rot structure has the full
                % contents of the latest file so start with everything before the last
                % file. Add the new stuff to the last value of the previous stuff
                tot_rot_gyro = [tot_rot_gyro(before_iFile,:); tot_rot_gyro(before_iFile(end),:)+rot.tot_rot_gyro];
                tot_rot_acc  = [tot_rot_acc(before_iFile); rot.tot_rot_acc(end)+rot.tot_rot_acc];
                tot_pressure = [tot_pressure(before_iFile); rot.pressure];
                tot_time =     [tot_time(before_iFile); rot.time];
                tot_file_num = [tot_file_num(before_iFile); rot.file_num];
            end

        elseif isempty(tot_file_num)
            % If file_num is empty, you're at the beginning of the
            % deployment
                tot_rot_gyro = rot.tot_rot_gyro;
                tot_rot_acc  = rot.tot_rot_acc;
                tot_pressure = rot.pressure;
                tot_time = rot.time;
                tot_file_num = rot.file_num;
        end

    elseif iFile>last_file_idx
        % For everything after, concatenate to what you had before
        tot_rot_gyro = [tot_rot_gyro; tot_rot_gyro(end,:)+rot.tot_rot_gyro];
        tot_rot_acc  = [tot_rot_acc; tot_rot_acc(end,:)+rot.tot_rot_acc];
        tot_pressure = [tot_pressure; rot.pressure];
        tot_time =     [tot_time; rot.time];
        tot_file_num = [tot_file_num; rot.file_num];
    end

    % Save the index of the last file you added
    last_file_idx = iFile;

    % Save the length of the last file you added
    last_file_length = numel(rot.time);

    clear rot;
end

save([rotDataDir 'Latest_rot_acc_count.mat'], ...
    'tot_rot_gyro', ...
    'tot_rot_acc',...
    'tot_time',...
    'tot_pressure',...
    'last_file_idx',...
    'last_file_length',...
    'tot_file_num');

%% Plot data
if ~isempty(tot_time)

    % Calculate fall rate, dpdt
    dt = diff(tot_time)*24*3600;
    dt = [dt(1); dt];
    dp = diff(tot_pressure);
    dp = [dp(1); dp];
    fall_rate = dp(:)./dt(:);

    % Clear figure
    clf;

    % Calculate the gyro value to plot
    gyro_value = tot_rot_gyro(:,3)/pi/2;
    gyro_value = gyro_value(~isnan(gyro_value));

    % Plot total gyro counts
    h(1) = plot(datetime(tot_time,'ConvertFrom','datenum'), gyro_value,'o','linewidth',2,'color','b','DisplayName','Rot by gyro [dn]');
    hold on
    if sum(fall_rate<0.001)>10
        h(2) = plot(datetime(tot_time(fall_rate<0.001),'ConvertFrom','datenum'), gyro_value(fall_rate<0.001),'x','linewidth',2,'color','c','DisplayName','Rot by gyro [up]');
    end
    set(gca,'XTickLabelRotation',45)
    hold off;
    grid on;
    xlabel('time');
    ylabel('Number of rotations');

    % Add title with rotation count
    title(['MODfish: Rotation count = ' num2str(round(gyro_value(end))) '  _ _ _ ']);

    % Add legend
    str_legend={'Rot by acc [dn]','Rot by acc [up]','Rot by gyro [dn]','Rot by gyro [up]'};
    hl = legend('Location','NorthWest');
    set(hl,'Fontsize',20);
    set(gca,'Fontsize',20);
    set(gca,'XLim',datetime([tot_time(end)-12/24 tot_time(end)],'ConvertFrom','datenum'))

    % Print last rotation count
    fprintf("Last rotation count %i\r\n",round(gyro_value(end)))

    % "Neutral is negative"
    annotation(gcf,'textbox',...
    [0.274088541666665 0.165666666666666 0.437174479166667 0.0462962962962963],...
    'String',{'"Neutral is Negative" - Set fin angle to neutral to make the rotation count go down.'});

else
    % If there is no new data, print previous rotation count
    fprintf("No new data\n")
    fprintf("Last rotation count %i\r\n",round(gyro_value(end)))
end


clear