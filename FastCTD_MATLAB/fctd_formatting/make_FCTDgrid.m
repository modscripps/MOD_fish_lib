
function [FCTDgrid] = make_FCTDgrid(FCTDall,fctd_mat_dir)
%% Grid data and save
% NC May 2023 - grid both up and downcast data
% Get downcasts and upcasts

disp('Creating FCTDgrid')

% Grid the downcasts and upcasts
FCTDdown = FastCTD_GridData(FCTDall,'downcast');
FCTDup = FastCTD_GridData(FCTDall,'upcast');

if (~isempty(FCTDup) & ~isempty(FCTDdown))
    % Merge all the structure fields for the down and up casts
    field_list = fields(FCTDdown);
    for iF=1:length(field_list)
        % skip tGrid but do all the others
        if ~strcmp(field_list{iF},'tGrid')
            FCTDboth.(field_list{iF}) = [FCTDdown.(field_list{iF}),...
                FCTDup.(field_list{iF})];
        end
    end

    % Now sort the merged profiles by dnum
    [~,iSort] = sort(FCTDboth.time);
    field_list = fields(FCTDboth);
    for iF=1:length(field_list)
        % Don't sort depth by dnum, but do make it a Nx1 field
        if strcmp(field_list{iF},'depth')
            FCTDgrid.(field_list{iF}) = FCTDboth.(field_list{iF})(:,1);
        else
            FCTDgrid.(field_list{iF}) = FCTDboth.(field_list{iF})(:,iSort);
        end
    end


    % Save grid
    % TODO: Add a warning to set fctd_mat_dir if you want to save the data
    save(fullfile(fctd_mat_dir,'FCTDgrid'),'FCTDgrid','FCTDup','FCTDdown')
else
    disp('No profiles yet')
    FCTDgrid=[];
end %end (~isempty(FCTDup) & ~isempty(FCTDdown))
end %end make_FCTDgrid