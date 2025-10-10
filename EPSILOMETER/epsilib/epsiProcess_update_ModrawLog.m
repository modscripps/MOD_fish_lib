function [ModrawLog] = epsiProcess_update_ModrawLog(file_list_all,ModrawLog,Meta_Data)
% [ModrawLog] = epsiProcess_update_ModrawLog(file_list_all,ModrawLog)
%
% INPUTS
%       file_list_all - a structure of files in a directory, made by calling
%                       something like dir(fullfile(data_path,*.modraw));
%       ModrawLog - a 4-column table with headers 'File_Name','Survey_Name','Cal_File_SN','DCal_SN'
%                   If you don't have a ModrawLog already started, use
%                   epsiProcess_make_ModrawLog(file_list_all) to start a
%                   new one.
             
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
        surveyflag_str      = surveyflag_str(1:find(uint8(surveyflag_str)==10,1,'first')); %Find the first newline - ('\n') has an ASCII value of 10
        surveyflag_name     = strsplit(surveyflag_str,'=');
        survey_name_in_file = surveyflag_name{2}(1:end-1);
    else
        survey_name_in_file = '';
    end

    % Find line that starts CTD calibration data from stored calibration file
    calfile_flag = contains(str,'SERIALNO');
    if calfile_flag
        calfile_str = str(strfind(str,'SERIALNO')+(0:100));
        calfile_str = calfile_str(1:find(uint8(calfile_str)==10,1,'first'));
        calfile_str_parts = strsplit(calfile_str,'=');
        calfile_SN = strrep(calfile_str_parts{2}(1:end-1),' ','');
    else
        calfile_SN = '';
    end

    % Find line that starts CTD calibration data from instrument
    dcal_flag = contains(str,'$DCAL');
    if dcal_flag
        dcal_str = str(strfind(str,'$DCAL')+(0:100));
        dcal_str = dcal_str(1:find(uint8(dcal_str)==10,1,'first'));
        dcal_str_parts = strsplit(dcal_str,'SERIAL NO.');
        dcal_SN = strrep(dcal_str_parts{2}(1:end-2),' ','');
    else
        dcal_SN = '';
    end

    % Append row to table
    ModrawLog = [ModrawLog;
        {string(file_list_all(i).name), survey_name_in_file, calfile_SN, dcal_SN}];

    % Save modraw log
    save(fullfile(Meta_Data.paths.raw_data,'ModrawLog'))

end %End loop through all files
