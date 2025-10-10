function [ModrawLog] = epsiProcess_make_ModrawLog(file_list_all,Meta_Data)
% [ModrawLog] = epsiProcess_make_ModrawLog(file_list_all)
%
% INPUTS
%       file_list_all - a structure of files in a directory, made by calling
%                       something like dir(fullfile(data_path,*.modraw));
%
% If you have an existing ModrawLog you want to append, use
% epsiProcess_update_ModrawLog.

% Preallocate an empty table with descriptive headers
ModrawLog = table('Size',[0 4], ...
                  'VariableTypes',{'string','string','string','string'}, ...
                  'VariableNames',{'File_Name','Survey_Name','Cal_File_SN','DCal_SN'});

ModrawLog = epsiProcess_update_ModrawLog(file_list_all,ModrawLog,Meta_Data);