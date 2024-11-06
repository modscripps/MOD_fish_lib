function Meta_Data= read_MetaProcess(Meta_Data,filename)
% NC 10/9/24 - Added ability to read Setup file info copied to the end of a
% Meta_Data_Process file. I started copying the Setup file to the
% Meta_Data_Process file in RUN_Auto_process_data.m.
%

fid=fopen(filename);
str_meta_data_process=textscan(fid,'%s','Delimiter','\n');
str_meta_data_process=[str_meta_data_process{:}];
fclose(fid);
for f=1:length(str_meta_data_process)

    % For the first 30 lines, put data into Meta_Data.PROCESS
    if f<30 && ~isempty(str_meta_data_process{f}) && isempty(strfind(str_meta_data_process{f},'%'))
        eval(['Meta_Data.PROCESS.' str_meta_data_process{f}])

    elseif f>30 %After the first 30 line, start looking for data from Setup
        if strncmp(str_meta_data_process{f},'CTD.experiment',14)
            Meta_Data.experiment_name = str_meta_data_process{f}(16:end);
        elseif strncmp(str_meta_data_process{f},'CTD.cruise',10)
            Meta_Data.cruise_name = str_meta_data_process{f}(12:end);
        elseif strncmp(str_meta_data_process{f},'CTD.vehicle',11)
            Meta_Data.vehicle_name = str_meta_data_process{f}(13:end);
        elseif strncmp(str_meta_data_process{f},'CTD.fish_pc',11)
            Meta_Data.pressure_case_name = str_meta_data_process{f}(13:end);
        elseif strncmp(str_meta_data_process{f},'CTD.fishflag',12)
            Meta_Data.fishflag_name = str_meta_data_process{f}(14:end);
        elseif strncmp(str_meta_data_process{f},'CTD.ch1_sn',10)
            Meta_Data.AFE.t1.SN = str_meta_data_process{f}(12:end-1);
        elseif strncmp(str_meta_data_process{f},'CTD.ch2_sn',10)
            Meta_Data.AFE.t2.SN = str_meta_data_process{f}(12:end-1);
        elseif strncmp(str_meta_data_process{f},'CTD.ch3_sn',10)
            Meta_Data.AFE.s1.SN = str_meta_data_process{f}(12:end-1);
        elseif strncmp(str_meta_data_process{f},'CTD.ch4_sn',10)
            Meta_Data.AFE.s2.SN = str_meta_data_process{f}(12:end-1);
        end
    end
end

Meta_Data.PROCESS.filename = filename;

if isfield(Meta_Data.paths,'data')
    save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data')
end