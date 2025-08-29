function [repeat,obj]=epsi_class_yaml_MetaData_doesnot_exist(data_path,obj)
disp('--- epsi_class_yaml_MetaData_doesnot_exist.m ---')

Meta_Data = obj.Meta_Data;
save(fullfile(obj.Meta_Data.paths.data,'Meta_Data'),'Meta_Data');

                    % Check that Meta_Data has everything you need. If it's
                    % missing something, set repeat=0 and checkMD=[]
                    if  isdir(obj.Meta_Data.paths.process_library) && ...
                            isdir(obj.Meta_Data.paths.data) && ...
                            isdir(obj.Meta_Data.paths.calibrations.ctd) && ...
                            isclassfield(obj.Meta_Data.paths,'raw_data') && ...
                            isclassfield(obj.Meta_Data.PROCESS,'nfft')
                        Meta_Data = obj.Meta_Data;
                        %                         save(fullfile(obj.Meta_Data.paths.data,'Meta_Data'),'Meta_Data');

                        repeat = 0; %Stop repeating. Keep this Meta_Data.
                    else
                        repeat  =  1;
                    end

