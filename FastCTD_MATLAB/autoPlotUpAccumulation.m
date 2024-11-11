clear;
clc;
cmdWinDoc=com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
listeners = cmdWinDoc.getDocumentListeners;
%find text area part
jTextArea=listeners(5); %or listeners(3) or listeners (4) depending on matlab
%set colour of command window
jTextArea.setBackground(java.awt.Color(1,0.4,0.4));

%% input structure
rot_strut.raw_incoming = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/RAW_full_cruise/';
rot_strut.mat    = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/MAT_full_cruise/';
rot_strut.rot    = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/ROT_full_cruise/';
rot_strut.raw_copy   = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/RAW_full_cruise/raw/';

root_software='/Volumes/DEV1_HD/Users/Shared/Software_current_cruise/MOD_fish_lib/';
Meta_Data_process_file = 'MDP_motive_2024.txt';
input_struct.raw_dir   = rot_strut.raw_incoming;
Meta_Data_process_dir = fullfile(root_software,'/EPSILOMETER/Meta_Data_Process/');
input_struct.Meta_Data_process_file = fullfile(Meta_Data_process_dir,Meta_Data_process_file);

Meta_Data_process_file = input_struct.Meta_Data_process_file;
obj = epsi_class(input_struct.raw_dir,Meta_Data_process_file);

%%

f = figure(1);
f.Units = 'normalized';
set(f,'position',[0 0,0.6,0.6]);
clf;

ROT_timer = timer;
s = false;

ROT_timer.StartFcn = 'disp(''Plot ROTATION begins now!'');';
ROT_timer.TimerFcn = [...
    'if s, '...
    'stop(ROT_timer); '...
    'delete(ROT_timer); '...
    'else, '...
    'try,'...
    '[matData] = epsiProcess_convert_new_raw_to_mat(rot_strut,obj.Meta_Data);'...
    'catch err,'...
    'display_error_stack(err); '...
    'end;'...
    'try,'...
    'PlotUpAccumulation;'...
    'catch err,'...
    'display_error_stack(err); '...
    'end;'...
    'disp([''Done at '' datestr(now)]);'...
    'end;'];
ROT_timer.Period = 60;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
ROT_timer.BusyMode = 'drop';
ROT_timer.Name = 'ROT_timer';
ROT_timer.Tag = 'ROT_timer';
ROT_timer.StopFcn = 'clear(''rawDir'',''rawDirAway''); disp([datestr(now) '': Stopped ROT_timer'']);';
ROT_timer.ExecutionMode = 'fixedSpacing';
% ROT_timer.ExecutionMode = 'singleShot';
ROT_timer.TasksToExecute = Inf;
ROT_timer.ErrorFcn = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';

start(ROT_timer);