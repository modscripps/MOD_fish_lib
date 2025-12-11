clear
clc

% set color of command window
cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
listeners = cmdWinDoc.getDocumentListeners;
% find text area part
jTextArea=listeners(5); %or listeners(3) or listeners (4) depending on matlab
% set color
jTextArea.setBackground(java.awt.Color(1,0.4,0.4));

% input structure
% rot_struct.raw_incoming = '/Users/gunnar/Projects/motive/cruises/cruise2/fctd-twist-counting/incoming/';
% rot_struct.mat = '/Users/gunnar/Projects/motive/cruises/cruise2/fctd-twist-counting/mat/';
% rot_struct.rot = '/Users/gunnar/Projects/motive/cruises/cruise2/fctd-twist-counting/rot/';
% rot_struct.raw_copy = '/Users/gunnar/Projects/motive/cruises/cruise2/fctd-twist-counting/raw/';
% root_software = '/Users/gunnar/Projects/matlab/toolboxes_git/FromOtherDevelopers/MOD_fish_lib/';


rot_struct.raw_incoming = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/RAW_full_cruise_twist_counter/';
rot_struct.mat = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/MAT_full_cruise_twist_counter/';
rot_struct.rot = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/ROT_full_cruise_twist_counter/';
rot_struct.raw_copy = '/Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/RAW_full_cruise_twist_counter/';

root_software = '/Volumes/DEV1_HD/Users/Shared/Software_current_cruise/MOD_fish_lib/';


Meta_Data_process_file = 'MDP_motive_2025.txt';
input_struct.raw_dir = rot_struct.raw_incoming;
Meta_Data_process_dir = fullfile(root_software, '/EPSILOMETER/Meta_Data_Process/');
input_struct.Meta_Data_process_file = fullfile(Meta_Data_process_dir, Meta_Data_process_file);

Meta_Data_process_file = input_struct.Meta_Data_process_file;
obj = epsi_class(input_struct.raw_dir, Meta_Data_process_file);

[matData] = epsiProcess_convert_new_raw_to_mat_twist_counter(rot_struct,obj.Meta_Data);


f = figure(1);
clf(f)
f.Units = 'normalized';
set(f,'position',[0.6, 0, 0.4, 0.4]);

% start counting after spool swap 2025-21-11 03:40 UTC 
rot = GV_PlotUpAccumulation(rot_struct.mat, rot_struct.rot, '20251211_030000');

% ROT_timer = timer;
% s = false;
% 
% ROT_timer.StartFcn = 'disp(''Plot ROTATION begins now!'');';
% ROT_timer.TimerFcn = [...
%     'if s, '...
%     'stop(ROT_timer); '...
%     'delete(ROT_timer); '...
%     'else, '...
%     'try,'...
%     '[matData] = epsiProcess_convert_new_raw_to_mat_twist_counter(rot_struct,obj.Meta_Data);'...
%     'catch err,'...
%     'display_error_stack(err); '...
%     'end;'...
%     'try,'...
%     'PlotUpAccumulation;'...
%     'catch err,'...
%     'display_error_stack(err); '...
%     'end;'...
%     'disp([''Done at '' datestr(now)]);'...
%     'end;'];
% ROT_timer.Period = 10*60;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
% ROT_timer.BusyMode = 'drop';
% ROT_timer.Name = 'ROT_timer';
% ROT_timer.Tag = 'ROT_timer';
% ROT_timer.StopFcn = 'clear(''rawDir'',''rawDirAway''); disp([datestr(now) '': Stopped ROT_timer'']);';
% ROT_timer.ExecutionMode = 'fixedSpacing';
% % ROT_timer.ExecutionMode = 'singleShot';
% ROT_timer.TasksToExecute = Inf;
% ROT_timer.ErrorFcn = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';
% 
% start(ROT_timer);