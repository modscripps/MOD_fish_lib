function FCTDprofiles=FCTD_FindProfiles_AAA(FCTD,datapath)
%
% FCTD_FindProfiles_AAA is a re-write of Rob's H_2016_C_FCTDreader.m
% It now uses profile finding the same way that FastCTD_GridData does, with
% the stipulation that only pairs of up-down profiles are saved.
% The code is designed to read an FCTD data structure (from San's codes),
% sort the data into up and down profiles, and save the
% profiles into up and down data without gridding.
%
% Inputs:
%  FCTD: a data structure of FCTD data, which is the output of San's
%  codes, containing at minimum 1-D vectors of time, temperature,
%  pressure, and conductivity, as well as a header structure with fields
%  seialno and gm_time
%
%  datapath: a string path to the data file where the intermediate output of
%  the processing step will be saved in, such as 
%  '/Volumes/Andriatis_T7BL/Data/TFO/TFO1/FCTD/Matching/FCTDprofiles.mat'
%
%
% Outputs:
%  FCTDprofiles: a data structure containing structures 'up' and down';
%  Each one contains cell arrays for each variable - time, pressure,
%  temperature, and conductivty. Each entry in the cell array is the
%  timeseries of data for a profile
%
%
%
% Alex Andriatis
% 2021-02-01
FCTDprofiles=[];

if ~ischar(datapath)
    error('The path for saving the FCTD profiles is not a character array');
end

reqfields = {'time','temperature','conductivity','pressure'};

if ~all(isfield(FCTD,reqfields))
    error('One or more required fields is missing from the FCTD structure');
end

if ~(isfield(FCTD,'header'))
    error('No header file present in the FCTD structure');
end

if ~(isfield(FCTD,'GPS'))
    error('No GPS file present in the FCTD structure');
end

serial=NaN(1,length(FCTD.header.gm_time));
time=NaN(1,length(FCTD.header.serialno));
for i=1:length(FCTD.header.gm_time)
    time(i)=datenum(FCTD.header.gm_time{i});
    serial(i)=str2num(FCTD.header.serialno{i});
end

FCTD.latitude=FCTD.GPS.latitude;
FCTD.longitude=FCTD.GPS.longitude;

% Make a list of auxilary fields to carry through
allfields=fieldnames(FCTD);
% Keep only those that match the time dimension
NT = length(FCTD.time);
keepfields={};
for i=1:length(allfields)
    field = allfields{i};
    if size(FCTD.(field),1)==NT
        keepfields{end+1}=field;
    end
end
    



 
 
 disp('Finding downcasts');

FCTD = FastCTD_FindCasts(FCTD,'downcast');
down = logical(FCTD.drop);

disp('Finding upcasts');
FCTD = FastCTD_FindCasts(FCTD,'upcast');
up = logical(FCTD.drop);
drop = down-up;

if sum(down)==0 || sum(up)==0
    error('No up-down profile found');
end

I_drop = 1:length(FCTD.time);
t_down = FCTD.time(down);
I_drop = I_drop(down);
I_drop = I_drop((diff(t_down)*86400>30)); %Define individual drops as being spaced by at least 30 seconds between downcasts
I_drop = [0,I_drop,length(drop)];
drops = length(I_drop); % Number of potential profiles

nprofile=0;
for i=2:drops
    fprintf('Trying %i of %i profiles \n',i-1,drops-1);
    down_start = find(drop(I_drop(i-1)+1:I_drop(i))==1,1);
    down_start=down_start+I_drop(i-1);
    
    
    if i<drops
        up_start = find(drop(down_start:I_drop(i+1))==-1,1);
    else
        up_start = find(drop(down_start:end)==-1,1);
    end
    up_start = up_start+down_start-1;
    
    down_end = find(drop(down_start:up_start)<1,1);
    down_end = down_end+down_start-1;
    
    if i<drops
        up_end = find(drop(up_start:I_drop(i+1))>-1,1);
    else
        up_end = find(drop(up_start:end)>-1,1);
    end
    up_end = up_end+up_start-1;
    
    I_down = [down_start:down_end];
    I_up = [up_start:up_end];
    
    if isempty(I_down) || isempty(I_up)
        continue
    end
    
        
    if (max(FCTD.pressure(I_down))-min(FCTD.pressure(I_down)))<10 || (max(FCTD.pressure(I_up))-min(FCTD.pressure(I_up))<10)
        continue; % Profile must be at least 10 meters deep to count
    end
    
    if max(min(FCTD.pressure(I_down)),min(FCTD.pressure(I_up)))>min(max(FCTD.pressure(I_down)),max(FCTD.pressure(I_up)))
        continue; % Pressures of profiles must overlap
    end
    
    nprofile=nprofile+1;
    for n=1:length(keepfields)
        name = keepfields{n};
        FCTDprofiles.down.(name){nprofile}=FCTD.(name)(I_down,:);
        FCTDprofiles.up.(name){nprofile}=FCTD.(name)(I_up,:);
    end
    
    tmean = mean([FCTDprofiles.down.time{nprofile};FCTDprofiles.up.time{nprofile}]);
    [~,Itime] = closest(time,tmean);
    FCTDprofiles.header.serial(nprofile)=serial(Itime);
end

if length(FCTDprofiles.up.time)~=length(FCTDprofiles.down.time)
    error('Number of up and down profiles dont match')
end

disp('Saving...');
save(datapath,'-struct','FCTDprofiles','-v7.3');
disp(['Saved FCTD profiles in ' datapath]);

