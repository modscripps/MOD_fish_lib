%% This script is for up-down matching and salinity despiking of FASTCTD data.
%
% Alex Andriatis
% 2021-01-29
%
% 
%
% December 2021: Validation of processing code and testing


ccc;
addpath(genpath('~/Documents/MATLAB')); % Path to some needed MATLAB files. 

%% Define some paths for processing

% Path where FastCTD data is stored.
datapath = '/Volumes/Andriatis_T7BL/Data/TFO/TFO1/FCTD/MAT/';

% Processing path, where the outputs of this code will go
processpath = '/Volumes/Andriatis_T7BL/Data/TFO/TFO1/FCTD/Test_Processing/';

% Softwarepaths
% Adds the processing folder for despiking
softpath_processing = '/Users/alexandriatis/UCSD/Andriatis_Thesis/Projects/TFO/Processing/FCTD/TFO1/FCTD_DehystDespikingGridding_AAA_December2021';
addpath(genpath(softpath_processing));

% Adds SAN's matlab codes, which include useful scripts for fastCTD
% processing
softpath_san = '/Users/alexandriatis/UCSD/Andriatis_Thesis/Software/SAN/MATLAB_Version_MISO_BOB';
addpath(genpath(softpath_san));

figpath = fullfile(processpath,'figures');
mkdir(figpath);

experiment = 'TFO1'; % Prefix applied to all the data files

%% 1. Merge all Fast CTD data into one big file
tic
% Make a new folder for combined data 
mkdir(fullfile(processpath,'Timeseries'));
% Name the combined data
savepath = fullfile(processpath,'Timeseries',[experiment '_FCTD_full_raw.mat']);

recalculate=0; % Switch to 1 to remake data even if it exists.
if exist(savepath,'file')==2 && ~recalculate
    FCTDfull = load(savepath);
else
    if exist(fullfile(datapath,'FastCTD_MATfile_TimeIndex.mat'),'file')
        load(fullfile(datapath,'FastCTD_MATfile_TimeIndex.mat'));
        FCTDfull=[];
        for i=1:length(FastCTD_MATfile_TimeIndex.filenames)
            fname = [FastCTD_MATfile_TimeIndex.filenames{i} '.mat'];
            load(fullfile(datapath,fname));
            FCTDfull = FastCTD_MergeFCTD(FCTDfull,FCTD);
            fprintf('Finished %i of %i \n', i,length(FastCTD_MATfile_TimeIndex.filenames));
        end
    else
        fnames=dir(fullfile(datapath,'*.mat'));
        FCTDfull=[];
        for i=1:length(fnames)
            fname = fullfile(datapath,fnames(i).name);
            load(fname);
            FCTDfull = FastCTD_MergeFCTD(FCTDfull,FCTD);
            fprintf('Finished %i of %i \n', i,length(fnames));
        end
    end
        
    save(savepath,'-struct','FCTDfull','-v7.3');
    disp(['Saved timeseries FCTD mat structure in ' savepath]);
end
toc



%% Make some quick plots to make sure the data makes sense

plotFCTD_timeseries_AAA(FCTDfull,figpath);

%% !!!!! The epsi fastctd code doesn't inclue a "header" field. 
% % Here I make a "header" with a constant serial number over the whole
% % timeseries. 
% FCTDfull.header.gm_time{1}=datestr(mean(FCTDfull.time));
% FCTDfull.header.serialno{1}='1';


%% Step 1: Convert the FastCTD timeseries into pairs of profiles to be matched. 
% Replaces Rob's H_2016_C_FCTDreader.m
tic
mkdir(fullfile(processpath,'Profiles'));
savepath = fullfile(processpath,'Profiles',[experiment '_FCTD_profiles.mat']);
recalculate=0;
if exist(savepath,'file')==2 && ~recalculate
    FCTDprofiles = load(savepath);
else
    FCTDprofiles=FCTD_FindProfiles_AAA(FCTDfull,savepath);
end
toc


%% Plot the whole timeseries profile by profile

plotFCTD_profiles_AAA(FCTDprofiles,figpath)

%% Step 2: Clean up the profiles by applying some velocity and despiking criteria.
% This code introduces a new 'qc' variable into the header
% qc=0 is a good profile, which might be trimmed by this code
% qc=1 is a bad profile based on the cleaning criteria, which is left as-is
% but excluded from further processing


tic
mkdir(fullfile(processpath,'Profiles'));
savepath = fullfile(processpath,'Profiles',[experiment '_FCTDprofiles_qc.mat']);
recalculate=0;
if exist(savepath,'file')==2 && ~recalculate
    FCTDprofiles_qc = load(savepath);
else
    FCTDprofiles_qc=FCTD_CleanProfiles_AAA(FCTDprofiles,savepath);
end
toc

%% Remove profiles by hand if necessary
% qc=2 is a bad profile that is removed by hand
badprofiles=[];
FCTDprofiles_qc.header.qc(badprofiles)=2;


%% Step 3: Correct Hysteresis

tic
mkdir(fullfile(processpath,'Matching'));
savepath = fullfile(processpath,'Matching',[experiment '_FCTD_deHysteresis.mat']);
parampath = [experiment '_Parameters_deHysteresis.mat'];



recalculate=0; % Turn on for testing
if exist(savepath,'file')==2 && ~recalculate
    FCTDdehysteresis = load(savepath);
    Parameters = load(parampath);
else
    if exist(parampath,'file')
        Parameters = load(parampath);
        [FCTDdehysteresis,Parameters]=FCTD_deHysteresis_AAA(FCTDprofiles_qc,savepath,'parameters',Parameters,'figpath',figpath,'plotevery',1000);
    else
        [FCTDdehysteresis,Parameters]=FCTD_deHysteresis_AAA(FCTDprofiles_qc,savepath,'figpath',figpath,'plotevery',1000);
        save(parampath,'-struct','Parameters','-v7.3');
    end
end
toc


%% Step 4: Do repsonse-matching (salinity despiking)


tic
savepath = fullfile(processpath,'Matching',[experiment '_FCTD_responsematched.mat']);
recalculate=1;
if exist(savepath,'file')==2 && ~recalculate
    FCTD_responsematched = load(savepath);
else
    FCTD_responsematched=FCTD_responsematch_AAA(FCTDprofiles_qc,FCTDdehysteresis,char(savepath),'figpath',figpath,'plotevery',4000,'depthrange',[200 300]);
end
toc

%% Salinity Check!

tic
recalculate=1;
if recalculate
    FCTD_SalinityCheck_AAA(FCTDprofiles_qc,FCTDdehysteresis,FCTD_responsematched,'figpath',figpath,'plotevery',1000);
end
toc


%% Step 5: Check the result of the response matching code

tic
recalculate=1;
if recalculate
    FCTD_responsecheck_AAA(FCTDprofiles_qc,FCTD_responsematched,'figpath',figpath);
end
    

%% Step 6: Make a final corrected dataset

tic
mkdir(fullfile(processpath,'Corrected'));
savepath = fullfile(processpath,'Corrected',[experiment '_FCTD_corrected.mat']);
recalculate=1;
if exist(savepath,'file')==2 && ~recalculate
    FCTDcorrected = load(savepath);
else
    FCTDcorrected=FCTD_CorrectedProfiles_AAA(FCTDprofiles_qc,FCTD_responsematched,savepath,'parameters',Parameters,'figpath',figpath,'plotevery',1000);
end
toc


%% For some processing, need to unwrap the profiles back into a big timeseries
tic
savepath = fullfile(processpath,'Corrected',[experiment '_FCTD_corrected_timeseries.mat']);
recalculate=1;
if exist(savepath,'file')==2 && ~recalculate
    FCTD_corrected_timeseries = load(savepath);
else
    FCTD_corrected_timeseries = FCTD_ProfilestoTimeseries_AAA(FCTDcorrected,savepath);
end
toc


%% Before gridding, add lat and lon information

isgps = isfield(FCTDcorrected.up,'latitude') & isfield(FCTDcorrected.up,'longitude');
gpspath = ['/Volumes/Andriatis_T7BL/Data/TFO/TFO1/AssetTracking/ASSET_RECORD/ALL_STRUCT/SallyRide.mat'];
% If GPS data doesn't exist, leave gpspath blank, Gridding will use defualt
% values of 35 lat and 0 lon
savepath = fullfile(processpath,'Corrected',[experiment '_FCTD_corrected_gps.mat']);
tic
if exist(gpspath,'file') && ~isgps
    GPS = load(gpspath);
    FCTDcorrected = FCTD_AddPosition_AAA(FCTDcorrected,GPS,savepath);
   
end
toc


%% Step 6: Grid the processed data onto a uniform pressure grid

% Update this code to reflect the new GPS structure and the way I input
% data

% Pressure step for gridding
dz=0.25;

% Path to gsw_matlab
swpath = '/Users/aandriatis/Documents/MATLAB/gsw_matlab_v3_06_11';
addpath(genpath(swpath));

tic
mkdir(fullfile(processpath,'Gridded'));
savepath = fullfile(processpath,'Gridded',[experiment '_FCTD_gridded_raw.mat']);
recalculate=1;
if exist(savepath,'file')==2 && ~recalculate
    FCTD_gridded = load(savepath);
else
    FCTD_gridded = FCTD_grid_AAA(FCTDcorrected,savepath,'dz',dz);
end
toc


%% Step 7: Combine up and down grids

% Update this code

tic
savepath = fullfile(processpath,'Gridded',[experiment '_FCTD_gridded_combined_raw.mat']);
recalculate=1;
if exist(savepath,'file')==2 && ~recalculate
    FCTD_gridded_both = load(savepath);
else
    FCTD_gridded_both = FCTD_combine_grid_AAA(FCTD_gridded,savepath);
end
toc

%% Quick test plot
figure
contourf(mean(FCTD_gridded_both.time,1,'omitnan'),FCTD_gridded_both.depth,FCTD_gridded_both.sigma0)