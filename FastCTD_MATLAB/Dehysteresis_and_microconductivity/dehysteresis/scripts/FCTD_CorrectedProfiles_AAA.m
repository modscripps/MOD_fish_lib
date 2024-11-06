function FCTDcorrected = FCTD_CorrectedProfiles_AAA(FCTDprofiles,FCTDmatched,datapath,varargin)
%
% FCTD_responsematch_AAA is a re-write of Rob's
% K_2016C_FCTD_Response_Fleat.m
%
% It now works with a profile structure saved as cell arrays.
%
% Inputs:
%   FCTDprofiles: A strucure of profiles which is the output of
%   FCTD_FindProfiles_qc_AAA.
%
%   FCTDcorrected: A structure of hysteresis-corrected profiles which is
%   the output of FCTD_deHysteresis_AAA.
%  
%   datapath: A file path where the intermediate processing step is saved
%
%   Optional:
%
%       figpath: The string path to the folder where diagnostic figures from this
%           code will be saved. If not supplied, figures will not be
%           saved.
%               Example: 'figpath,'./figures/'
%
%       plotevery: A number, plots the hysteresis correction for every
%           plotevery profiles. To plot only the first profile, set
%           plotevery equal to the number of profiles. By defualt,
%           plotevery is set to zero and plots will not be made.
%               Example: 'plotevery',10
%Outputs:
%   FCTDcorrected: A strucutre of profiles after response matching, for
%   both raw and hysteresis-corrected profiles.
%
% Examples:
%
%   FCTDcorrected=FCTD_responsematch_AAA(FCTDprofiles_qc,FCTDdehyst,'./Matching/FCTD_responsematched.mat')
%   FCTDcorrected=FCTD_responsematch_AAA(FCTDprofiles_qc,FCTDdehyst,'./Matching/FCTD_responsematched.mat','figpath','./figures/','plotevery',10)
%
%
%
% Alex Andriatis
% 2021-02-01

% Define structure and transfer auxilary data

P=inputParser;
addRequired(P,'FCTDprofiles',@isstruct);

addRequired(P,'FCTDmatched',@isstruct);

validationFcn = @(x) validateattributes(x,{'char','string'},{});
addRequired(P,'datapath',validationFcn);

defaultFigpath = '';
validationFcn = @(x) validateattributes(x,{'char','string'},{});
addParameter(P,'figpath',defaultFigpath,validationFcn);

defaultPlotevery = 0;
addParameter(P,'plotevery',defaultPlotevery,@isnumeric);

defaultParameters = [];
addParameter(P,'Parameters',defaultParameters,@isstruct);

parse(P,FCTDprofiles,FCTDmatched,datapath,varargin{:});
FCTDprofiles = P.Results.FCTDprofiles;
FCTDmatched = P.Results.FCTDmatched;
datapath = P.Results.datapath;
figpath = P.Results.figpath;
plotevery = P.Results.plotevery;
Parameters = P.Results.Parameters;

% Turn on plotting
if ~isempty(figpath)
    tosave=1;
    fpath = fullfile(figpath,'Corrected');
    mkdir(fpath);
end
serials = unique(FCTDprofiles.header.serial);

%% Combine data into one structure

FCTDcorrected=FCTDprofiles;
directions = {'up','down'};
variables = {'temperature','conductivity','pressure'};
for i=1:length(directions)
    direction = directions{i};
    for j=1:length(variables)
        variable = variables{j};
        FCTDcorrected.(direction).(variable) = FCTDmatched.dehyst.(direction).(variable);
    end
end

%% Cut off the ends of the profiles
% The hysteresis and response matching codes create sections of false data

nprof = length(FCTDprofiles.down.time);
variables = fieldnames(FCTDcorrected.down);
for k=1:nprof
    for j=1:2
        direction = directions{j};
        Iserial = find(serials==FCTDcorrected.header.serial(k));
        
        start_cut=[];
        end_cut=[];
        if ~isempty(Parameters)
            lag = Parameters.TParams(Iserial,1);
            offset_T = Parameters.offset_T;
            offset_C = Parameters.offset_C;
            start_cut = 2*lag+offset_T(offset_T>=0)+offset_C(offset_C>=0);
            end_cut = offset_T(offset_T<=0)+offset_C(offset_C<=0);
        end
        transient=[];
        if isfield(FCTDmatched.dehyst.(direction),'transient')
            transient = FCTDmatched.dehyst.(direction).transient{k};
        end
        if isempty(start_cut)
            start_cut=0;
        end
        if isempty(end_cut)
            end_cut = 0;
        end
        if isempty(transient)
            transient=0;
        end
        start_cut = max(start_cut,transient);
        end_cut = max(end_cut,transient);
        
        for n=1:length(variables)
            variable = variables{n};
            FCTDcorrected.(direction).(variable){k} = FCTDcorrected.(direction).(variable){k}(start_cut+1:end-end_cut,:);
        end     
    end
end

%% Save the data

save(datapath,'-struct','FCTDcorrected','-v7.3');
disp(['Saved Corrected FCTD profiles in ' datapath]);


%% Plot profiles
% Plot tvp cvp and tvc
if plotevery
for k=1:plotevery:nprof
    
figure
    clf;
        ax(1)=subplot(2,3,1);
        plot(FCTDprofiles.up.temperature{k},FCTDprofiles.up.pressure{k});
        hold on;
        plot(FCTDprofiles.down.temperature{k},FCTDprofiles.down.pressure{k});
        xlabel('Temperature');
        ylabel('Pressure');
        title(['Raw T #' num2str(k)]);
        legend('Up','Down','Location','best');
        grid on;
        set(gca,'YDir','reverse');
        
    ax(2)=subplot(2,3,2);
        plot(FCTDprofiles.up.conductivity{k},FCTDprofiles.up.pressure{k});
        hold on;
        plot(FCTDprofiles.down.conductivity{k},FCTDprofiles.down.pressure{k});
        xlabel('Conductivity');
        ylabel('Pressure');
        title(['Raw C #' num2str(k)]);
        legend('Up','Down','Location','best');
        grid on;
        set(gca,'YDir','reverse');
        
    ax(3)=subplot(2,3,3);
        plot(FCTDprofiles.down.conductivity{k},FCTDprofiles.down.temperature{k});
        hold on;
        plot(FCTDprofiles.up.conductivity{k},FCTDprofiles.up.temperature{k});
        xlabel('Conductivity');
        ylabel('Temperautre');
        title(['Raw T-C #' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
        
    ax(4)=subplot(2,3,4);
        plot(FCTDcorrected.up.temperature{k},FCTDcorrected.up.pressure{k});
        hold on;
        plot(FCTDcorrected.down.temperature{k},FCTDcorrected.down.pressure{k});
        xlabel('Temperature');
        ylabel('Pressure');
        title(['Corrected T #' num2str(k)]);
        legend('Up','Down','Location','best');
        grid on;
        set(gca,'YDir','reverse');
        
    ax(5)=subplot(2,3,5);
        plot(FCTDcorrected.up.conductivity{k},FCTDcorrected.up.pressure{k});
        hold on;
        plot(FCTDcorrected.down.conductivity{k},FCTDcorrected.down.pressure{k});
        xlabel('Conductivity');
        ylabel('Pressure');
        title(['Corrected C #' num2str(k)]);
        legend('Up','Down','Location','best');
        grid on;
        set(gca,'YDir','reverse');
        
    ax(6)=subplot(2,3,6);
        plot(FCTDcorrected.down.conductivity{k},FCTDcorrected.down.temperature{k});
        hold on;
        plot(FCTDcorrected.up.conductivity{k},FCTDcorrected.up.temperature{k});
        xlabel('Conductivity');
        ylabel('Temperautre');
        title(['Corrected T-C #' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
 
        linkaxes(ax([1:2 4:5]),'y');
        linkaxes(ax([1 4]),'x');
        linkaxes(ax([2 5]),'x');
        
    if tosave
        fname = ['corrected_profile_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
    pause(1);
    
end
end
end

        
        



