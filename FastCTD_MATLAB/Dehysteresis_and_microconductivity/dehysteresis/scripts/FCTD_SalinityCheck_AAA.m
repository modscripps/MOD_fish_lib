function FCTD_SalinityCheck_AAA(FCTDprofiles,FCTDdehyst,FCTDmatched,varargin)
% Plots salinity profiles on a per-profile basis to check the effects of 
% the hysteresis and response matching codes
%
%
% Alex Andriatis
% 2021-02-01

P=inputParser;
addRequired(P,'FCTDprofiles',@isstruct);

addRequired(P,'FCTDdehyst',@isstruct);

addRequired(P,'FCTDmatched',@isstruct);

defaultFigpath = '';
validationFcn = @(x) validateattributes(x,{'char','string'},{});
addParameter(P,'figpath',defaultFigpath,validationFcn);

defaultPlotevery = 0;
addParameter(P,'plotevery',defaultPlotevery,@isnumeric);

defaultProfile = 0;
addParameter(P,'profile',defaultProfile,@isnumeric);

parse(P,FCTDprofiles,FCTDdehyst,FCTDmatched,varargin{:});
FCTDprofiles = P.Results.FCTDprofiles;
FCTDdehyst = P.Results.FCTDdehyst;
FCTDmatched = P.Results.FCTDmatched;
figpath = P.Results.figpath;
plotevery = P.Results.plotevery;
profile = P.Results.profile;

% Turn on saving
if ~isempty(figpath)
    tosave=1;
    fpath = fullfile(figpath,'SalinityCheck');
    mkdir(fpath);
end

nprof = length(FCTDprofiles.down.time);
% Turn on multiple profiles or just one
if profile>0
    toplot = profile;
elseif plotevery>0
    toplot = 1:plotevery:nprof;
else
    toplot=1;
end

%% Plot profiles

for k=toplot
    
% Plot the t-c profile
figure
clf;
    ax(1)=subplot(1,4,1)
        plot(FCTDprofiles.down.conductivity{k},FCTDprofiles.down.temperature{k});
        hold on;
        plot(FCTDprofiles.up.conductivity{k},FCTDprofiles.up.temperature{k});
        xlabel('Conductivity');
        ylabel('Temperautre');
        title(['T-C Raw: #' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
       
    ax(2)=subplot(1,4,2)
        plot(FCTDdehyst.down.conductivity{k},FCTDdehyst.down.temperature{k});
        hold on;
        plot(FCTDdehyst.up.conductivity{k},FCTDdehyst.up.temperature{k});
        xlabel('Conductivity');
        ylabel('Temperautre');
        title(['T-C deHyst: #' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;

    ax(3)=subplot(1,4,3)
        plot(FCTDmatched.raw.down.conductivity{k},FCTDmatched.raw.down.temperature{k});
        hold on;
        plot(FCTDmatched.raw.up.conductivity{k},FCTDmatched.raw.up.temperature{k});
        xlabel('Conductivity');
        ylabel('Temperautre');
        title(['T-C deSpiked: #' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
        
    ax(4)=subplot(1,4,4)
        plot(FCTDmatched.dehyst.down.conductivity{k},FCTDmatched.dehyst.down.temperature{k});
        hold on;
        plot(FCTDmatched.dehyst.up.conductivity{k},FCTDmatched.dehyst.up.temperature{k});
        xlabel('Conductivity');
        ylabel('Temperautre');
        title(['T-C deHyst and deSpiked: #' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
 
        linkaxes(ax);
        
    if tosave
        fname = ['salinitycheck_TC_profile_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
    pause(0.5);

% Plot the s-p profile
figure
clf;
    ax(1)=subplot(1,4,1)
        SP = gsw_SP_from_C(FCTDprofiles.down.conductivity{k}*10,FCTDprofiles.down.temperature{k},FCTDprofiles.down.pressure{k});
        plot(SP,FCTDprofiles.down.pressure{k});
        hold on;
        SP = gsw_SP_from_C(FCTDprofiles.up.conductivity{k}*10,FCTDprofiles.up.temperature{k},FCTDprofiles.up.pressure{k});
        plot(SP,FCTDprofiles.up.pressure{k});
        xlabel('Salinity');
        ylabel('Pressure');
        title(['S-P Raw: #' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
        
    ax(2)=subplot(1,4,2)
        SP = gsw_SP_from_C(FCTDdehyst.down.conductivity{k}*10,FCTDdehyst.down.temperature{k},FCTDdehyst.down.pressure{k});
        plot(SP,FCTDdehyst.down.pressure{k});
        hold on;
        SP = gsw_SP_from_C(FCTDdehyst.up.conductivity{k}*10,FCTDdehyst.up.temperature{k},FCTDdehyst.up.pressure{k});
        plot(SP,FCTDdehyst.up.pressure{k});
        xlabel('Salinity');
        ylabel('Pressure');
        title(['S-P deHyst: #' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;

    ax(3)=subplot(1,4,3)
        SP = gsw_SP_from_C(FCTDmatched.raw.down.conductivity{k}*10,FCTDmatched.raw.down.temperature{k},FCTDmatched.raw.down.pressure{k});
        plot(SP,FCTDmatched.raw.down.pressure{k});
        hold on;
        SP = gsw_SP_from_C(FCTDmatched.raw.up.conductivity{k}*10,FCTDmatched.raw.up.temperature{k},FCTDmatched.raw.up.pressure{k});
        plot(SP,FCTDmatched.raw.up.pressure{k});
        xlabel('Salinity');
        ylabel('Pressure');
        title(['S-P deSpiked: #' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
        
    ax(4)=subplot(1,4,4)
        SP = gsw_SP_from_C(FCTDmatched.dehyst.down.conductivity{k}*10,FCTDmatched.dehyst.down.temperature{k},FCTDmatched.dehyst.down.pressure{k});
        plot(SP,FCTDmatched.dehyst.down.pressure{k});
        hold on;
        SP = gsw_SP_from_C(FCTDmatched.dehyst.up.conductivity{k}*10,FCTDmatched.dehyst.up.temperature{k},FCTDmatched.dehyst.up.pressure{k});
        plot(SP,FCTDmatched.dehyst.up.pressure{k});
        xlabel('Salinity');
        ylabel('Pressure');
        title(['S-P deHyst and deSpiked: #' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
 
        linkaxes(ax);
        set(ax,'YDir','Reverse');
        
    if tosave
        fname = ['salinitycheck_SP_profile_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
    pause(0.5);
    
end
end