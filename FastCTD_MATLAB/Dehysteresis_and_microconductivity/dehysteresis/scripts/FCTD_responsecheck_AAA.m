function FCTD_responsecheck_AAA(FCTDprofiles,FCTDmatched,varargin)
%
% FCTD_responsecheck_AAA is a re-write of Rob's
% L_2016B_ResponseCorrectionVerificationFLEAT.m
%
% It now works with a profile structure saved as cell arrays.
%
% Inputs:
%   FCTDprofiles: A strucure of profiles which is the output of
%   FCTD_FindProfiles_AAA.
%
%   FCTDmatched: A structure of hysteresis-corrected and response-matched profiles which is
%   the output of FCTD_responsematch_AAA.
%
% Alex Andriatis
% 2021-02-01

P=inputParser;
addRequired(P,'FCTDprofiles',@isstruct);

addRequired(P,'FCTDmatched',@isstruct);

defaultFigpath = '';
validationFcn = @(x) validateattributes(x,{'char','string'},{});
addParameter(P,'figpath',defaultFigpath,validationFcn);

defaultPlotevery = 0;
addParameter(P,'plotevery',defaultPlotevery,@isnumeric);

parse(P,FCTDprofiles,FCTDmatched,varargin{:});
FCTDprofiles = P.Results.FCTDprofiles;
FCTDmatched = P.Results.FCTDmatched;
figpath = P.Results.figpath;

%% The rest of this code is just used for plotting some diagnostics

%% Check input structures
reqnames = {'temperature','conductivity','pressure'};
directions = {'up','down'};

for i=1:2
    if ~all(isfield(FCTDprofiles,directions{i}))
        error(['The FCTD profile file does not have the ' directions{i} ' direction']);
    end
    if ~all(isfield(FCTDprofiles.(directions{i}),reqnames))
        error(['One or more required fields in FCTD is missing from the ' directions{i} ' direction']);
    end
end

if ~all(isfield(FCTDmatched,'dehyst'))
    error(['The FCTD matched profile file does not have the dehyst field']);
end
if ~all(isfield(FCTDmatched,'header'))
    error(['The FCTD matched profile file does not have a header field']);
end
for j=1:2
    if ~all(isfield(FCTDmatched.dehyst,directions{i}))
        error(['The FCTD matched profile file does not have the ' directions{i} ' direction in the dehyst field']);
    end
    if ~all(isfield(FCTDmatched.dehyst.(directions{i}),reqnames))
        error(['One or more required fields in FCTDmatched is missing from the ' directions{i} ' direction in the dehyst field']);
    end
end

nprof = length(FCTDprofiles.down.time);
fNy = round(1/(2*mean(diff(FCTDprofiles.down.time{1})*86400))); % Nyquist frequency of sample rate in cycles per second. Nomially 8 for 16 Hz sampling.
dt = 1/(2*fNy);
serials = unique(FCTDmatched.header.serial);

% Turn on plotting
if ~isempty(figpath)
    tosave=1;
    fpath = fullfile(figpath,'ResponseCheck_Deep');
    mkdir(fpath);
end


%% What to call the two kinds of profiles

profiletypes={'raw','corrected'};
%% Compute Average spectra for the raw and corrected profiles


Spectra=[];
% Compute Spectra
% Loop over raw and corrected
for i=1:2
%for i=1:1 % Testing
    switch i
        case 1
            FCTD = FCTDprofiles;
        case 2
            FCTD = FCTDmatched.dehyst;
    end
    profiletype=profiletypes{i};
    
    % Loop over up and down
    for j=1:2
%    for j=1:1 % Testing
        direction = directions{j};
        data = FCTD.(direction);
        % Loop over each profile
        for k=1:nprof
%        for k=1:1 % Testing

            T = data.temperature{k};
            C = data.conductivity{k}; 
            
            % Calculate spectra and cross-spectrum of dT/dt and dC/dt
            dTdt = diff(T);
            dCdt = diff(C);
            
            [f,t]=ss_fft_AAA(dTdt,dt);
            [~,c]=ss_fft_AAA(dCdt,dt);
            
            CRspec = c.*conj(t);
            Tspec = abs(t).^2;
            Cspec = abs(c).^2;
            
            Spectra.(profiletype).(direction).f{k}=f;
            Spectra.(profiletype).(direction).Tspec{k}=Tspec;
            Spectra.(profiletype).(direction).Cspec{k}=Cspec;
            Spectra.(profiletype).(direction).CRspec{k}=CRspec;

        end
    end
end
           
% Average Spectra
Spectra_mean=[];
% Loop over serial number

for l=1:length(serials)
    disp(['Calculating average spectra for CTD serial number ' num2str(serials(l))]);
    Iserial = find(FCTDmatched.header.serial == serials(l) & FCTDmatched.header.qc==0);
    % Loop over raw and corrected
    for i=1:2
        profiletype=profiletypes{i};
        % Loop over up and down
        for j=1:2
    %    for j=1:1 % Testing
    %        j=1; % Testing
            direction = directions{j};

            Specs = Spectra.(profiletype).(direction);

            % Need to figure out what the longest spectrum record is, use that
            % frequency base for averaging

            Max=0;
            fgrid=[];

            for k=1:length(Iserial)
            %for k=1:1 % Testing
            %k=1; % Testing
                speclen = length(Specs.f{Iserial(k)});
                if speclen>Max
                    Max=speclen;
                    fgrid=Specs.f{Iserial(k)};
                end
            end

            Tspecavg = NaN(Max,nprof);
            Cspecavg = NaN(Max,nprof);
            CRspecavg = NaN(Max,nprof);

            % Interpolate spectra onto the same grid, then nanmean to get the
            % average spectra

            for k=1:length(Iserial)
    %        for k=1:1 % Testing
    %        k=1; % Testing

                f=Specs.f{Iserial(k)};
                Tspec=Specs.Tspec{Iserial(k)};
                Cspec=Specs.Cspec{Iserial(k)};
                CRspec=Specs.CRspec{Iserial(k)};

                Tspec = interp1(f,Tspec,fgrid);
                Cspec = interp1(f,Cspec,fgrid);
                CRspec = interp1(f,CRspec,fgrid);

                Tspecavg(:,k)=Tspec;
                Cspecavg(:,k)=Cspec;
                CRspecavg(:,k)=CRspec;
            end

            Tspecavg = mean(Tspecavg,2,'omitnan');
            Cspecavg = mean(Cspecavg,2,'omitnan');
            CRspecavg = mean(CRspecavg,2,'omitnan');

            % Coherence
            Coherence = abs(CRspecavg).^2./Tspecavg./Cspecavg;

            % Phase, degrees
            Phase = angle(CRspecavg)*180/pi;

            % Transfer Function
            Transfer = abs(CRspecavg)./Tspecavg;

            Spectra_mean.(profiletype).(direction).f{l} = fgrid;
            Spectra_mean.(profiletype).(direction).Tspec{l} = Tspecavg;
            Spectra_mean.(profiletype).(direction).Cspec{l} = Cspecavg;
            Spectra_mean.(profiletype).(direction).CRspec{l} = CRspecavg;
            Spectra_mean.(profiletype).(direction).Coherence{l} = Coherence;
            Spectra_mean.(profiletype).(direction).Phase{l} = Phase;
            Spectra_mean.(profiletype).(direction).Transfer{l} = Transfer;
        end
    end
end


%% Rob's plots checking of the coherence, phase, and transfer function worked


for l=1:length(serials)
    disp(['Displaying average spectra for CTD serial number ' num2str(serials(l))]);
    figure
        clf;
        subplot(2,2,1);
            data=Spectra_mean.raw.up;
            plot(data.f{l},data.Coherence{l});
            hold on;
            data=Spectra_mean.corrected.up;
            plot(data.f{l},data.Coherence{l});
            title('T C Up Coherence Squared');
            xlabel('Frequency [Hz]');
            legend('Raw','Corrected','Location','best');
            grid on;

        subplot(2,2,2);
            data=Spectra_mean.raw.down;
            plot(data.f{l},data.Coherence{l});
            hold on;
            data=Spectra_mean.corrected.down;
            plot(data.f{l},data.Coherence{l});
            title('T C Down Coherence Squared');
            xlabel('Frequency [Hz]');
            legend('Raw','Corrected','Location','best');
            grid on;

        subplot(2,2,3);
            data=Spectra_mean.raw.up;
            plot(data.f{l},data.Phase{l});
            hold on;
            data=Spectra_mean.corrected.up;
            plot(data.f{l},data.Phase{l});
            title('T C Up Phase (degrees)');
            xlabel('Frequency [Hz]');
            legend('Raw','Corrected','Location','best');
            grid on;

        subplot(2,2,4);
            data=Spectra_mean.raw.down;
            plot(data.f{l},data.Phase{l});
            hold on;
            data=Spectra_mean.corrected.down;
            plot(data.f{l},data.Phase{l});
            title('T C Down Phase (degrees)');
            xlabel('Frequency [Hz]');
            legend('Raw','Corrected','Location','best');
            grid on;
            
            suplabel(['Serial # ' num2str(serials(l))],'t');
            
    if tosave
        fname = ['check_TC_Coherence_Phase_serial_' num2str(serials(l))];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
    pause(0.5)
    

    figure
        clf;
        subplot(2,2,1);
            data=Spectra_mean.raw.up;
            plot(data.f{l},data.Transfer{l});
            hold on;
            data=Spectra_mean.corrected.up;
            plot(data.f{l},data.Transfer{l});
            title('T C Up Transfer Function');
            xlabel('Frequency [Hz]');
            legend('Raw','Corrected','Location','best');
            grid on;

        subplot(2,2,2);
            data=Spectra_mean.raw.down;
            plot(data.f{l},data.Transfer{l});
            hold on;
            data=Spectra_mean.corrected.down;
            plot(data.f{l},data.Transfer{l});
            title('T C Down Transfer Function');
            xlabel('Frequency [Hz]');
            legend('Raw','Corrected','Location','best');
            grid on;  

        subplot(2,2,3);
            data=Spectra_mean.raw.up;
            plot(data.f{l},data.Tspec{l});
            hold on;
            data=Spectra_mean.raw.down;
            plot(data.f{l},data.Tspec{l});
            data=Spectra_mean.corrected.up;
            plot(data.f{l},data.Tspec{l});
            data=Spectra_mean.corrected.down;
            plot(data.f{l},data.Tspec{l});
            title('Log of dT Spectra');
            legend('Raw Up','Raw Down','Corrected Up','Corrected Down');
            set(gca, 'YScale', 'log')

        subplot(2,2,4);
            data=Spectra_mean.raw.up;
            plot(data.f{l},data.Cspec{l});
            hold on;
            data=Spectra_mean.raw.down;
            plot(data.f{l},data.Cspec{l});
            data=Spectra_mean.corrected.up;
            plot(data.f{l},data.Cspec{l});
            data=Spectra_mean.corrected.down;
            plot(data.f{l},data.Cspec{l});
            title('Log of dC Spectra');
            legend('Raw Up','Raw Down','Corrected Up','Corrected Down');
            set(gca, 'YScale', 'log')
            
            suplabel(['Serial # ' num2str(serials(l))],'t');
            
    if tosave
        fname = ['check_TC_Transfer_Spectra_serial_' num2str(serials(l))];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
    pause(0.5)
end