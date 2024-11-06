function FCTDmatched=FCTD_responsematch_AAA(FCTDprofiles,FCTDdehyst,datapath,varargin)
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
%   FCTDmatched=FCTD_responsematch_AAA(FCTDprofiles_qc,FCTDdehyst,'./Matching/FCTD_responsematched.mat')
%   FCTDmatched=FCTD_responsematch_AAA(FCTDprofiles_qc,FCTDdehyst,'./Matching/FCTD_responsematched.mat','figpath','./figures/','plotevery',10)
%
%
%
% Alex Andriatis
% 2021-02-01

% Define structure and transfer auxilary data
FCTDmatched=[];

P=inputParser;
addRequired(P,'FCTDprofiles',@isstruct);

addRequired(P,'FCTDdehyst',@isstruct);

validationFcn = @(x) validateattributes(x,{'char','string'},{});
addRequired(P,'datapath',validationFcn);

defaultFigpath = '';
validationFcn = @(x) validateattributes(x,{'char','string'},{});
addParameter(P,'figpath',defaultFigpath,validationFcn);

defaultPlotevery = 0;
addParameter(P,'plotevery',defaultPlotevery,@isnumeric);

defaultDepthRange = 0;
addParameter(P,'depthrange',defaultDepthRange,@isnumeric);

parse(P,FCTDprofiles,FCTDdehyst,datapath,varargin{:});
FCTDprofiles = P.Results.FCTDprofiles;
FCTDdehyst = P.Results.FCTDdehyst;
datapath = P.Results.datapath;
figpath = P.Results.figpath;
plotevery = P.Results.plotevery;
depthrange = P.Results.depthrange;


%% Checking the input files:

reqnames = {'temperature','conductivity','pressure'};
directions = {'up','down'};
profiletypes={'raw','dehyst'};

if ~isfield(FCTDprofiles,'header')
    error(['The FCTD profile file does not have a header']);
end

if ~isfield(FCTDdehyst,'header')
    error(['The FCTDdehyst file does not have a header']);
end


for i=1:2
    if ~all(isfield(FCTDprofiles,directions{i}))
        error(['The FCTD profile file does not have the ' directions{i} ' direction']);
    end
    if ~all(isfield(FCTDprofiles.(directions{i}),reqnames))
        error(['One or more required fields in FCTD is missing from the ' directions{i} ' direction']);
    end
end

for i=1:2
    if ~all(isfield(FCTDdehyst,directions{i}))
        error(['The FCTD hysteresis-corrected profile file does not have the ' directions{i} ' direction']);
    end
    if ~all(isfield(FCTDdehyst.(directions{i}),reqnames))
        error(['One or more required fields in FCTDdehyst is missing from the ' directions{i} ' direction']);
    end
end

nprof = length(FCTDprofiles.down.time);
fNy = round(1/(2*mean(diff(FCTDprofiles.down.time{1})*86400))); % Nyquist frequency of sample rate in cycles per second. Nomially 8 for 16 Hz sampling.
dt = 1/(2*fNy);

% Turn on plotting
if ~isempty(figpath)
    tosave=1;
    fpath = fullfile(figpath,'ResponseMatching_Deep');
    mkdir(fpath);
end

%% Not changing the header so carry it forward

FCTDmatched.header = FCTDprofiles.header;

%% Separate the data by serial number of FASTCTD for the optimization
serials = unique(FCTDdehyst.header.serial);

%% 1. For each profile, get the specrta of p, t, c, and the cross-spectrum

Spectra=[];

% Loop over raw and corrected
for i=1:2
    switch i
        case 1
            FCTD = FCTDprofiles;
        case 2
            FCTD = FCTDdehyst;
    end
    profiletype=profiletypes{i};
    
    % Loop over up and down
    for j=1:2
        direction = directions{j};
        
        % Loop over each profile
        for k=1:nprof

            T = FCTD.(direction).temperature{k};
            C = FCTD.(direction).conductivity{k};
            P = FCTD.(direction).pressure{k};
            
            if length(depthrange)==1
                I=find(P>=depthrange);
            else
                I=find(P>=depthrange(1) & P<=depthrange(end));
            end
            if length(I)<16
                continue;
            end
             
            % Calculate spectra and cross-spectrum of dT/dt and dC/dt
            dTdt = diff(T(I));
            dCdt = diff(C(I));
            
            % Maybe it's worth removing the trend too to reduce spectral
            % leakage? But near the top and bottom of the profile, the
            % instrument is moving slowly, so dT/dt is close to zero
            % already. Removing any sort of mean or trend would actually
            % shift the end points away from zero, which is something I'd
            % like to avoid.
            
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
           


%% 2. Average the spectra according to serial number
Spectra_mean=[];
for l=1:length(serials)
    disp(['Calculating average spectra for CTD serial number ' num2str(serials(l))]);
    Iserial = find(FCTDdehyst.header.serial == serials(l) & FCTD.header.qc==0); % Only use good data
    
    % Loop over raw and corrected
    for i=1:2
        profiletype=profiletypes{i};
        % Loop over up and down
        for j=1:2
            direction = directions{j};

            Specs = Spectra.(profiletype).(direction);

            % Need to figure out what the longest spectrum record is, use that
            % frequency base for averaging

            Max=0;
            fgrid=[];

            for k=1:length(Iserial)
                I = Iserial(k);
                speclen = length(Specs.f{I});
                if speclen>Max
                    Max=speclen;
                    fgrid=Specs.f{I};
                end
            end
            speclen = length(fgrid);

            Tspecavg = NaN(speclen,nprof);
            Cspecavg = NaN(speclen,nprof);
            CRspecavg = NaN(speclen,nprof);

            % Interpolate spectra onto the same grid, then nanmean to get the
            % average spectra

            for k=1:length(Iserial)
                I = Iserial(k);
                
                f=Specs.f{I};
                Tspec=Specs.Tspec{I};
                Cspec=Specs.Cspec{I};
                CRspec=Specs.CRspec{I};
                if length(f)<16
                    continue
                end

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

            Spectra_mean.(profiletype).(direction).f{l} = fgrid;
            Spectra_mean.(profiletype).(direction).Tspec{l} = Tspecavg;
            Spectra_mean.(profiletype).(direction).Cspec{l} = Cspecavg;
            Spectra_mean.(profiletype).(direction).CRspec{l} = CRspecavg;
        end
    end
end


%% Rob's plots of the average spectra
if plotevery %switch to turn on plotting
for l=1:length(serials)
    figure
    clf;
        subplot(2,2,1);
            data=Spectra_mean.raw.up;
            plot(data.f{l},abs(data.CRspec{l}.^2)./data.Tspec{l}./data.Cspec{l},'b');
            hold on;
            data=Spectra_mean.dehyst.up;    
            plot(data.f{l},abs(data.CRspec{l}.^2)./data.Tspec{l}./data.Cspec{l},'r');
            title(['T C Up Coherence Squared: Serial ' num2str(serials(l))]);
            xlabel('frequency Hz');
            legend('Raw','deHyst','Location','best');
            grid on;

        subplot(2,2,2);
            data=Spectra_mean.raw.down;
            plot(data.f{l},abs(data.CRspec{l}.^2)./data.Tspec{l}./data.Cspec{l},'b');
            hold on;
            data=Spectra_mean.dehyst.down;    
            plot(data.f{l},abs(data.CRspec{l}.^2)./data.Tspec{l}./data.Cspec{l},'r');
            title(['T C Down Coherence Squared: Serial ' num2str(serials(l))]);
            xlabel('frequency Hz');
            legend('Raw','deHyst','Location','best');
            grid on;

        subplot(2,2,3);
            data=Spectra_mean.raw.up;
            plot(data.f{l},angle(data.CRspec{l})*180/pi,'b');
            hold on;
            data=Spectra_mean.dehyst.up;
            plot(data.f{l},angle(data.CRspec{l})*180/pi,'r');
            title(['T C phase, degrees, Up: Serial ' num2str(serials(l))]);
            xlabel('frequency Hz');
            legend('Raw','deHyst','Location','best'); 
            grid on;

         subplot(2,2,4);
            data=Spectra_mean.raw.down;
            plot(data.f{l},angle(data.CRspec{l})*180/pi,'b');
            hold on;
            data=Spectra_mean.dehyst.down;
            plot(data.f{l},angle(data.CRspec{l})*180/pi,'r');
            title(['T C phase, degrees, Down: Serial ' num2str(serials(l))]);
            xlabel('frequency Hz');
            legend('Raw','deHyst','Location','best');
            grid on;

    if tosave
        fname = ['matching_TC_coherence_phase_' num2str(serials(l))];
        saveas(gcf,fullfile(fpath,fname),'png');
    end
    pause(0.5);

    figure
    clf;
        subplot(2,2,1);
            data=Spectra_mean.raw.up;
            plot(data.f{l},abs(data.CRspec{l})./data.Tspec{l},'b');
            hold on;
            data=Spectra_mean.raw.down;
            plot(data.f{l},abs(data.CRspec{l})./data.Tspec{l},'r');
            title(['T C Transfer Function: Serial ' num2str(serials(l))]);
            xlabel('frequency Hz');
            legend('Up','Down','Location','best'); 
            grid on;

        subplot(2,2,2);
            data=Spectra_mean.dehyst.up;
            plot(data.f{l},abs(data.CRspec{l})./data.Tspec{l},'b');
            hold on;
            data=Spectra_mean.dehyst.down;
            plot(data.f{l},abs(data.CRspec{l})./data.Tspec{l},'r');
            title(['T C Transfer Function, Hysteresis Corrected: Serial ' num2str(serials(l))]);
            xlabel('frequency Hz');
            legend('Up','Down','Location','best'); 
            grid on;

        subplot(2,2,3);
            data=Spectra_mean.raw.up;
            loglog(data.f{l},data.Tspec{l},'b');
            hold on;
            data=Spectra_mean.raw.down;
            loglog(data.f{l},data.Tspec{l},'r');        
            data=Spectra_mean.dehyst.up;
            loglog(data.f{l},data.Tspec{l},'m');        
            data=Spectra_mean.dehyst.down;
            loglog(data.f{l},data.Tspec{l},'c');         
            title(['dT/dt Spectra: Serial ' num2str(serials(l))]);
            xlabel('frequency Hz');
            legend('Raw Up','Raw Down','deHyst Up','deHyst Down','Location','best');
            grid on;

        subplot(2,2,4);
            data=Spectra_mean.raw.up;
            loglog(data.f{l},data.Cspec{l},'b');
            hold on;
            data=Spectra_mean.raw.down;
            loglog(data.f{l},data.Cspec{l},'r');        
            data=Spectra_mean.dehyst.up;
            loglog(data.f{l},data.Cspec{l},'m');        
            data=Spectra_mean.dehyst.down;
            loglog(data.f{l},data.Cspec{l},'c');         
            title(['dC/dt Spectra: Serial ' num2str(serials(l))]);
            xlabel('frequency Hz');
            legend('Raw Up','Raw Down','deHyst Up','deHyst Down','Location','best');
            grid on;
            
    if tosave
        fname = ['matching_TC_transfer_spectra_' num2str(serials(l))];
        saveas(gcf,fullfile(fpath,fname),'png');
    end
    pause(0.5);
end
end

%% Calculate polynomial coefficients of the Gain and Phase

% Loop over serial numbers
for l=1:length(serials)
    % Loop over raw and corrected
    for i=1:2
        profiletype=profiletypes{i};   
        % Loop over up and down
        for j=1:2
            direction = directions{j};

            data = Spectra_mean.(profiletype).(direction);
            f = data.f{l};
            
            Gain = abs(data.CRspec{l}).^2./data.Tspec{l}.^2;
            Gain = sqrt(Gain);
            Gain = reshape(Gain,[1,length(Gain)]);
            
            Phase = angle(data.CRspec{l});
            Phase = reshape(Phase,[1,length(Phase)]);            
            for p=2:length(Phase)
                while Phase(p)-Phase(p-1)>pi/32
                    Phase(p)=Phase(p)-pi/32;
                end
                while Phase(p)-Phase(p-1)<-pi/32
                    Phase(p)=Phase(p)+pi/32;
                end
            end
            
            % Make Gain an even function of frequency
            GainArray = [flip(Gain(2:end)) Gain];
            % Make Phase an odd function of frequency
            PhaseArray = [flip(-Phase(2:end)) Phase];
            
            fArray = [flip(-f(2:end)) f];

            % Fit a polynomial to the Gain
            [Again,~]=polyfit(fArray,GainArray,20);
            GainFit=polyval(Again,f);

            % Fit a polynomial to the Phase
            [Aphase,~]=polyfit(fArray,PhaseArray,20);
            PhaseFit=polyval(Aphase,f);

            Spectra_mean.(profiletype).(direction).Again{l} = Again;
            Spectra_mean.(profiletype).(direction).GainFit{l} = GainFit;
            Spectra_mean.(profiletype).(direction).Aphase{l} = Aphase;
            Spectra_mean.(profiletype).(direction).PhaseFit{l} = PhaseFit;

            % Some simple plots to check that the fits are working along the
            % way
            if plotevery
            figure
                clf;
                subplot(2,1,1)
                    plot(f,Gain);
                    hold on;
                    plot(f,GainFit);
                    ylabel('Gain');
                    legend('Gain','PolyFit','Location','best');
                    title(['Fits to Gain and Phase for ' profiletype ' profiles in the ' direction ' direction for CTD serial number ' num2str(serials(l))]);

                subplot(2,1,2)
                    plot(f,Phase);
                    hold on;
                    plot(f,PhaseFit);
                    ylabel('Phase');
                    legend('Phase','PolyFit','Location','best');
                    
              if tosave
                fname = ['matching_GainFit_PhaseFit_' num2str(serials(l)) '_' profiletype '_' direction];
                saveas(gcf,fullfile(fpath,fname),'png');
              end
              pause(0.5);
            end
        end
    end
end


%% Rob's plots of gain and phase fits
if plotevery
% Loop over serial numbers
for l=1:length(serials)
    figure
    clf;
        subplot(2,2,1)
            data=Spectra_mean.raw.up;
            plot(data.f{l},data.GainFit{l});
            hold on;
            data=Spectra_mean.raw.down;
            plot(data.f{l},data.GainFit{l});

            ylabel('Gain');
            xlabel('Frequency (Hz)');
            title(['Raw Gain: Serial ' num2str(serials(l))]);
            legend('Up','Down','Location','best');

        subplot(2,2,2)
            data=Spectra_mean.dehyst.up;
            plot(data.f{l},data.GainFit{l});
            hold on;
            data=Spectra_mean.dehyst.down;
            plot(data.f{l},data.GainFit{l});

            ylabel('Gain');
            xlabel('Frequency (Hz)');
            title(['deHyst Gain: Serial ' num2str(serials(l))]);
            legend('Up','Down','Location','best');       


        subplot(2,2,3)
            data=Spectra_mean.raw.up;
            plot(data.f{l},data.PhaseFit{l}*180/pi);
            hold on;
            data=Spectra_mean.raw.down;
            plot(data.f{l},data.PhaseFit{l}*180/pi);

            ylabel('Phase');
            xlabel('Frequency (Hz)');
            title(['Raw Phase: Serial ' num2str(serials(l))]);
            legend('Up','Down','Location','best');

        subplot(2,2,4)
            data=Spectra_mean.dehyst.up;
            plot(data.f{l},data.PhaseFit{l}*180/pi);
            hold on;
            data=Spectra_mean.dehyst.down;
            plot(data.f{l},data.PhaseFit{l}*180/pi);

            ylabel('Phase');
            xlabel('Frequency (Hz)');
            title(['deHyst Phase: Serial ' num2str(serials(l))]);
            legend('Up','Down','Location','best');
        if tosave
            fname = ['matching_GainFit_PhaseFit_' num2str(serials(l))];
            saveas(gcf,fullfile(fpath,fname),'png');
        end
        pause(0.5);
end
end
        
        
 %% Correcting the profiles
 
 % Loop over raw and corrected
for i=1:2

    switch i
        case 1
            FCTD = FCTDprofiles;
        case 2
            FCTD = FCTDdehyst;
    end
    profiletype=profiletypes{i};
    
    % Loop over up and down
    for j=1:2
  
        direction = directions{j};
        
        data = FCTD.(direction);
        
        % Loop over each profile
        for k=1:nprof

            Iserial = find(serials==FCTDdehyst.header.serial(k));

            % 1. Calculate FFT of Temperature
            
            T = data.temperature{k};
            C = data.conductivity{k};
            P = data.pressure{k};
            
            % Two different modes of calculating the response matching -
            % using dT and T
            
            % I've found that using dT eliminates the need to correct
            % endpoint transients resulting from phase shifts across
            % discontinuous data. It also doesn't require low-pass
            % filtering, I think for a similar reason. However, it comes at
            % the cost of random temperature offsets introduced due to the
            % variance of the first temperature data point, used for
            % reintegration. 
            
            % Therefore, for response matching after dehysteresis, I
            % recoomend using the dT correction method. For response
            % matching before hysteresis, use the T correction with
            % low-pass filtering and transient cutoffs.
            
            % After some discussion, keeping a low-pass filter is still
            % good becasuse it addresses the issue of decreasing coherence
            % between temperature and conductivity at high frequencies. If
            % you're trying to apply a correction to match one to another
            % but aren't sure how the two are related you're creating bogus
            % data anyway. It is also adviseable to keep the end-point
            % transients since the data isn't periodic. 
            
            % Also, maybe I want to detrend anyway before taking the
            % spectrum so that I'm as close to being centered around zero
            % as possible. But near the top and bottom of the profile, the
            % instrument is moving slowly, so dT/dt is close to zero
            % already. Removing any sort of mean or trend would actually
            % shift the end points away from zero, which is something I'd
            % like to avoid.
            
            dTdt = diff(T);
            [f,t]=ds_fft_AAA(dTdt,dt);
            
            % 2. Fit the Gain and Phase
            
            GainFit=polyval(Spectra_mean.(profiletype).(direction).Again{Iserial},f);
            PhaseFit=polyval(Spectra_mean.(profiletype).(direction).Aphase{Iserial},f);
            
            % Normalize the gain so that mean Temp gradient is unaltered 
            GainFit=GainFit./GainFit(1);
            
            % Reverse sign of negative freq phases
            [~,Ineg]=ds_f_ind_AAA(f);
            PhaseFit(Ineg)=-PhaseFit(Ineg);
            
            % 3. Correct temperature gradients
            t=t.*GainFit'.*exp(1i.*PhaseFit');
            
            % 4. Transform back into the time domain
            
            dT_corrected = real(ifft(t));
            Tcorrected = cumsum([T(1); dT_corrected]);
            
            % 5. Low-pass filter
            % Rob used a filter at 1/3 the nyquist frequency, I'll keep
            % that for now. Seems to be about where T-C coherence <0.9
            
            [b,a] = butter(4,1/3);
            Tcorrected = filtfilt(b,a,Tcorrected);
            Ccorrected = filtfilt(b,a,C);
            Pcorrected = filtfilt(b,a,P);
            
            % End point transient. When the timeseries is not zero at both
            % ends, changing the phase in fourier space causes large
            % amplitude changes at the ends of the data. This effect is
            % reduced by detrending the timeseries, and reduced further by
            % adjusting dT instead of T. It is best to cut off all these
            % "transient" points in the final data. Here, I find what the
            % maximum transient effect is by converting the phase shift
            % at each frequency into a time shift, and thus an index shift,
            % and save the maximum number of points. The cutting is done
            % later once the CTD data is combined with auxiliary data so
            % that everything is cut to the same length.
            
            tau = PhaseFit./(f); % Timeshift at each frequency
            tauind = tau/dt; % Indicie shift at each frequency;
            transient = ceil(max(abs(tauind(2:end))));
                        
            FCTDmatched.(profiletype).(direction).temperature{k} = Tcorrected;
            FCTDmatched.(profiletype).(direction).conductivity{k} = Ccorrected;
            FCTDmatched.(profiletype).(direction).pressure{k} = Pcorrected;
            FCTDmatched.(profiletype).(direction).transient{k} = transient;
            
        end
    end
end

save(datapath,'-struct','FCTDmatched','-v7.3');
disp(['Saved Response-Matched FCTD profiles in ' datapath]);

%% Make some plots
if plotevery
for k=1:plotevery:nprof
figure
clf;
    ax(1) = subplot(2,2,1);
        plot(FCTDprofiles.down.temperature{k},-FCTDprofiles.down.pressure{k});
        hold on;
        plot(FCTDdehyst.down.temperature{k},-FCTDdehyst.down.pressure{k});
        plot(FCTDmatched.dehyst.down.temperature{k},-FCTDmatched.dehyst.down.pressure{k});
        xlabel('Temperature');
        ylabel('Pressure');
        title(['Temperature Profile Down: Profile ' num2str(k)]);
        legend('Raw','deHyst','Matched','Location','best');
        grid on;
        
    ax(2) = subplot(2,2,2);
        plot(FCTDprofiles.up.temperature{k},-FCTDprofiles.up.pressure{k});
        hold on;
        plot(FCTDdehyst.up.temperature{k},-FCTDdehyst.up.pressure{k});
        plot(FCTDmatched.dehyst.up.temperature{k},-FCTDmatched.dehyst.up.pressure{k});
        xlabel('Temperature');
        ylabel('Pressure');
        title(['Temperature Profile Up: Profile ' num2str(k)]);
        legend('Raw','deHyst','Matched','Location','best');
        grid on;        

   ax(3) = subplot(2,2,3);
        plot(FCTDprofiles.down.conductivity{k},-FCTDprofiles.down.pressure{k});
        hold on;
        plot(FCTDdehyst.down.conductivity{k},-FCTDdehyst.down.pressure{k});
        plot(FCTDmatched.dehyst.down.conductivity{k},-FCTDmatched.dehyst.down.pressure{k});
        xlabel('Conductivity');
        ylabel('Pressure');
        title(['Conductivity Profile Down: Profile ' num2str(k)]);
        legend('Raw','deHyst','Matched','Location','best');
        grid on;
        
   ax(4) = subplot(2,2,4);
        plot(FCTDprofiles.up.conductivity{k},-FCTDprofiles.up.pressure{k});
        hold on;
        plot(FCTDdehyst.up.conductivity{k},-FCTDdehyst.up.pressure{k});
        plot(FCTDmatched.dehyst.up.conductivity{k},-FCTDmatched.dehyst.up.pressure{k});
        xlabel('Conductivity');
        ylabel('Pressure');
        title(['Conductivity Profile Up: Profile ' num2str(k)]);
        legend('Raw','deHyst','Matched','Location','best');
        grid on;        

        linkaxes(ax(1:2));
       linkaxes(ax(3:4));

    if tosave
        fname = ['matching_TC_P_profile_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
    pause(0.5);
       
figure
clf;
    ax(1) = subplot(2,2,1);
        plot(FCTDprofiles.up.temperature{k},-FCTDprofiles.up.pressure{k});
        hold on;
        plot(FCTDprofiles.down.temperature{k},-FCTDprofiles.down.pressure{k});
        xlabel('Temperature');
        ylabel('Pressure');
        title(['Temperature Profile Raw: Profile ' num2str(k)]);
        legend('Up','Down','Location','best');
        grid on;
        
    ax(2) = subplot(2,2,2);
        plot(FCTDmatched.dehyst.up.temperature{k},-FCTDmatched.dehyst.up.pressure{k});
        hold on;
        plot(FCTDmatched.dehyst.down.temperature{k},-FCTDmatched.dehyst.down.pressure{k});
        xlabel('Temperature');
        ylabel('Pressure');
        title(['Temperature Profile Corrected: Profile ' num2str(k)]);
        legend('Up','Down','Location','best');
        grid on;
        
   ax(3) = subplot(2,2,3);
        plot(FCTDprofiles.up.conductivity{k},-FCTDprofiles.up.pressure{k});
        hold on;
        plot(FCTDprofiles.down.conductivity{k},-FCTDprofiles.down.pressure{k});
        xlabel('Conductivity');
        ylabel('Pressure');
        title(['Conductivity Profile Raw: Profile ' num2str(k)]);
        legend('Up','Down','Location','best');
        grid on;
        
   ax(4) = subplot(2,2,4);
        plot(FCTDmatched.dehyst.up.conductivity{k},-FCTDmatched.dehyst.up.pressure{k});
        hold on;
        plot(FCTDmatched.dehyst.down.conductivity{k},-FCTDmatched.dehyst.down.pressure{k});
        xlabel('Conductivity');
        ylabel('Pressure');
        title(['Conductivity Profile Corrected: Profile ' num2str(k)]);
        legend('Up','Down','Location','best');
        grid on;      

        linkaxes(ax(1:2));
        linkaxes(ax(3:4));
       
    if tosave
        fname = ['matching_TC_P_updown_profile_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
    pause(0.5)
figure
clf;
    subplot(1,4,1)
        plot(FCTDprofiles.down.conductivity{k},FCTDprofiles.down.temperature{k});
        hold on;
        plot(FCTDprofiles.up.conductivity{k},FCTDprofiles.up.temperature{k});
        xlabel('Conductivity');
        ylabel('Temperautre');
        title(['T C Profile Raw: Profile ' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
        
    subplot(1,4,2)
        plot(FCTDdehyst.down.conductivity{k},FCTDdehyst.down.temperature{k});
        hold on;
        plot(FCTDdehyst.up.conductivity{k},FCTDdehyst.up.temperature{k});
        xlabel('Conductivity');
        ylabel('Temperautre');
        title(['T C Profile De-Hysteresis: Profile ' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;

    subplot(1,4,3)
        plot(FCTDmatched.raw.down.conductivity{k},FCTDmatched.raw.down.temperature{k});
        hold on;
        plot(FCTDmatched.raw.up.conductivity{k},FCTDmatched.raw.up.temperature{k});
        xlabel('Conductivity');
        ylabel('Temperautre');
        title(['T C Profile Raw, Response-Matched: Profile ' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
        
    subplot(1,4,4)
        plot(FCTDmatched.dehyst.down.conductivity{k},FCTDmatched.dehyst.down.temperature{k});
        hold on;
        plot(FCTDmatched.dehyst.up.conductivity{k},FCTDmatched.dehyst.up.temperature{k});
        xlabel('Conductivity');
        ylabel('Temperautre');
        title(['T C Profile De-Hysteresis and Response-Matched: Profile ' num2str(k)]);
        legend('Down','Up','Location','best');
        grid on;
 
        linkaxes;
        
    if tosave
        fname = ['matching_TC_profile_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
    pause(0.5);
end
end