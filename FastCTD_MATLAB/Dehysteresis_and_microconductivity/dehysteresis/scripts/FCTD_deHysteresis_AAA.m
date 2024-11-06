function [FCTDdehyst,Parameters]=FCTD_deHysteresis_AAA(FCTD,datapath,varargin)
% This code calculates the parameters to be used for the hysteresis
% calculation. It is a rewrite of Rob Pinkel's
% I_2016_D_deHysteresisFLEAT2016.m
%
% For now, while I figure out how to use gradient descent, the code samples
% the parameter space to determine the best combination of parameters
%
% Inputs:
%
%   Required:
%       FCTD: A strucure of qc'd profiles which is the output of
%           FCTD_CleanProfiles_AAA. 
%
%       datapath: A string with the filename for the output FCTDdehyst
%           structure, such as './Matching/FCTD_dehyst.mat'
%
%   Optional:
%       parameters: A matrix of parameter inputs to the dehystereis code,
%           usualy made from a previous run of the code. If not supplied
%           the code will re-calculate optimal parameters.
%               Example: 'parameters',Parameters_deHysteresis
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
%
% Outputs:
%
%   FCTDdehyst: A structure of profiles which are hysteresis-corrected
%
%   Parameters: A structure of numeric values as the parameters used in the
%    dehysteresis adjustment
%
% Examples:
%
%   [FCTDdehyst,Parameters]=FCTD_deHysteresis_AAA(FCTDprofiles_qc,'./Matching/FCTD_dehyst.mat')
%   FCTDdehyst=FCTD_deHysteresis_AAA(FCTDprofiles_qc,'./Matching/FCTD_dehyst.mat','parameters',Parameters_deHysteresis,'figpath','./figures/','plotevery',10)
%
% Alex Andriatis
% 2021-02-01

%Testing
% ccc;
% datapath = '/Volumes/Andriatis_T7BL/Data/TFO/TFO1/FCTD';
% FCTD = load(fullfile(datapath,'Matching','FCTDprofiles.mat'));
% Parameters.TParams=[4, 0.075; 4, 0.085];
% Parameters.serial=[536,537];
% Parameters.offset_T=1;
% Parameters.CParams=[0.0316, 4.00, 0.1269; 0.0316, 4.00, 0.1269];
% Parameters.offset_C=0;

P=inputParser;
addRequired(P,'FCTD',@isstruct);

validationFcn = @(x) validateattributes(x,{'char','string'},{});
addRequired(P,'datapath',validationFcn);

defaultParameters = [];
addParameter(P,'parameters',defaultParameters,@isstruct);

defaultFigpath = '';
validationFcn = @(x) validateattributes(x,{'char','string'},{});
addParameter(P,'figpath',defaultFigpath,validationFcn);

defaultPlotevery = 0;
addParameter(P,'plotevery',defaultPlotevery,@isnumeric);

parse(P,FCTD,datapath,varargin{:});
FCTD = P.Results.FCTD;
datapath = P.Results.datapath;
Parameters = P.Results.parameters;
figpath = P.Results.figpath;
plotevery = P.Results.plotevery;

% Define output structure and carry over auxilary variables.
% This is fine since I'm not changing the length of profiles here.
FCTDdehyst.header = FCTD.header;

% Turn on plotting
if ~isempty(figpath)
    tosave=1;
    fpath = fullfile(figpath,'deHysteresis');
    mkdir(fpath);
end

% Specify number of profiles to use in the optimization
nrand = 1000;


%% Properties of the data file:

reqnames = {'time','temperature','conductivity','pressure'};

directions = {'up','down'};

for i=1:2
    if ~all(isfield(FCTD,directions{i}))
        error(['The FCTD profile file does not have the ' directions{i} ' direction']);
    end
    if ~all(isfield(FCTD.(directions{i}),reqnames))
        error(['One or more required fields is missing from the ' directions{i} ' direction']);
    end
end

if ~isfield(FCTD,'header')
	error('The FCTD profile file does not have a header with the serial number')
end

nprof = length(FCTD.down.time);
fNy = round(1/(2*mean(diff(FCTD.down.time{1})*86400))); % Nyquist frequency of sample rate in cycles per second. Nomially 8 for 16 Hz sampling.

%% Separate the data by serial number of FASTCTD for the optimization
serials = unique(FCTD.header.serial);

%% 1. Correct the temperature cell relative to the pressure cell.
% Correct the temperatrue for the offset between pressure and temperature sensors
%
% Because the construction of the CTDs is the same, it's safe to treat the
% serial numbers the same.
%
% This code applies a time offset to correct the temperature, based on the
% difference in position between the temperature and pressure cells. The
% geometry of the fastctd changes between the way up and the way down, and
% aslo depends on the fall rate. A faster fall, such as at the top of the
% profile, will result in a larger pressure difference than towards the
% bottom of the profile. A positive offset means that the temperatrue cell
% is shallower than the pressure cell on the way down. It is assumed that
% the opposite orientation is true on the way up. If the angle of attack is
% steeper on the way down than the way up, the correction will apply a warm
% bias on the way down, and a cold bias on the way up. 
%
% For each profile, load the pressure and temperature
% Compute the cost of the mismatch:
%   1. Grid the up and down profiles on the same pressure grid
%   2. Compute the mean square error between up and down
% Try a variety of offset parameters from 0 to... 12? 4? 4 is probably
% enough. 4 = .25seconds * 4m/s = 1m for 16Hz sampling.
% Pick the offset that results in the least error
% Correct temperature using this offset measure


optimize=1;
if isfield(Parameters,'offset_T')
    if ~isnan(Parameters.offset_T)
        optimize=0;
    end
end
if optimize
    % Subsample the profiles to keep computation time reasonable -
    % use a random 500 profiles
    Isub = randi([1 nprof],[1,nrand]);
    Isub = unique(Isub);
    I_qc=FCTD.header.qc(Isub);
    Isub(I_qc~=0)=[];

    data=[];
    data.up.pressure = FCTD.up.pressure(Isub);
    data.up.var = FCTD.up.temperature(Isub);
    data.down.pressure = FCTD.down.pressure(Isub);
    data.down.var = FCTD.down.temperature(Isub);

    offsetrange = [-4:4];
    Jrecord=NaN(1,length(offsetrange));
    
    for i=1:length(offsetrange)
        Jrecord(i)=cost_Offset_P(data,offsetrange(i));
    end
    [~,I] = min(Jrecord);
    Parameters.offset_T=offsetrange(I);
    
    % Plot the offset cost function
    figure
        plot(offsetrange,Jrecord);
        xlabel('Offset T');
        ylabel('Up-Down T vs P error')
        title('Missmatch between up and down T-P profiles with T offset');
        hold on
        plot(offsetrange(I),Jrecord(I),'pm','MarkerFaceColor','m','MarkerSize',20);

    if tosave
        fname = ['deHyst_temperature_offset'];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
end

for n=1:2
    direction = directions{n};
    for k=1:nprof
        T = FCTD.(direction).temperature{k};
        offset = Parameters.offset_T;
        FCTDdehyst.(direction).temperature{k}=deHyst_Offset_P(T,offset);
    end
    FCTDdehyst.(direction).pressure=FCTD.(direction).pressure;
end


%% 2. Apply a temperature hysteresis correction.
% The corrected temperature at time t is given by a combination of the
% current temperature and the temperatrue at a prior time, with a weighting
% factor. The lag and weighting are found by minimizing the difference
% between up and downcasts in dT/dz vs T space.

% "Dumb" parameter space optimization of temperature hysteresis correction
%
% After many attempts at making gradient decent work, using both Newton's
% method and vanilla gradient descent, I give up. The parameters are too
% poorly conditioned for this to work reliably over a range of datasets,
% requiring hand-picked parameters to run the gradient descent code.
%
% Instead, let's just loop over the parameter space and do a ton of
% computations to find approximately a good parameter
%
% For temperature, we'll loop over different value of lagg and alpha.

for i=1:length(serials)
    optimize=1;
    if isfield(Parameters,'TParams') & isfield(Parameters,'serial')
        I=find(Parameters.serial == serials(i));
        if ~isempty(I)
            if I<=size(Parameters.TParams,1)
                if all(~isnan(Parameters.TParams(I,:)))
                    optimize=0;
                end
            end
        end
    end
    if optimize

        % Number of sample parameters
        numel=50;

        % A range of possible lags
        laggs = [1:15];

        % A range of possible time constants
        taulim = [0.01,0.2];
        taus = logspace(log10(taulim(1)),log10(taulim(2)),numel);

        disp(['Optimizing parameters for CTD serial number ' num2str(serials(i))]);

        Iserial = find(FCTD.header.serial == serials(i) & FCTD.header.qc==0);
        
        vars = {'pressure','temperature'};
        data=[];
        for n=1:length(directions)
            direction = directions{n};
            for m=1:length(vars)
                var = vars{m};
                data.(direction).(var)=FCTDdehyst.(direction).(var)(Iserial);
            end
        end
        
        % Subsample the profiles to keep computation time reasonable -
        % use a random 500 profiles
        nsub = length(data.up.pressure);
        Isub = randi([1 nsub],[1,nrand]);
        Isub = unique(Isub);

        vars = {'pressure','temperature'};
        for n=1:length(directions)
            direction = directions{n};
            for m=1:length(vars)
                var = vars{m};
                data.(direction).(var)=data.(direction).(var)(Isub);
            end
        end

        costs = NaN(length(laggs),length(taus));
        for j=1:length(laggs)
            TParams(1) = laggs(j);
            disp(['Trying lagg = ' num2str(TParams(1))]);
            for k=1:length(taus)
                TParams(2) = taus(k);
                %disp(['Trying tau = ' num2str(TParams(2))]);
                costs(j,k)=cost_T_serial(data,TParams);
            end
        end

        Jmin = min(costs,[],'all');
        [I,J]=find(costs==Jmin,1);
        lag = laggs(I);
        tau = taus(J);
        Parameters.TParams(i,:)=[lag, tau];
        Parameters.serial(i) = serials(i);

        figure
            [C,h]=contourf(laggs,taus,log10(costs'),linspace(min(log10(costs),[],'all'),max(log10(costs),[],'all'),100));
            h.LineColor='none';
            hold on;
            plot(laggs(I),taus(J),'pm','MarkerFaceColor','m','MarkerSize',20);
            xlabel('lag');
            ylabel('tau');
            title(['Temperature hysteresis parameter space for CTD Serial ' num2str(serials(i))]);
            c = colorbar;
            c.Label.String = 'Up-Down RMSE';
            
            if tosave
                fname = ['deHyst_temperature_parameters_' num2str(serials(i))];
                saveas(gcf,fullfile(fpath,fname),'fig');
                saveas(gcf,fullfile(fpath,fname),'png');
            end

        lag = Parameters.TParams(i,1);
        tau = Parameters.TParams(i,2);

        laggs = [max(lag-2,1):lag+2];
        taus = linspace(0.8*tau,1.2*tau,numel);

        costs = NaN(length(laggs),length(taus));
        for j=1:length(laggs)
             TParams(1) = laggs(j);
             disp(['Trying lagg = ' num2str(TParams(1))]);
            for k=1:length(taus)
                TParams(2) = taus(k);
                %disp(['Trying tau = ' num2str(TParams(2))]);
                costs(j,k)=cost_T_serial(data,TParams);
            end
        end
        Jmin = min(costs,[],'all');
        [I,J]=find(costs==Jmin,1);
        lag = laggs(I);
        tau = taus(J);
        Parameters.TParams(i,:)=[lag, tau];

        figure
            [C,h]=contourf(laggs,taus,log10(costs'),linspace(min(log10(costs),[],'all'),max(log10(costs),[],'all'),100));
            h.LineColor='none';
            hold on;
            plot(laggs(I),taus(J),'pm','MarkerFaceColor','m','MarkerSize',20);
            xlabel('lag');
            ylabel('tau');
            title(['Temperature hysteresis parameter space for CTD Serial ' num2str(serials(i))]);
            c = colorbar;
            c.Label.String = 'Up-Down RMSE';
            
            if tosave
                fname = ['deHyst_temperature_parameters_zoom_' num2str(serials(i))];
                saveas(gcf,fullfile(fpath,fname),'fig');
                saveas(gcf,fullfile(fpath,fname),'png');
            end
    end
end

%% Correct temperature in the up and down casts

for n=1:2
    direction = directions{n};
    data = FCTDdehyst.(direction);
    for i=1:nprof
        Iserial = find(serials==FCTD.header.serial(i));
      
        T = data.temperature{i};
        TParams = Parameters.TParams(Iserial,:);
        
        That = deHyst_T(T,TParams);
        
        FCTDdehyst.(direction).temperature{i} = That;
    end
end


%% 3. Correct the conductivity cell relative to the pressure cell.
% Just as with the spatial offset between the pressure and temperature
% cells, the conductivity cell is offset from the temperature and pressure
% cells.  Since we adjust temp relative to pressure, adjust the
% conductivity relative to pressure too. 

optimize=1;
if isfield(Parameters,'offset_C')
    if ~isnan(Parameters.offset_C)
        optimize=0;
    end
end
if optimize
    % Subsample the profiles to keep computation time reasonable -
    % use a random 500 profiles
    Isub = randi([1 nprof],[1,nrand]);
    Isub = unique(Isub);
    I_qc=FCTD.header.qc(Isub);
    Isub(I_qc~=0)=[];

    data=[];
    data.up.pressure = FCTD.up.pressure(Isub);
    data.up.var = FCTD.up.conductivity(Isub);
    data.down.pressure = FCTD.down.pressure(Isub);
    data.down.var = FCTD.down.conductivity(Isub);

    offsetrange = [-4:4];
    Jrecord=NaN(1,length(offsetrange));
    
    for i=1:length(offsetrange)
        Jrecord(i)=cost_Offset_P(data,offsetrange(i));
    end
    [~,I] = min(Jrecord);
    Parameters.offset_C=offsetrange(I);
    
    % Plot the offset cost function
    figure
        plot(offsetrange,Jrecord);
        xlabel('Offset C');
        ylabel('Up-Down C vs P error')
        title('Missmatch between up and down C-P profiles with C offset');
        hold on
        plot(offsetrange(I),Jrecord(I),'pm','MarkerFaceColor','m','MarkerSize',20);

    if tosave
        fname = ['deHyst_conductivity_offset'];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
end

for n=1:2
    direction = directions{n};
    for k=1:nprof
        C = FCTD.(direction).conductivity{k};
        offset = Parameters.offset_C;
        FCTDdehyst.(direction).conductivity{k}=deHyst_Offset_P(C,offset);
    end
end



%% 4. Apply a conductivity hysteresis correction.
% The corrected conductivity at time t is given the conductivity at time C
% plus a correction factor. The correction factor at time t is given some
% factor of the conductivity at the previous time plus a factor based on
% the temperature at the previous time. Factors are chosen that minimize
% the difference between the up and downcasts in dC/dz vs C space.

for i=1:length(serials)
    optimize=1;
    if isfield(Parameters,'CParams')
        I=find(Parameters.serial == serials(i));
        if ~isempty(I)
            if I<=size(Parameters.CParams,1)
                if all(~isnan(Parameters.CParams(I,:)))
                    optimize=0;
                end
            end
        end
    end
    if optimize

        % The range of parameters
        numel=10;

        gammalim = [0.01 0.1];
        gammas = logspace(log10(gammalim(1)),log10(gammalim(2)),numel);

        TimeCs = linspace(0,20,numel);

        alphalim = [0.001 0.5];
        alphas = logspace(log10(alphalim(1)),log10(alphalim(2)),numel);
        
        serial = serials(i);
        disp(['Optimizing Conductivity parameters for CTD serial number ' num2str(serial)]);
        Iserial = find(FCTD.header.serial == serial  & FCTD.header.qc==0);
        
        vars = {'time','pressure','conductivity','temperature'};
        data=[];
        for n=1:length(directions)
            direction = directions{n};
            for m=1:length(vars)
                var = vars{m};
                data.(direction).(var)=FCTDdehyst.(direction).(var)(Iserial);
            end
        end
        
        % Subsample the profiles to keep computation time reasonable -
        % use a random 500 profiles
        nsub = length(data.up.pressure);
        Isub = randi([1 nsub],[1,nrand]);
        Isub = unique(Isub);

        vars = {'time','pressure','temperature','conductivity'};
        for n=1:length(directions)
            direction = directions{n};
            for m=1:length(vars)
                var = vars{m};
                data.(direction).(var)=data.(direction).(var)(Isub);
            end
        end

        lag = Parameters.TParams(i,1);

        costs = NaN(length(gammas),length(TimeCs),length(alphas));
        for j=1:length(gammas)
            disp(['Trying gamma = ' num2str(gammas(j))]);
            for k=1:length(TimeCs)
                for l=1:length(alphas)
                    aC=4*fNy*alphas(l)*TimeCs(k)/(1+4*fNy*TimeCs(k)); 
                    bC=1-2*aC/alphas(l);
                    costs(j,k,l) = cost_C_serial(data,aC,bC,gammas(j),lag);
                end
            end
        end

        Jmin = min(costs,[],'all');
        [I,J,K]=ind2sub(size(costs),find(costs==Jmin,1));
        gamma = gammas(I);
        TimeC = TimeCs(J);
        alpha = alphas(K);
        Parameters.CParams(i,:)=[gamma, TimeC, alpha];

        figure
            subplot(3,1,1)
                z = log10(squeeze(costs(:,:,K))');
                [C,h]=contourf(gammas,TimeCs,z,linspace(min(z,[],'all'),max(z,[],'all'),100));
                h.LineColor='none';
                hold on;
                plot(gammas(I),TimeCs(J),'pm','MarkerFaceColor','m','MarkerSize',20);
                xlabel('gamma');
                ylabel('TimeC');
                title(['Conductivity hysteresis parameter space for CTD Serial ' num2str(serials(i))]);
                c = colorbar; 
                c.Label.String = 'Up-Down RMSE';

            subplot(3,1,2)
                z = log10(squeeze(costs(:,J,:))');
                [C,h]=contourf(gammas,alphas,z,linspace(min(z,[],'all'),max(z,[],'all'),100));
                h.LineColor='none';
                hold on;
                plot(gammas(I),alphas(K),'pm','MarkerFaceColor','m','MarkerSize',20);
                xlabel('gamma');
                ylabel('alpha');
                title(['Conductivity hysteresis parameter space for CTD Serial ' num2str(serials(i))]);
                c = colorbar;
                c.Label.String = 'Up-Down RMSE';     

             subplot(3,1,3)
                z = log10(squeeze(costs(I,:,:))');
                [C,h]=contourf(TimeCs,alphas,z,linspace(min(z,[],'all'),max(z,[],'all'),100));
                h.LineColor='none';
                hold on;
                plot(TimeCs(J),alphas(K),'pm','MarkerFaceColor','m','MarkerSize',20);
                xlabel('TimeC');
                ylabel('alpha');
                title(['Conductivity hysteresis parameter space for CTD Serial ' num2str(serials(i))]);
                c = colorbar;
                c.Label.String = 'Up-Down RMSE';  

                
            if tosave
                fname = ['deHyst_conductivity_parameters_' num2str(serials(i))];
                saveas(gcf,fullfile(fpath,fname),'fig');
                saveas(gcf,fullfile(fpath,fname),'png');
            end

        gamma = Parameters.CParams(i,1);
        TimeC = Parameters.CParams(i,2);
        alpha = Parameters.CParams(i,3);

        gammas = linspace(0.8*gamma,1.2*gamma,numel);
        TimeCs = linspace(0.8*TimeC,1.2*TimeC,numel);
        alphas = linspace(0.8*alpha,1.2*alpha,numel);

        lag = Parameters.TParams(i,1);

        costs = NaN(length(gammas),length(TimeCs),length(alphas));
        for j=1:length(gammas)
            disp(['Trying gamma = ' num2str(gammas(j))]);
            for k=1:length(TimeCs)
                for l=1:length(alphas)
                    aC=4*fNy*alphas(l)*TimeCs(k)/(1+4*fNy*TimeCs(k)); 
                    bC=1-2*aC/alphas(l);
                    costs(j,k,l) = cost_C_serial(data,aC,bC,gammas(j),lag);
                end
            end
        end

        Jmin = min(costs,[],'all');
        [I,J,K]=ind2sub(size(costs),find(costs==Jmin,1));
        gamma = gammas(I);
        TimeC = TimeCs(J);
        alpha = alphas(K);
        Parameters.CParams(i,:)=[gamma, TimeC, alpha];

        figure
            subplot(3,1,1)
                z = log10(squeeze(costs(:,:,K))');
                [C,h]=contourf(gammas,TimeCs,z,linspace(min(z,[],'all'),max(z,[],'all'),100));
                h.LineColor='none';
                hold on;
                plot(gammas(I),TimeCs(J),'pm','MarkerFaceColor','m','MarkerSize',20);
                xlabel('gamma');
                ylabel('TimeC');
                title(['Conductivity hysteresis parameter space for CTD Serial ' num2str(serials(i))]);
                c = colorbar;
                c.Label.String = 'Up-Down RMSE';

            subplot(3,1,2)
                z = log10(squeeze(costs(:,J,:))');
                [C,h]=contourf(gammas,alphas,z,linspace(min(z,[],'all'),max(z,[],'all'),100));
                h.LineColor='none';
                hold on;
                plot(gammas(I),alphas(K),'pm','MarkerFaceColor','m','MarkerSize',20);
                xlabel('gamma');
                ylabel('alpha');
                title(['Conductivity hysteresis parameter space for CTD Serial ' num2str(serials(i))]);
                c = colorbar;
                c.Label.String = 'Up-Down RMSE';     

             subplot(3,1,3)
                z = log10(squeeze(costs(I,:,:))');
                [C,h]=contourf(TimeCs,alphas,z,linspace(min(z,[],'all'),max(z,[],'all'),100));
                h.LineColor='none';
                hold on;
                plot(TimeCs(J),alphas(K),'pm','MarkerFaceColor','m','MarkerSize',20);
                xlabel('TimeC');
                ylabel('alpha');
                title(['Conductivity hysteresis parameter space for CTD Serial ' num2str(serials(i))]);
                c = colorbar;
                c.Label.String = 'Up-Down RMSE';
                
            if tosave
                fname = ['deHyst_conductivity_parameters_zoom_' num2str(serials(i))];
                saveas(gcf,fullfile(fpath,fname),'fig');
                saveas(gcf,fullfile(fpath,fname),'png');
            end
    end
end



%% Conductivity correction
for n=1:2
    direction = directions{n};
    data = FCTDdehyst.(direction);
    
    for k=1:nprof
        Iserial = find(serials==FCTD.header.serial(k));
        lag = Parameters.TParams(Iserial,1);
        gamma = Parameters.CParams(Iserial,1);
        TimeC = Parameters.CParams(Iserial,2);
        alpha = Parameters.CParams(Iserial,3);
        
        aC=4*fNy*alpha*TimeC/(1+4*fNy*TimeC); 
        bC=1-2*aC/alpha;

        C = data.conductivity{k};
        T = data.temperature{k};
    
        Ccorrected = deHyst_C(T,C,aC,bC,gamma,lag);
        
        FCTDdehyst.(direction).conductivity{k} = Ccorrected;
    end
end    

%% Save the corrected profles

save(datapath,'-struct','FCTDdehyst','-v7.3');
disp(['Saved hysteresis-corrected profiles in ' datapath]);

%% Plot diagnostics
if plotevery %switch to turn on plotting
%% Compare raw up and down casts with corrected up and down casts

for k=1:plotevery:nprof

Pup = FCTD.up.pressure{k};
Cup = FCTD.up.conductivity{k};
Cup_dehyst = FCTDdehyst.up.conductivity{k};
Tup = FCTD.up.temperature{k};
Tup_dehyst = FCTDdehyst.up.temperature{k};


Pdown = FCTD.down.pressure{k};
Cdown = FCTD.down.conductivity{k};
Cdown_dehyst = FCTDdehyst.down.conductivity{k};
Tdown = FCTD.down.temperature{k};
Tdown_dehyst = FCTDdehyst.down.temperature{k};

dCdzup = diff(Cup)./diff(Pup);
Cupmid = Cup(2:end); % Temperatrue at the top of a dT/dz step
dCdzdown = diff(Cdown)./diff(Pdown);
Cdownmid = Cdown(1:end-1); % Temperature at the top of a dT/dz step

dTdzup = diff(Tup)./diff(Pup);
Tupmid = Tup(2:end); % Temperatrue at the top of a dT/dz step
dTdzdown = diff(Tdown)./diff(Pdown);
Tdownmid = Tdown(1:end-1); % Temperature at the top of a dT/dz step

dCdzup_dehyst = diff(Cup_dehyst)./diff(Pup);
Cupmid_dehyst = Cup_dehyst(2:end); % Temperatrue at the top of a dT/dz step
dCdzdown_dehyst = diff(Cdown_dehyst)./diff(Pdown);
Cdownmid_dehyst = Cdown_dehyst(1:end-1); % Temperature at the top of a dT/dz step

dTdzup_dehyst = diff(Tup_dehyst)./diff(Pup);
Tupmid_dehyst = Tup_dehyst(2:end); % Temperatrue at the top of a dT/dz step
dTdzdown_dehyst = diff(Tdown_dehyst)./diff(Pdown);
Tdownmid_dehyst = Tdown_dehyst(1:end-1); % Temperature at the top of a dT/dz step

figure
    subplot(2,1,1)
    plot(Tupmid,dTdzup); hold on;
    plot(Tdownmid,dTdzdown);
    title(['dT/dz vs T Raw Profile ' num2str(k)]);
    legend('Up','Down');
    ylabel('dT/dz');
    xlabel('T');
    grid on;    
    
    subplot(2,1,2)
    plot(Tupmid_dehyst,dTdzup_dehyst); hold on;
    plot(Tdownmid_dehyst,dTdzdown_dehyst);
    title(['dT/dz vs T deHyst Profile ' num2str(k)]);
    legend('Up','Down');
    ylabel('dT/dz');
    xlabel('T');
    grid on;    

    linkaxes;
    
    if tosave
        fname = ['deHyst_dTdz_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'png');
    end

figure
    subplot(2,1,1)
    plot(Cupmid,dCdzup); hold on;
    plot(Cdownmid,dCdzdown);
    title(['dC/dz vs C Raw Profile ' num2str(k)]);
    legend('Up','Down');
    ylabel('dC/dz');
    xlabel('C');
    grid on;    
    
    subplot(2,1,2)
    plot(Cupmid_dehyst,dCdzup_dehyst); hold on;
    plot(Cdownmid_dehyst,dCdzdown_dehyst);
    title(['dC/dz vs C deHyst Profile ' num2str(k)]);
    legend('Up','Down');
    ylabel('dC/dz');
    xlabel('C');
    grid on;    

    linkaxes;
    
    if tosave
        fname = ['deHyst_dCdz_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
end


%% Plot up vs down profiles for T and C
for k=1:plotevery:nprof

Tup = FCTD.up.temperature{k};
Tup_dehyst = FCTDdehyst.up.temperature{k};
Pup = FCTD.up.pressure{k};

Tdown = FCTD.down.temperature{k};
Tdown_dehyst = FCTDdehyst.down.temperature{k};
Pdown = FCTD.down.pressure{k};


Cup = FCTD.up.conductivity{k};
Cup_dehyst = FCTDdehyst.up.conductivity{k};
Cdown = FCTD.down.conductivity{k};
Cdown_dehyst = FCTDdehyst.down.conductivity{k};

    
figure
    subplot(1,2,1)
    plot(Tup,Pup); hold on;
    plot(Tdown,Pdown);
    set(gca,'YDir','reverse');
    title(['T vs P Raw Profile ' num2str(k)]);
    legend('Up','Down');
    xlabel('T');
    ylabel('Pressure');    
    
    subplot(1,2,2)
    plot(Tup_dehyst,Pup); hold on;
    plot(Tdown_dehyst,Pdown);
    set(gca,'YDir','reverse');
    title(['T vs P deHyst Profile ' num2str(k)]);
    legend('Up','Down');
    xlabel('T');
    ylabel('Pressure');
    
    linkaxes;

    if tosave
        fname = ['deHyst_TvP_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end    
    
figure
    subplot(1,2,1)
    plot(Cup,Pup); hold on;
    plot(Cdown,Pdown);
    set(gca,'YDir','reverse');
    title(['C vs P Raw Profile ' num2str(k)]);
    legend('Up','Down');
    xlabel('C');
    ylabel('Pressure');    
    
    subplot(1,2,2)
    plot(Cup_dehyst,Pup); hold on;
    plot(Cdown_dehyst,Pdown);
    set(gca,'YDir','reverse');
    title(['C vs P deHyst Profile ' num2str(k)]);
    legend('Up','Down');
    xlabel('C');
    ylabel('Pressure');
    
    linkaxes;
    
    if tosave
        fname = ['deHyst_CvP_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end 
end

%% Plot up vs down profiles for T - C
for k=1:plotevery:nprof
    

Tup = FCTD.up.temperature{k};
Tup_dehyst = FCTDdehyst.up.temperature{k};
Pup = FCTD.up.pressure{k};

Tdown = FCTD.down.temperature{k};
Tdown_dehyst = FCTDdehyst.down.temperature{k};
Pdown = FCTD.down.pressure{k};


Cup = FCTD.up.conductivity{k};
Cup_dehyst = FCTDdehyst.up.conductivity{k};
Cdown = FCTD.down.conductivity{k};
Cdown_dehyst = FCTDdehyst.down.conductivity{k};

    
figure
    clf
    subplot(1,2,1)
    plot(Cup,Tup); hold on;
    plot(Cdown,Tdown);
    title(['T-C Raw Profile ' num2str(k)]);
    legend('Up','Down');
    xlabel('C');
    ylabel('T'); 
    
    subplot(1,2,2)
    plot(Cup_dehyst,Tup_dehyst); hold on;
    plot(Cdown_dehyst,Tdown_dehyst);
    title(['T-C DeHyst Profile ' num2str(k)]);
    legend('Up','Down');
    xlabel('C');
    ylabel('T'); 

    linkaxes;
    
    if tosave
        fname = ['deHyst_TvC_' num2str(k)];
        saveas(gcf,fullfile(fpath,fname),'fig');
        saveas(gcf,fullfile(fpath,fname),'png');
    end
end
end
