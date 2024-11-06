function DataGrid = FastCTD_GridData(FCTD,varargin)
% Version of FastCTD_GridData used on TFO Seamounts. Adds uConductivity and
% fluorometer data (lines 255- )
%
%  DataGrid = FCTDMakeGrid(FCTD);
%
%  INPUTS:
%       FCTD - structure with timeseries of time, pressure, temperature,
%       etc
%
%   Temperature, Conductivity, Pressure, Salinity and Potential Density
%   will be placed on a grid of temperature versus depth
%
%  DataGrid = FCTDMakeGrid(FCTD,downcast,todo);
%   in addition to the variables described above, one can specify other
%   variable by defining TODO using a cell structure of strings
%
%  DataGrid is a structure
%
% Written by Jody Klymak
% Updated 2011 07 14 by San Nguyen
% used 20 point median filter to smooth the pressure

%26 Aug 2024: MHA has redone the short term response matching code.

if ~isfield(FCTD,'pressure')
    disp(FCTD);
    DataGrid=[];
    return;
end;

vars2Grid = {'pressure','temperature','conductivity','altDist','fluor','drop','longitude','latitude','chi','chi2','w'};
vars2Grid_default = {'pressure','temperature','conductivity','density','salinity'};
pickDownCast = true;
zInterval = 0.5;
zMin = 0;
zMax = 2000;

persistent argsNameToCheck;
if isempty(argsNameToCheck);
    argsNameToCheck = {'VarsToGrid','upcast','downcast','zMin','zMax','zInterval'};
end

index = 1;
n_items = nargin-1;
while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        error('MATLAB:FastCTD_GridData:wrongOption','Incorrect option specified: %s',varargin{index});
    end
    
    switch i
        case 1 % varsToGrid
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            vars2Grid = varargin{index+1};
            if iscellstr(vars2Grid)
                for i = 1:length(vars2Grid)
                    if ~isfield(FCTD,vars2Grid{i})
                        error('MATLAB:FastCTD_GridData:wrongVar2Grid','Wrong variable to grid: %s', vars2Grid{i});
                    else
                        switch lower(vars2Grid{i})
                            case vars2Grid_default
                                vars2Grid{i} = lower(vars2Grid{i});
                                continue;
                            otherwise
                                error('MATLAB:FastCTD_GridData:wrongVar2Grid','Wrong variable to grid: %s', vars2Grid{i});
                        end
                    end
                end
            elseif ischar(vars2Grid)
                switch lower(vars2Grid)
                    case vars2Grid_default
                        vars2Grid = {lower(vars2Grid)};
                        continue;
                    otherwise
                        error('MATLAB:FastCTD_GridData:wrongVar2Grid','Wrong variable to grid: %s', vars2Grid);
                end
            else
                error('MATLAB:FastCTD_GridData:wrongVar2Grid','Variable to grid must be specified as a cell strings of variables: pressure, temperature, conductivity');
            end
            
            index = index +2;
            n_items = n_items-2;
        case 2 % upcast
            pickDownCast = false;
            
            index = index + 1;
            n_items = n_items - 1;
        case 3 % downcast
            pickDownCast = true;
            
            index = index + 1;
            n_items = n_items - 1;
        case 4 % zMin
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            zMin = varargin{index+1};
            if ~isnumeric(zMin)
                error('MATLAB:FastCTD_GridData:zMinNumeric','zMin must be numeric');
            elseif length(zMin) > 1 || isempty(zMin)
                error('MATLAB:FastCTD_GridData:zMinScalar','zMin must be a scalar');
            end
            index = index +2;
            n_items = n_items-2;
        case 5 % zMax
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            zMax = varargin{index+1};
            if ~isnumeric(zMax)
                error('MATLAB:FastCTD_GridData:zMaxNumeric','zMin must be numeric');
            elseif length(zMax) > 1 || isempty(zMax)
                error('MATLAB:FastCTD_GridData:zMinScalar','zMax must be a scalar');
            end
            index = index +2;
            n_items = n_items-2;
        case 6 % zInterval
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            zInterval = varargin{index+1};
            if ~isnumeric(zInterval)
                error('MATLAB:FastCTD_GridData:zMaxNumeric','zMin must be numeric');
            elseif length(zInterval) > 1 || isempty(zInterval)
                error('MATLAB:FastCTD_GridData:zMinScalar','zMax must be a scalar');
            elseif zInterval <= 0
                error('MATLAB:FastCTD_GridData:zMinZero','zInterval must be greater than zero');
            end
            index = index +2;
            n_items = n_items-2;
    end
end

if zMin > zMax
    zTemp = zMin;
    zMin = zMax;
    zMax = zTemp;
end
clear zTemp;

if zInterval > (zMax-zMin)
    error('MATLAB:FastCTD_GridData:zIntervalTooBig','zInterval must less than the range of z');
end

if pickDownCast
    FCTD = FastCTD_FindCasts(FCTD);
else
    FCTD = FastCTD_FindCasts(FCTD,'upcast');
end

if ~isfield(FCTD,'drop');
    DataGrid = [];
    return;
end;

zMin = zMin - zInterval/2;
zMax = zMax + zInterval/2;
DataGrid.depth=(zMin:zInterval:zMax)';



% % time correction for the Conductivity cell
% df = 2.2;
% dpha = 1.05;
%
% dt = dpha/df/2/pi;
% t = FCTD.time-dt/24/3600;
% conductivity = FCTD.conductivity;
%
% good = find(FCTD.time>0);
% if isempty(good);
%     DataGrid = [];
%     return;
% end;
%
% [it,ind] = unique(FCTD.time(good));
% good = good(ind);
% if length(good) > 3
%     conductivity(good) = interp1(FCTD.time(good),FCTD.conductivity(good),t(good));
% end
% FCTD.conductivity = conductivity;

drops = unique(FCTD.drop);
drops = drops(drops>0);

num = 0;
for i=1:length(drops);
    ind = find(FCTD.drop==drops(i));
    if max(FCTD.pressure(ind))-min(FCTD.pressure(ind))>10
        num = num+1;
    end
end

% allocate space to grid data
DataGrid.time = NaN(1,num);
for i = 1:length(vars2Grid)
    DataGrid.(vars2Grid{i}) = NaN(length(DataGrid.depth)-1,num);
end

% Add longitude and latitude
if isfield(FCTD,'GPS')
    FCTD.longitude = FCTD.GPS.longitude;
    FCTD.latitude = FCTD.GPS.latitude;
end

% load correction factors for FCTD
FCTD_SalCorr = load('FCTD_SalinityCorrectionFactors_toCond.mat');
% FCTD.depth = sw_dpth(FCTD.pressure,median(FCTD.GPS.latitude));
% FCTD.depth = sw_dpth(FCTD.pressure,20);
% for i = 1:length(vars2Grid)
%     FCTD.(vars2Grid{i}) = nanmedfilt1(FCTD.(vars2Grid{i}),10);
% end
FCTD_SalCorr.GainPFit = FCTD_SalCorr.GainPFit_Dn;
FCTD_SalCorr.PhsPFit = FCTD_SalCorr.PhsPFit_Dn;

% Loop through the drops and concatenate them
num = 0;
for i=1:length(drops)
    ind = find(FCTD.drop==drops(i));
    if max(FCTD.pressure(ind))-min(FCTD.pressure(ind))>10
        for j=1:length(vars2Grid)
            %ALB add try/catch because it break in the fluor during TFO seamount 2023 
            try
                myFCTD.(vars2Grid{j}) = FCTD.(vars2Grid{j})(ind);
            catch
                myFCTD.(vars2Grid{j})=nan.*ind;
            end
            % ALB end of hack
        end;
        %         myFCTD.temperature = FCTD.temperature(ind);
        %         myFCTD.pressure = FCTD.pressure(ind);
        %         myFCTD.conductivity = FCTD.conductivity(ind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BEGIN MHA EDITS OF THIS FUNCTION 8/26/2024
use_old_code=1; %Use 1 here to use San's legacy response matching code. It must provide the corrected variables myFCTD.pressure, temperature, conductivity and depth.
if use_old_code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin legacy code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do salinity despiking corrections
        good_ind = ~isnan(myFCTD.pressure);
        npts = sum(good_ind);
        df = FCTD_SalCorr.f_Ny/floor(npts/2);
        % creating the frequency axis
        myFCTD.f = (0:npts-1)'*df;
        myFCTD.f(myFCTD.f>FCTD_SalCorr.f_Ny) = myFCTD.f(myFCTD.f>FCTD_SalCorr.f_Ny)-2*FCTD_SalCorr.f_Ny;
        FCTD_SalCorr.GainFit = polyval(FCTD_SalCorr.GainPFit,myFCTD.f);
        FCTD_SalCorr.GainFit = FCTD_SalCorr.GainFit/FCTD_SalCorr.GainFit(1); % need to normalize Gain
        FCTD_SalCorr.PhsFit = polyval(FCTD_SalCorr.PhsPFit,myFCTD.f);
        
        % FFT
        myFCTD.T = fft(myFCTD.temperature(good_ind),npts);
        myFCTD.C = fft(myFCTD.conductivity(good_ind),npts);
        myFCTD.P = fft(myFCTD.pressure(good_ind),npts);
        
        % correct the conductivity
        myFCTD.CCorr = myFCTD.C.*FCTD_SalCorr.GainFit.*exp(-1i*FCTD_SalCorr.PhsFit);
        % Low Pass filter
        myFCTD.TCorr = myFCTD.T.*FCTD_SalCorr.LPfilter(myFCTD.f);
        myFCTD.CCorr = myFCTD.CCorr.*FCTD_SalCorr.LPfilter(myFCTD.f);
        myFCTD.PCorr = myFCTD.P.*FCTD_SalCorr.LPfilter(myFCTD.f);
        
        % get back to physical units
        myFCTD.tCorr = real(ifft(myFCTD.TCorr));
        myFCTD.cCorr = real(ifft(myFCTD.CCorr));
        myFCTD.pCorr = real(ifft(myFCTD.PCorr));
        
        % patching up the pressure because the pressure doesn't have a
        % phaseshift
        myFCTD.pCorr(1:13) = NaN;%myFCTD.pressure(1:13);
        myFCTD.pCorr(end-12:end) = NaN;%myFCTD.pressure(end-12:end);
        
        myFCTD.pressure = myFCTD.pCorr;
        myFCTD.temperature = myFCTD.tCorr;
        myFCTD.conductivity = myFCTD.cCorr;
        myFCTD.depth = sw_dpth(myFCTD.pressure,20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%end legacy code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%begin new MHA code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%0: Set parameters
        %MHA low pass filter.  Default (aggressive) values are sharpness = 2, fc=4.
        sharpness=2; %defines how quickly the filter goes from 1 to 0 around fc.
        fc=2; %cutoff freq.  The filter is 0.5 at fc.

        use_phys=0; %1 for physical model; 0 for polyfit
        zs=1/16; %sample interval in s - fixed at 1/16 sec.
        zc=4; %cutoff interval in s
        n_extra=16*zc*2; %get a factor of two more than the filter length on either side.

        %1 to compute spectra and polyval parameters on the fly; 0 to skip
        %and use defaults.
        compute_spectra=0;

        %Set to zero to ignore thermal mass correction; 1 to perform it.
        do_thermal_mass_correction = 0;
        % do_thermal_mass_correction = 0;

        % Defaults from Lueck and Picklo are alpha = 0.02 and beta = 0.1.
        %These are the parameters from Lueck and Picklo. 
        CTpar.alfa = 0.02; % SWIMS was 0.023
        CTpar.beta = 0.10; % SWIMS was 1/6.5

        %physical model parameters
        L=-.25/16;
        tau=0.07;

        N=4; %order of polynomial fit
        %For the downcast I looked at in TFO NESMA RR2410, pv_gain was:
        %pv_gain = 1x5 double
        %   -0.0015    0.0277   -0.1520    0.0892    1.0000
        %And for the upcast they were very similar:
        %pv_gain = 1x5 double
        %   -0.0012    0.0241   -0.1450    0.1041    1.0000
        %For now these are hard coded in.  They can be changed and should be loaded
        %in!

        %These coefficients and alpha need to be computed for each
        %cruise/sn combo.  This routine will compute them if
        %compute_spectra is set to 1.

        %polyfit/polyval parameters for gain and phase (radians)
        pv_gain=[-0.0015    0.0277   -0.1520    0.0892    1.0000];
        pv_ph=[-0.0048    0.0661   -0.3121    0.8792         0];
        %alpha = dC/dT
        alpha=9.7;

%1. Compute spectra so we can evaluate alpha
        %Routine to compute t, c and tc spec of tdata and cdata
        nfft=256;
        samplefreq=16;    %1/dt;
        df=samplefreq/nfft; % elementary frequency bandwidth
        f=(df:df:df*nfft/2)'; % frequency vector for spectra

        %physical model
        H=1./ ( (1-2*pi*i*f*tau).*exp(i*2*pi*f*L));



        if compute_spectra
        i1=1:length(myFCTD.temperature);
        it=i1(500:end-500); %choose a random range right now - THIS WILL CHOKE SOMETIMES
        lag=0;
        tdata=detrend(myFCTD.temperature(it+lag));
        cdata=detrend(myFCTD.conductivity(it));

        [Pt,fdum] = pwelch(diffs(tdata),nfft,[],nfft,samplefreq,'psd'); %Set NFFT equal to the window length.  Verified that kdum is equal to k computed above.
        Pt=Pt.';
        Pt=Pt(2:length(Pt)); %delete f=0; normalization is already correct.

        [Pc,fdum] = pwelch(diffs(cdata),nfft,[],nfft,samplefreq,'psd'); %Set NFFT equal to the window length.  Verified that kdum is equal to k computed above.
        Pc=Pc.';
        Pc=Pc(2:length(Pc)); %delete f=0; normalization is already correct.

        %Now compute cross spectrum with same dofs etc.
        [Pct,fdum] = cpsd(diffs(cdata),diffs(tdata),nfft,[],nfft,samplefreq,'psd'); %Set NFFT equal to the window length.  Verified that kdum is equal to k computed above.
        Pct=Pct.';
        Pct=Pct(2:length(Pct)); %delete f=0; normalization is already correct.
        %Maybe this is the best way to compute alpha...
        if1=find(f < 1);
        alpha=sqrt(nanmean(Pt(if1)) ./ nanmean(Pc(if1)));

        %Next I'll try a polyfit
        x=f;
        y=abs(Pct)./Pc/alpha;
                
        pv_gain=polyfit(x,y,N);
        %force low-freq val to be 1
        pv_gain(end)=1;
        
        % Phase
        x=f;
        y=atan2(imag(Pct),real(Pct));
                
        pv_ph=polyfit(x,y,N);
        %force low-freq val to be 0
        pv_ph(end)=0;
    
    end
%2. Load in San's corrections, just for comparison.  The only bit that is
%strictly needed here is the creation of the myFCTD.f frequency vector.
    good_ind = ~isnan(myFCTD.pressure);
    npts = sum(good_ind);
    df = FCTD_SalCorr.f_Ny/floor(npts/2);
    % creating the frequency axis
    myFCTD.f = (0:npts-1)'*df;
    myFCTD.f(myFCTD.f>FCTD_SalCorr.f_Ny) = myFCTD.f(myFCTD.f>FCTD_SalCorr.f_Ny)-2*FCTD_SalCorr.f_Ny;
    FCTD_SalCorr.GainFit = polyval(FCTD_SalCorr.GainPFit,myFCTD.f);
    FCTD_SalCorr.GainFit = FCTD_SalCorr.GainFit/FCTD_SalCorr.GainFit(1); % need to normalize Gain
    FCTD_SalCorr.PhsFit = polyval(FCTD_SalCorr.PhsPFit,myFCTD.f);
%3. Make the MHAcorr structure.

    %MHA low pass filter - tanh with half power point at fc.
    MHAcorr.sharpness=sharpness;
    MHAcorr.fc=fc;
    lp3=1/2-1/2*tanh(sharpness*(f - fc));

    ifpos=find(myFCTD.f>0);
    
    if use_phys==1
        %USE PHYSICAL MODEL
        %amplitude
        MHAcorr.G=interp1(-f,abs(H),myFCTD.f); %neg part
        MHAcorr.G(ifpos)=interp1(f,abs(H),myFCTD.f(ifpos)); %pos part
    
        %Same for phase (but make minus at -f)
        MHAcorr.ph=interp1(-f,-atan2(imag(H),real(H)),myFCTD.f); %neg part
        MHAcorr.ph(ifpos)=interp1(f,atan2(imag(H),real(H)),myFCTD.f(ifpos)); %pos part
        MHAcorr.method='physical';
        MHAcorr.L=L;
        MHA.tau=tau;
    else
        %USE POLYFIT
        %amplitude
        MHAcorr.G=interp1(-f,polyval(pv_gain,f),myFCTD.f); %neg part
        MHAcorr.G(ifpos)=interp1(f,polyval(pv_gain,f),myFCTD.f(ifpos)); %pos part
    
        %Same for phase (but make minus at -f)
        MHAcorr.ph=interp1(-f,-polyval(pv_ph,f),myFCTD.f); %neg part
        MHAcorr.ph(ifpos)=interp1(f,polyval(pv_ph,f),myFCTD.f(ifpos)); %pos part
        MHAcorr.method='polyfit';
        MHAcorr.pv_gain=pv_gain;
        MHAcorr.pv_ph=pv_ph;
    end
        
    %amplitude of new MHA LP filter
    MHAcorr.LP=interp1(-f,lp3,myFCTD.f); %neg part
    MHAcorr.LP(ifpos)=interp1(f,lp3,myFCTD.f(ifpos)); %pos part
    
    %Fill in zero freq which got missed in the interpolation
    
    MHAcorr.G=fftshift(NANinterp(fftshift(MHAcorr.G)));
    MHAcorr.ph=fftshift(NANinterp(fftshift(MHAcorr.ph)));
    MHAcorr.LP=fftshift(NANinterp(fftshift(MHAcorr.LP)));

% 4. Apply the spectral corrections.

        % FFT
        myFCTD.T = fft(myFCTD.temperature(good_ind),npts);
        myFCTD.C = fft(myFCTD.conductivity(good_ind),npts);
        myFCTD.P = fft(myFCTD.pressure(good_ind),npts);
        
        % correct the conductivity
        myFCTD.CCorr = myFCTD.C.*FCTD_SalCorr.GainFit.*exp(-1i*FCTD_SalCorr.PhsFit);
        % Low Pass filter
        myFCTD.TCorr = myFCTD.T.*FCTD_SalCorr.LPfilter(myFCTD.f);
        myFCTD.CCorr = myFCTD.CCorr.*FCTD_SalCorr.LPfilter(myFCTD.f);
        myFCTD.PCorr = myFCTD.P.*FCTD_SalCorr.LPfilter(myFCTD.f);
        
        % get back to physical units
        myFCTD.tCorr = real(ifft(myFCTD.TCorr));
        myFCTD.cCorr = real(ifft(myFCTD.CCorr));
        myFCTD.pCorr = real(ifft(myFCTD.PCorr));

        %----------MHA versions-------------
        %1. Create low-pass version of t,c,p

        %myFCTD.temperature, cond, pressure are now indexed into the entire FCTD record by 
        %ind.  good_ind is the non-nan part of this.  I'd like to grab an extra few points on either end of this 
        %to avoid edge effects in filtering.

        %ALB 09/02/2024 When running real time n_extra is bothering the
        %last bit of the "running profile"
        % I think it is fine to set it to 0 when
        % ind(end)+n_extra>length(FCTD timeseries)
        % It will fix itself when start a new cast. 
        % 
        if (ind(end)+n_extra)>length(FCTD.temperature)
            n_extra=0;
        end
        t_tmp=FCTD.temperature((ind(1)-n_extra):(ind(end)+n_extra));
        c_tmp=FCTD.conductivity((ind(1)-n_extra):(ind(end)+n_extra));
        p_tmp=FCTD.pressure((ind(1)-n_extra):(ind(end)+n_extra));

        %Filter.  We perform the response matching on only the high-passed signal to avoid problems with diff and reintegration. 
        %Then filter the extended records
        [b,a]=MHAButter(zs,zc); 
        tlow_tmp=filtfilt(b,a,t_tmp);
        clow_tmp=filtfilt(b,a,c_tmp);
        plow_tmp=filtfilt(b,a,p_tmp);

        %and truncate them back to the correct size.
        myFCTD.tlow=tlow_tmp(n_extra + (1:length(ind)));
        myFCTD.clow=clow_tmp(n_extra + (1:length(ind)));
        myFCTD.plow=plow_tmp(n_extra + (1:length(ind)));

        %Form the high-passed records on which we will perform the response
        %correction.  
        myFCTD.thigh=myFCTD.temperature(good_ind) - myFCTD.tlow;
        myFCTD.chigh=myFCTD.conductivity(good_ind) - myFCTD.clow;
        myFCTD.phigh=myFCTD.pressure(good_ind) - myFCTD.plow;

        % FFT
        myFCTD.T = fft(myFCTD.thigh,npts);
        myFCTD.C = fft(myFCTD.chigh,npts);
        myFCTD.P = fft(myFCTD.phigh,npts);

        % correct the conductivity
        myFCTD.CCorrMHA = myFCTD.C.*MHAcorr.G.*exp(-1i*MHAcorr.ph);
        % Low Pass filter all.
        myFCTD.TCorrMHA = myFCTD.T.*MHAcorr.LP;
        myFCTD.CCorrMHA = myFCTD.CCorrMHA.*MHAcorr.LP;
        myFCTD.PCorrMHA = myFCTD.P.*MHAcorr.LP;
        
        % get back to physical units, add back on low-pass signals.
        myFCTD.tCorrMHA = real(ifft(myFCTD.TCorrMHA)) + myFCTD.tlow;
        myFCTD.cCorrMHA = real(ifft(myFCTD.CCorrMHA)) + myFCTD.clow;
        myFCTD.pCorrMHA = real(ifft(myFCTD.PCorrMHA)) + myFCTD.plow;

%5. Finally, apply the thermal mass correction.
%This is from the MP processing script proc_CTD_MP52.m
% Thermal Mass algorithm is from SeaBird SeaSoft-Win32 manual (see
% Module12_AdvancedDataProcessing.pdf in dehysteresis folder)

%Sample frequency - do not change.
CTpar.freq = 16; % sample rate (Hz) for SBE49

T=myFCTD.tCorrMHA;
C=myFCTD.cCorrMHA;

%This code computes ctm, the conductivity thermal mass error.
    % compute/initialize temp diffs, cond corrections
    dTp = T;
    dTp(2:end) = diff(T);
    dTp(1) = dTp(2);
    dcdt = 0.1 * (1 + 0.006*(T-20)); %This is the expression for dcdt from SBE manual!
    ctm = 0*dTp;
    % a,b
    aa = 2 * CTpar.alfa / (2 + CTpar.beta/CTpar.freq);
    bb = 1 - (2*aa/CTpar.alfa);
    % compute corrections
    for i=2:length(C)
        ctm(i) = -1.0*bb*ctm(i-1) + aa*dcdt(i)*dTp(i);
    end
%    CC = C + ctm;


        %Put the corrected MHA fields into the expected fields of myFCTD
        %for passing to the next stage of processing (gridding etc).
        myFCTD.pressure = myFCTD.pCorrMHA;
        myFCTD.temperature = myFCTD.tCorrMHA;
        myFCTD.conductivity = myFCTD.cCorrMHA + do_thermal_mass_correction.*ctm; %ADD THERMAL MASS ERROR IF do_thermal_mass_correction is set to 1
        myFCTD.depth = sw_dpth(myFCTD.pressure,20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%end new code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

end
        num = num+1;
        DataGrid.time(num) = nanmean(FCTD.time(ind));
        for j=1:length(vars2Grid);
            DataGrid.(vars2Grid{j})(:,num) = bindata1d(DataGrid.depth,...
                myFCTD.depth, ...
                myFCTD.(vars2Grid{j}));
        end;
    end;
end;

DataGrid.depth = midpoints(DataGrid.depth);

if ~isfield(DataGrid,'temperature') || ~isfield(DataGrid,'pressure') || ~isfield(DataGrid,'conductivity')
    return;
end

DataGrid.salinity = sw_salt(DataGrid.conductivity*10/sw_c3515,DataGrid.temperature,DataGrid.pressure);
DataGrid.density = sw_pden(DataGrid.salinity,DataGrid.temperature,DataGrid.pressure,0);

% gridding in time
mintime = nanmin(DataGrid.time);
maxtime = nanmax(DataGrid.time);

DataGrid.tGrid.time = mintime:2*nanmedian(diff(DataGrid.time)):maxtime; % every minute
DataGrid.tGrid.depth = DataGrid.depth;


% allocate space to grid data
for i = 1:length(vars2Grid)
    DataGrid.tGrid.(vars2Grid{i}) = NaN(length(DataGrid.tGrid.depth),length(DataGrid.tGrid.time)-1);
end

for i = 1:length(DataGrid.tGrid.depth)
    for j = 1:length(vars2Grid)
        DataGrid.tGrid.(vars2Grid{j})(i,:) = bindata1d(DataGrid.tGrid.time,...
            DataGrid.time, DataGrid.(vars2Grid{j})(i,:));
    end
end

DataGrid.tGrid.salinity = sw_salt(DataGrid.tGrid.conductivity*10/sw_c3515,DataGrid.tGrid.temperature,DataGrid.tGrid.pressure);
DataGrid.tGrid.density = sw_pden(DataGrid.tGrid.salinity,DataGrid.tGrid.temperature,DataGrid.tGrid.pressure,0);

DataGrid.tGrid.time = midpoints(DataGrid.tGrid.time);

return;
end

function FCTD = FastCTD_FindCasts(FCTD,varargin)
% FCTD = FastCTD_FindCasts(FCTD);
%   finds when the FastCTD is going up and when it is going down
%   if 'UPCAST' is specified then the profile will search for upcast but by
%   defalt it would search for downcasts
%
%   In this data set we assume the data is in descending order of time
%
% Written by Jody Klymak
% Updated 2011 07 14 by San Nguyen
% Updated 2012 09 29 by San Nguyen



if ~isfield(FCTD,'pressure')
    return;
end;

% use 20 point median filter to smooth out the pressure field
p = medfilt1(FCTD.pressure,256);
% try to smooth out the data a bit
dp = conv2(diff(conv2(p,ones(256,1)/256,'same'),1,1)',ones(1,256)/256,'same');

%downLim = 0.1;
downLim = 0.025;
downCast = true;

persistent argsNameToCheck;
if isempty(argsNameToCheck);
    argsNameToCheck = {'downLim','threshold','upcast','downcast'};
end

index = 1;
n_items = nargin-1;
while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        error('MATLAB:FastCTD_FindCasts:wrongOption','Incorrect option specified: %s',varargin{index});
    end
    
    switch i
        case 1 % downLim
            if n_items == 1
                error('MATLAB:FastCTD_FindCasts:missingArgs','Missing input arguments');
            end
            downLim = varargin{index+1};
            if downLim == 0
                error('MATLAB:FastCTD_FindCasts:downLim0Err','The threshold cannot be zero!');
            end
            index = index +2;
            n_items = n_items-2;
        case 2 % downLim
            if n_items == 1
                error('MATLAB:FastCTD_FindCasts:missingArgs','Missing input arguments');
            end
            downLim = varargin{index+1};
            if downLim == 0
                error('MATLAB:FastCTD_FindCasts:downLim0Err','The threshold cannot be zero!');
            end
            index = index +2;
            n_items = n_items-2;
        case 3 % upcast
            downCast = false;
            
            index = index + 1;
            n_items = n_items - 1;
        case 4 % downcast
            downCast = true;
            
            index = index + 1;
            n_items = n_items - 1;
    end
end

%defining the threshold for going up and down
if downCast && downLim < 0
    downLim = -downLim;
elseif (~downCast) && downLim > 0% going up
    downLim = -downLim;
end

if downCast
    dn = find(dp>downLim);
else
    dn = find(dp<downLim);
end

% find all indices of going down

if isempty(dn)
    return;
end;

dn = [0, dn];


% find jumps in indices to indicate a start of a profile
startdown = dn(find(diff(dn)>1)+1);

if isempty(startdown);
    return;
end;

dn = dn(2:end);
FCTD.drop = 0*FCTD.time;

if dn(1)<startdown(1)
    startdown=[dn(1) startdown];
end;

if startdown(end)<dn(end);
    startdown = [startdown dn(end)];
end;


for i=1:(length(startdown)-1);
    in = intersect(startdown(i):startdown(i+1)-1,dn);
    FCTD.drop(in) = i;
end;
end
