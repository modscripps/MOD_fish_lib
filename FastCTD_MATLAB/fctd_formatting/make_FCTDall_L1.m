function [FCTDall,FCTDold] = make_FCTDall_L1(FCTDall,fctd_mat_dir)
% 1. Divide into down and upcasts
% 2. Add temperature/conductivity response matching code here instead of
% where it was in FastCTD_GridData. We should save the corrected
% temperature and conductivity in the FCTDall structure.
% 3. Add microconductivity chi calculation at the FCTDall level.

apply_response_matching_code = 0;

fprintf('Creating FCTDall_L1 from FCTDall_L0 and saving to %s\n',fctd_mat_dir);

%% Define downcasts and upcasts
% Make space for downcasts and upcasts
FCTDall.down = nan(size(FCTDall.time,1),size(FCTDall.time,2));
FCTDall.up = nan(size(FCTDall.time,1),size(FCTDall.time,2));
FCTDall.cast = nan(size(FCTDall.time,1),size(FCTDall.time,2));

% Find downcasts
FCTDall = FastCTD_FindCasts(FCTDall);
FCTDall.down(FCTDall.drop>0) = FCTDall.drop(FCTDall.drop>0);

% Find upcasts
FCTDall = FastCTD_FindCasts(FCTDall,'upcast');
FCTDall.up(FCTDall.drop>0) = FCTDall.drop(FCTDall.drop>0);

% Make a 'cast' field wheren downcasts are negative and upcasts are
% positive
FCTDall.cast(FCTDall.down>0) = -FCTDall.down(FCTDall.down>0);
FCTDall.cast(FCTDall.up>0) = FCTDall.up(FCTDall.up>0);

% Remove 'drop' since it is EITHER downcasts or upcasts depending on which
% way you called FastCTD_FindCasts last, and we now have the more specific
% fields 'down', 'up', and 'cast'.
FCTDall = rmfield(FCTDall,'drop');

% Create a copy of FCTDall so that you don't have to write over the
% original data. You can come back to FCTDold if you need to compare/debug.
FCTDold = FCTDall;

%% load correction factors for FCTD
FCTD_SalCorr = load('FCTD_SalinityCorrectionFactors_toCond.mat');
FCTD_SalCorr.GainPFit = FCTD_SalCorr.GainPFit_Dn;
FCTD_SalCorr.PhsPFit = FCTD_SalCorr.PhsPFit_Dn;

%% Define variables in FCTDall to carry forward
vars2grid_list = get_FCTD_fields(FCTDall);

% Loop through the casts and correct T and C
casts_u = unique(FCTDall.cast);
casts_n = casts_u(~isnan(casts_u));
[~,idx] = sort(abs(casts_n));
casts = casts_n(idx);

num = 0;
for iCast=1:length(casts)
    fprintf('\nProcessing cast %4.0f of %4.0f\n',iCast,length(casts));
    ind = find(FCTDall.cast==casts(iCast));

    myFCTD.time = FCTDall.time(ind);
    for j=1:length(vars2grid_list)
        %ALB add try/catch because it break in the fluor during TFO seamount 2023
        try
            myFCTD.(vars2grid_list{j}) = FCTDall.(vars2grid_list{j})(ind,:);
        catch
            myFCTD.(vars2grid_list{j})=nan.*ind;
        end
        % ALB end of hack
    end

    % If chi_param exists, add it
    if isfield(FCTDall,'chi_param')
        myFCTD.chi_param = FCTDall.chi_param;
    end

    if ~apply_response_matching_code
        fprintf('  Flag set to NOT apply response matching code.\n');
    elseif apply_response_matching_code
        if max(FCTDall.pressure(ind))-min(FCTDall.pressure(ind))>10 

            %         myFCTD.temperature = FCTD.temperature(ind);
            %         myFCTD.pressure = FCTD.pressure(ind);
            %         myFCTD.conductivity = FCTD.conductivity(ind);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %BEGIN MHA EDITS OF THIS FUNCTION 8/26/2024
            use_old_code=0; %Use 1 here to use San's legacy response matching code. It must provide the corrected variables myFCTD.pressure, temperature, conductivity and depth.
            if use_old_code
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Begin legacy code
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp(' Applying legacy response matching code')
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

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %end legacy code
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %begin new MHA code
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp(' Applying new response matching code')
                %0: Set parameters
                %MHA low pass filter.  Default (aggressive) values are sharpness = 2, fc=4.
                sharpness=.3; %defines how quickly the filter goes from 1 to 0 around fc.
                fc=.2; %cutoff freq.  The filter is 0.5 at fc.

                use_phys=0; %1 for physical model; 0 for polyfit
                zs=1/16; %sample interval in s - fixed at 1/16 sec.
                zc=10; %cutoff interval in s
                n_extra=16*zc*2; %get a factor of two more than the filter length on either side.

                %1 to compute spectra and polyval parameters on the fly; 0 to skip
                %and use defaults.
                compute_spectra=1;

                %Set to zero to ignore thermal mass correction; 1 to perform it.
                do_thermal_mass_correction = 1;
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

                    [Pt,fe] = pwelch(diffs(tdata),nfft,[],nfft,samplefreq,'psd'); %Set NFFT equal to the window length.  Verified that kdum is equal to k computed above.
                    Pt=Pt.';
                    Pt=Pt(2:length(Pt)); %delete f=0; normalization is already correct.

                    [Pc,~] = pwelch(diffs(cdata),nfft,[],nfft,samplefreq,'psd'); %Set NFFT equal to the window length.  Verified that kdum is equal to k computed above.
                    Pc=Pc.';
                    Pc=Pc(2:length(Pc)); %delete f=0; normalization is already correct.

                    %Now compute cross spectrum with same dofs etc.
                    [Pct,~] = cpsd(diffs(cdata),diffs(tdata),nfft,[],nfft,samplefreq,'onesided'); %Set NFFT equal to the window length.  Verified that kdum is equal to k computed above.
                    Pct=Pct.';
                    Pct=Pct(2:length(Pct)); %delete f=0; normalization is already correct.
                    %Maybe this is the best way to compute alpha...
                    if1=find(f < 1);% 1 Hz
                    alpha=sqrt(mean(Pt(if1),'omitmissing') ./ mean(Pc(if1),'omitmissing'));

                    %Next I'll try a polyfit
                    x=f;
                    y=abs(Pct)./Pc/alpha;
                    %ALB
                    % y=abs(Pct).^2./(Pc.*Pt);

                    pv_gain=polyfit(x,y,N);
                    %force low-freq val to be 1
                    pv_gain(end)=1;
                    % test=polyval(pv_gain,f);
                    % Phase
                    x=f;
                    y=atan2(imag(Pct),real(Pct));

                    pv_ph=polyfit(x,y,N);
                    %force low-freq val to be 0
                    pv_ph(end)=0;

                end %end if compute_spectra

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
                end %end if use_phys

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

                if (ind(end)+n_extra)>length(FCTDall.temperature)
                    n_extra=0;
                end
                t_tmp=FCTDall.temperature((ind(1)-n_extra):(ind(end)+n_extra));
                c_tmp=FCTDall.conductivity((ind(1)-n_extra):(ind(end)+n_extra));
                p_tmp=FCTDall.pressure((ind(1)-n_extra):(ind(end)+n_extra));

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
                for jj=2:length(C)
                    ctm(jj) = -1.0*bb*ctm(jj-1) + aa*dcdt(jj)*dTp(jj);
                end
                %    CC = C + ctm;


                %Put the corrected MHA fields into the expected fields of myFCTD
                %for passing to the next stage of processing (gridding etc).
                myFCTD.pressure = myFCTD.pCorrMHA;
                myFCTD.temperature = myFCTD.tCorrMHA;
                myFCTD.conductivity = myFCTD.cCorrMHA + do_thermal_mass_correction.*ctm; %ADD THERMAL MASS ERROR IF do_thermal_mass_correction is set to 1

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %end new code
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end %end if use old_code

        else
            fprintf('  Cast not long enough to apply thermal mass correction.\n')
        end %end if there are at least 2000 samples

    end %end if apply_response_matching_code

    % Replace data in FCTDall
    for iField=1:length(vars2grid_list)
            FCTDall.(vars2grid_list{iField})(ind,:) = myFCTD.(vars2grid_list{iField});
    end

    % Add OUT for debugging
    %FCTDall.OUT = myFCTD.OUT;

    clear myFCTD
end %end loop through casts


%% Process microconductivity data to make chi fields

if isfield(FCTDall,'uConductivity')
    FCTDall.chi = nan(size(FCTDall.time,1),size(FCTDall.time,2));
    FCTDall.eps_chi = nan(size(FCTDall.time,1),size(FCTDall.time,2));
    FCTDall.chi_tot = nan(size(FCTDall.time,1),size(FCTDall.time,2));

    % Check for gaps in data > 1 second. If there are gaps longer than 1
    % second, divide the data and process it in pieces to avoid errors
    % caused by interpolating over large periods.
    idx = find(diff(FCTDall.time)>days(seconds(1)));
    if isempty(idx)
        ranges = [1 length(FCTDall.pressure)];
    else
        if isscalar(idx)
            ranges = [1          idx;...
                      idx+1,     length(FCTDall.pressure)];
        elseif numel(idx)>1
            ranges = [1               idx(1);...
                      idx(1:end-1)+1, idx(2:end);...
                      idx(end)+1,     length(FCTDall.pressure)];
        end
    end

    for iR=1:length(ranges)
        inRange = ranges(iR,1):ranges(iR,2);

        fprintf('\nAdding chi to FCTDall range %3.0f of %3.0f.\n',iR,length(ranges));
        myFCTD.time = FCTDall.time(inRange);
        for j=1:length(vars2grid_list)
            try
                myFCTD.(vars2grid_list{j}) = FCTDall.(vars2grid_list{j})(inRange,:);
            catch
                myFCTD.(vars2grid_list{j})=nan.*inRange;
            end

        end

        if length(myFCTD.pressure)>2000 % ALB in the case we are using MHA new code we need at least 2000 samples, which is about 2 minutes (2000/16Hz). Since we now do this at the level of FCTDall and not in individual .mat files, this should always be okay.

            % Calculate chi if uConductivity exists
            if ~isfield(myFCTD,'chi_param')
                myFCTD.chi_param=FCTD_DefaultChiParam;
                myFCTD.chi_param.min_spd=0.01; %TFO RR2410 upcast have a slow last part of up cast.
                myFCTD.chi_param.plotit = 1;
                myFCTD.chi_param.fs = 320;
            end

            % Convert microconductivity to chi
            myFCTD = add_chi_microMHA_v3(myFCTD,myFCTD.chi_param);

        else
            fprintf('  Not enough data in range to process microconductivity.\n')
        end
    
        % Replace data in FCTDall
        for iField=1:length(vars2grid_list)
            FCTDall.chi(inRange,:) = myFCTD.chi;
            FCTDall.eps_chi(inRange,:) = myFCTD.eps_chi;
            FCTDall.chi_tot(inRange,:) = myFCTD.chi_tot;
        end
    end %end loop through range
end %end if there is microconductivity data

% Output data and save concatenated file
save(fullfile(fctd_mat_dir,'FCTDall_L1'),'FCTDall','-v7.3')