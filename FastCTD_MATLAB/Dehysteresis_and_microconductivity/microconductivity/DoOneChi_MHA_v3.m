function out=DoOneChi_MHA_v3(FCTD,tlim,chi_param)
    %function out=DoOneChi_MHA_v3(FCTD,tlim,chi_param)
    %Compute chi for the specified time period in tlim using the parameters in chi_param.  If the latter is not specified, defaults will be used.
    %The FCTD structure must have had micro data processed by
    %add_chi_microMHA_v2.m.
    %
    %v2 is by MHA August 2024. Following AP00 appendix, we now:
    %
    %Compute conductivity gradient spectra
    %Compute dCdT from SBE formula (a function of temperature)
    %Correct spectrum for the antialias filter
    %Integrate observed corrected spectrum from 1 cpm to kmax=fn/speed =
    %50/speed.
    %
    %This is a highly underestimated estimate of chi, by a factor that depends on
    %epsilon.  To correct for this one needs Rrho and N2 and to assume a constant
    %gamma following AP00.  
    
    %v3: February 2025.
    %Change to integration to 12.5 cpm or 50 Hz / w cpm.
    %
    %output: 
    %    out.chi_stupid=chi_stupid; %Disregard.  This is simply the variance
    %    computed in the time domain.
    %    out.chi1=chi1; %The integral of the temperature gradient spectrum from
    %                   k=1 cpm up to kmax.  
    %    out.w=w;       %vertical velocity used for the estimate.
    %    out.pres=p;    %Mean pressure of the estimate.
    %    out.k=k;       %wavenumber vector, cpm
    %    out.P=pcorr;   %The wavenumber spectrum of conductivity gradient,
    %                   \Phi_dCdz [S/m^2 / cpm].  To obtain temperature
    %                   gradient wavenumber spectrum \Phi_TG(k), multiply by
    %                   dTdC^2.
    %    out.dTdC=dTdC; %dTdC in K / (S/m), computed from SBE formula.  
    %    out.kmax=kmax; %output max wavenumber used in cpm.  
    
    if nargin < 3
        chi_param=FCTD_DefaultChiParam;
        chi_param.mfile = 'add_chi_microMHA_v3.m';
    end
    
    %Extract params from params struct
    fs=chi_param.fs;
    fslow=chi_param.fslow;
    nsec=chi_param.nsec;
    D=chi_param.D;
    min_spd=chi_param.min_spd;
    
    % kmax=chi_param.kmax;
    % kmin=chi_param.kmin;
    % dt_sec=chi_param.dt_sec;
    % dTdC=chi_param.dTdC;
    % gain=chi_param.gain;
    % offset=chi_param.offset;
    % plotit=chi_param.plotit;
    
    %Get indices
    i1=find(FCTD.microtime>tlim(1)& FCTD.microtime<tlim(2));
    i2=find(FCTD.time>tlim(1)&FCTD.time<tlim(2));
    
    % NC added 2/19/25 - can't do anything if indices are empty
    if isempty(i1) || isempty(i2)
        out.chi_stupid=NaN;
        out.chi1=NaN;
        out.w=NaN;
        out.pres=NaN;
    else 
    
    %Get indexed data.
    data=FCTD.ucon_corr(i1); %microconductivity data in S/m, sampled at 320Hz
    datac=FCTD.conductivity(i2); %This is SBE C in S/m, sampled at 16 Hz
    datat=FCTD.temperature(i2); %This is SBE T in deg C
    
    %MHA change 8/2024: let's compute dCdT from the SBE formula
    dCdT_SBE=0.1*(1+0.006*(mean(datat,"omitmissing") - 20));
    dTdC=1./ dCdT_SBE;
    
    %Mean fall rate
    w=mean(FCTD.dPdt(i2),'omitmissing');
    spd=abs(w); %speed
    
    %mean P
    p=mean(FCTD.pressure(i2),'omitmissing');
    
    %set new kmax lims to avoid noise at high wavenumber
    fn=fs/2; %Nyquist
    % fmax=fn/2; %integrate to half the Nyquist based on observed noise spectra.  
    % ALB the sinc2 has a first dip at 50Hz Should we force the integration to stop
    % at 50Hz
    % fmax=fn/2; %integrate to half the Nyquist based on observed noise spectra.  
    fmax=50; %integrate to half the Nyquist based on observed noise spectra.  
    % This also avoids the need to correct for the difference versus derivative operator.
    
    %February 2025: use 12.5 cpm
    kmin=1;
    kmax1=fmax / spd;
    kmax=min(12.5,kmax1); %Use 12.5 cpm, but do not in any case integrate to greater than 50 Hz.
    
    %Let's do everything from here on out on wavenumber spectra of spatial quantities.
    
    %High-freq quantities
    WINDOW=fs*nsec;
    ks=fs/spd;
    dt=1/fs;
    dz=spd*dt;
    
    %low-freq quantities
    
    % dtlow=1/fslow;
    % kslow=fslow/spd;
    % dzlow=spd*dtlow;
    % WINDOW_low=fslow*nsec;
    
    
    %escape if the data is not long enough or if data is all nan
    if length(data)<WINDOW || all(isnan(data))
        out.chi_stupid=NaN;
        out.chi1=NaN;
        out.w=NaN;
        out.pres=NaN;
    
    else
    
        %Compute chi from the standard formula as 6*D*(RMS dTdz)^2
        chi1=nan;
        chi_stupid=6*D*dTdC^2*std(diff(data)*fs/w,[],'omitmissing').^2;
    
        %Compute the wavenumber spectrum of dCdz from microconductivity.
        [P,k] = pwelch(detrend(diff(NANinterp(data))/dz),WINDOW,[],WINDOW,ks,'psd'); %Set NFFT equal to the window length
    
        %Antialias filter of microconductivity electronics in SBE7 per AP00
        sinc2=sinc_nopi(pi*k./max(k)).^2;
    
        %Corrected spectrum
        pcorr=P(:)./sinc2(:);
        
        % if FCTD.time>datenum('11/26/2024 05:08:22')
        %     loglog(k,sinc2)
        %     hold on
        %     loglog(k,P)
        %     loglog(k,pcorr)
        %     loglog([kmax kmax],[1e-10 1e-2],'k--')
        %     loglog([1/spd 1/spd],[1e-10 1e-2],'k--')
        %     hold off
        %     grid on
        %     ylim([1e-10 1e2])
        %     pause
        % end
        
    
        if spd > min_spd
            ichi=find(k > kmin & k < kmax);
            dk=k(2)-k(1); %freq res
            chi1=6*D*dTdC^2*sum(pcorr(ichi),'omitmissing')*dk;
        end
    
        %Compute the wavenumber spectrum of dCdz from sbe conductivity.
    %    if length(datat)>WINDOW_low && length(datat)>WINDOW_low
    %        [Pc,kc] = pwelch(diff(NANinterp(datac))/dzlow,WINDOW_low,[],WINDOW_low,kslow,'psd'); %Set NFFT equal to window length
    
            %Compute the wavenumber spectrum of dTdz from sbe temperature.
    %        [Pt,kt] = pwelch(diff(NANinterp(datat))/dzlow,WINDOW_low,[],WINDOW_low,kslow,'psd'); %Set NFFT equal to window length
    %    end
    
    
        out.chi_stupid=chi_stupid;
        out.chi1=chi1;
        out.w=w;
        out.pres=p;
        out.k=k;
        out.P=pcorr;
        out.dTdC=dTdC; %2/2025 change.  Output dTdC.
        out.kmax=kmax;
        out.idx_microtime = i1;
        out.idx_time = i2;
        out.fs = fs;
        out.ks = ks;
        out.WINDOW = WINDOW;
    end %end if window is too short
    end %end if indices are empty
    