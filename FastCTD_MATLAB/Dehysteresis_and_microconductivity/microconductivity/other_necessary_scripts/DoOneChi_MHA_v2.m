function out=DoOneChi_MHA_v2(FCTD,tlim,chi_param)
%function out=DoOneChi_MHA_v2(FCTD,tlim,chi_param)
%Compute chi for the specified time period in tlim using the parameters in chi_param.  If the latter is not specified, defaults will be used.
%The FCTD structure must have had micro data processed by
%add_chi_microMHA_v2.m.
%
%v2 is by MHA August 2024. Following AP00 appendix, we now:
%
%Compute conductivity gradient spectra
%Compute dCdT from SEB formula (function of temperature)
%Correct spectrum for the antialias filter
%Integrate observed corrected spectrum from 1 cpm to kmax=fn/speed =
%160/speed.
%
%This is a highly underestimated spectrum, by a factor that depends on
%epsilon.  To correct for this one needs Tz and N2 and to assume a constant
%gamma following AP00.  I will do that using the gridded data.

if nargin < 3
    chi_param=FCTD_DefaultChiParam;
end

%Extract params from params struct
fs=chi_param.fs;
fslow=chi_param.fslow;
nsec=chi_param.nsec;
D=chi_param.D;
kmax=chi_param.kmax;
kmin=chi_param.kmin;
min_spd=chi_param.min_spd;
dt_sec=chi_param.dt_sec;
dTdC=chi_param.dTdC;
gain=chi_param.gain;
offset=chi_param.offset;
plotit=chi_param.plotit;

%Get indices
i1=find(FCTD.microtime>tlim(1)& FCTD.microtime<tlim(2));
i2=find(FCTD.time>tlim(1)&FCTD.time<tlim(2));



%Get indexed data.
data=FCTD.ucon_corr(i1); %microconductivity data in S/m, sampled at 320Hz
datac=FCTD.conductivity(i2); %This is SBE C in S/m, sampled at 16 Hz
datat=FCTD.temperature(i2); %This is SBE T in deg C

%MHA change 8/2024: let's compute dCdT from the SBE formula
dCdT_SBE=0.1*(1+0.006*(nanmean(datat) - 20));
dTdC=1./ dCdT_SBE;

%Mean fall rate
w=nanmean(FCTD.dPdt(i2));
spd=abs(w); %speed

%mean P
p=nanmean(FCTD.pressure(i2));

%set new kmax lims to avoid noise at high wavenumber
fn=fs/2; %Nyquist
fmax=fn/2; %integrate to half the Nyquist based on observed noise spectra.  
% This also avoids the need to correct for the difference versus derivative operator.

kmin=1;
kmax=fmax / spd;

%Let's do everything from here on out on wavenumber spectra of spatial quantities.

%High-freq quantities
WINDOW=fs*nsec;
ks=fs/spd;
dt=1/fs;
dz=spd*dt;

%low-freq quantities

kslow=fslow/spd;
dtlow=1/fslow;
dzlow=spd*dtlow;
WINDOW_low=fslow*nsec;


%escape if the data is not long enough
if length(data)<WINDOW
    out.chi_stupid=NaN;
    out.chi1=NaN;
    out.w=NaN;
    out.pres=NaN;

else

    %Compute chi from the standard formula as 6*D*(RMS dTdz)^2
    chi1=nan;
    chi_stupid=6*D*dTdC^2*nanstd(diff(data)*fs/w).^2;

    %Compute the wavenumber spectrum of dCdz from microconductivity.
    [P,k] = pwelch(diff(NANinterp(data))/dz,WINDOW,[],WINDOW,ks,'psd'); %Set NFFT equal to the window length

    %Antialias filter of microconductivity electronics in SBE7 per AP00
    sinc2=sinc(pi*k./max(k)).^2;

    %Corrected spectrum
    pcorr=P./sinc2.';

    if spd > min_spd
        ichi=find(k > kmin & k < kmax);
        dk=k(2)-k(1); %freq res
        chi1=6*D*dTdC^2*nansum(pcorr(ichi))*dk;
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

end
