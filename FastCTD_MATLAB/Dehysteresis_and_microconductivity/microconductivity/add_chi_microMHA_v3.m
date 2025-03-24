function FCTD=add_chi_microMHA_v3(FCTD,chi_param)
%function FCTD=add_chi_microMHA_v3(FCTD)
%Matthew Alford
%March 2019
%Version 0.2
%
%As of 15 February 2025, v3 calls DoOneChi_MHA_v3.m.
%
%Major changes since 0.1:
%1) chi1 and chi2 are now switched; ie chi1 is the integrated quantity and
%chi2 is stupid (time domain).
%2) Chi calc and plotting are functionized.
%3) Default chi_params structure is created and passed.
%
%Unpack microconductivity data, remove preemphasis filter, loop through in
%time windows of specified length and compute chi in those windows by integrating the wavenumber spectrum.
%
%In v3 I also follow the appendix of Alford and Pinkel 2000b (AP00b) and compute the density1 ratio
%Rrho and N2 from the data, and then estimate the total chi and epsilon_chi
%that is consistent with stratified turbulence.
%
%The essence of the method is that we measure the level of the temperature
%gradient spectrum out to a wavenumber kmax, which is *not* near the
%Batchelor wavenumber.  So it is under-resolved by a factor that increases
%the greater epsilon is, which is not known.  With knowledge of Rrho and
%N2 and by assuming gamma = 0.2, we can compute Jb from chi_pe, and then
%use that to determine epsilon which in turn determines how under-resolved
%the chi estimate was in the first place.  AP00b computed the expression
%for eps_chi in closed form, which we use here.  The constant c was
%recomputed in Feb 2025 as MHA found AP00b's value of 2e-5 was 2.5 times
%too small.
%
%Note this estimate is uncertain and seems to be somewhat high when N2 is low.  It should be more
%carefully evaluated using the gridded data and used with caution. In
%particular, a hard-wired Rrho = 5 is used at present which should be
%improved upon for the best estimates.  The dependence of eps_chi on Rrho
%is weak so this is probably a minor effect.
%
%See ExamineFCTD-chi-2025.ipynb and ExamineFCTD-chi-2025-2.ipynb where I
%diagnosed these issues in the previous version _2 of these codes and
%developed and tested the eps_chi calculations going into the v3 code.
%
%MHA Feb 2025
%
%Output: The following fields are added to the FCTD structure:
% chi (an underestimated chi equal to 6D\int_(1 cpm)^(kmax) where kmax is
%        12.5 cpm or 50 Hz / w, whichever is smaller (12.5 cpm unless w < 4
%        m/s). Units: K^2/s.
%chi_tot (an estimate of the total chi following AP00b.  Units: K^2/s.
%eps_chi (an estimate of epsilon from the observed chi following AP00b).
%   Units: W/kg
%chi_meta: metadata structure.
%chi_param: parameters used
%
%Dependencies:
%CenteredConv
%remove_sbe_preemphasisMHA.m
%DisplayProgress
%DoOneChi_MHA_v3
%
%

%%%%%%%%%%%%%%%%%%%%
%Set parameters
%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    chi_param=FCTD_DefaultChiParam;
    chi_param.mfile = 'add_chi_microMHA_v3.m';
end

chi_meta.units='degC^2/s';
%8/2024 change: change description to match actual v2 fields
chi_meta.description='chi1 is 6D times the integral of the temperature gradient wavenumber spec from kmin to kmax.';
chi_meta.mfiles='add_chi_microMHA_v3.m, DoOneChi_MHA_v3.m';
%Create a higher-freq time vector
% Check for gaps in data > 1 second before making a new time array
idx = find(diff(FCTD.time)>days(seconds(1)));
if ~isempty(idx)
    disp('Skipping this cast')
    disp('Gaps in data > 1 second found. Check data used for uconductivity chi calculation.')
    disp('add_chi_micro_MHA_v3 is outputting nans for this cast')

    FCTD.chi_param=chi_param;
    FCTD.chi_meta=chi_meta;
    FCTD.chi = nan.*FCTD.time;
    FCTD.eps_chi = nan.*FCTD.time;
    FCTD.chi_tot = nan.*FCTD.time;

else


    FCTD.microtime=linspace(FCTD.time(1),FCTD.time(end),chi_param.fs./chi_param.fslow*length(FCTD.time))';

    i1=1:length(FCTD.microtime);
    %Also compute fall rate in m/s
    FCTD.dPdt=CenteredConv(diffs(FCTD.pressure)*chi_param.fslow,1,4*chi_param.fslow); %smooth over a few sec

    %Now unpack micro and express as raw volts - convert from counts
    if chi_param.fs==160
        FCTD.ucon=reshape(FCTD.uConductivity.',10*length(FCTD.time),1)/2^16; %pre SOM/epsi - 160 Hz @ 16 bits
    elseif chi_param.fs==320
        FCTD.ucon=reshape(FCTD.uConductivity.',20*length(FCTD.time),1)/2^24; %post SOM/epsi - 320 Hz @ 24 bits
    end
    %reshape takes them columnwise so transpose

    data=FCTD.ucon(i1); %This is the uncalibrated voltage

    %
    plot(FCTD.time,FCTD.uConductivity,'b')
    hold on
    plot(FCTD.microtime,FCTD.ucon,'r');
    a = [];

    %This is the function that does all the work.
    datac=remove_sbe_preemphasisMHA(NANinterp(data),chi_param.fs)*chi_param.gain+chi_param.offset;

    %Store corrected output field with the scaling and offset applied
    FCTD.ucon_corr=nan*FCTD.ucon;
    FCTD.ucon_corr(i1)=datac;

    %Now loop through and compute spectra for each block

    %chi_param.plotit=0;
    tout=FCTD.time(1):(chi_param.dt_sec/24/3600)/8:FCTD.time(end);
    chi_all=nan*tout; %This is chi by simple time-domain RMS'ing, in bins
    chi2_all=nan*tout; %This is chi by integrating the wavenumber spectrum
    w_all=nan*tout;
    disp(sprintf('  Calculating chi for %7.0f scans',length(tout)))
    for c=1:length(tout)

        %DisplayProgress(c,10000)
        tlim=tout(c)+[-1 1].*chi_param.dt_sec/2 /24/3600; %Make a time window surrounding that center time

        %Use v3
        %    out=DoOneChi_MHA_v2(FCTD,tlim,chi_param);
        out=DoOneChi_MHA_v3(FCTD,tlim,chi_param);

        w_all(c)=out.w;
        chi2_all(c)=out.chi_stupid;
        chi_all(c)=out.chi1;

        % NC - debugging 3/5/25
        OUT{c} = out;

        % Print progress status
        if chi_param.plotit %2/2025 change: only output status if plotit is set to 1
            if mod(c,20)==0 && mod(c,100)~=0
                fprintf(num2str(c))
            elseif mod(c,100)==0
                fprintf([num2str(c) '\n'])
            elseif c==length(tout)
                fprintf('. \n')
            else
                fprintf('.')
            end
        end

    end

    % NC - debugging 3/5/25
    FCTD.OUT = OUT;

    %% Now output to the FCTD grid
    FCTD.chi=interp1(tout,chi_all,FCTD.time);
    %FCTD.chi2=interp1(tout,chi2_all,FCTD.time); %Feb 2025: let's not output chi_stupid
    %which is, well, stupid.

    FCTD.w=interp1(tout,w_all,FCTD.time);

    FCTD.chi_param=chi_param;
    FCTD.chi_meta=chi_meta;


    %% Now compute the fully resolved quantities.
    % Calculate salinity and density, to be deleted later. Sometimes FCTD may already have these. In case of that, call them salinity1 and density1 so those fields can safely be removed. This function is for dealing with microstructure. It is not the function in which we want to calculate salinity and density. That belongs somewhere else.
    FCTD.salinity1=sw_salt(FCTD.conductivity*10/sw_c3515,FCTD.temperature,FCTD.pressure); %
    FCTD.density1 = sw_pden(FCTD.salinity1,FCTD.temperature,FCTD.pressure,0)-1000;

    %Then compute Tz and N2
    nt=500;
    rhoo=1026;
    g=9.8;
    FCTD.n2=g/rhoo*diffs(CenteredConv(FCTD.density1,1,nt))./diffs(CenteredConv(FCTD.pressure,1,nt));
    FCTD.Tz=diffs(CenteredConv(FCTD.temperature,1,nt))./diffs(CenteredConv(FCTD.pressure,1,nt));
    FCTD.Sz=diffs(CenteredConv(FCTD.salinity1,1,nt))./diffs(CenteredConv(FCTD.pressure,1,nt));

    FCTD.alpha = sw_alpha(FCTD.salinity1,FCTD.temperature,FCTD.pressure);
    FCTD.beta = sw_beta(FCTD.salinity1,FCTD.temperature,FCTD.pressure);

    FCTD.Rrho=FCTD.alpha .* FCTD.Tz ./ FCTD.beta ./ FCTD.Sz;

    %Let's now be brave and compute eps_chi following Alford and Pinkel 2000b
    gamma=0.2;
    g=9.8;
    rhoo=1026;
    Rrho=5; %choose a constant value for now.  This is known to be incorrect, though the dependency on Rrho is weak
    % - can redo better in gridding.

    %alpha=1.2e-4; %This should properly be computed from in situ T,S

    c=5e-5; %Note AP00 used 2.5e-5; I find 5e-5 for kmax = 12.5 cpm
    p=-0.4; %power law value from AP00.

    FCTD.chi_pe_meas=g*FCTD.alpha/rhoo.*FCTD.chi./FCTD.n2.*(1+1./Rrho.^2);
    FCTD.eps_chi=(2*gamma*c).^(-1/(p+1)) .* FCTD.chi_pe_meas.^(1./(p+1));

    %Then we can compute chi now that we know epsilon

    %Fraction resolved
    FCTD.r=c.*FCTD.eps_chi.^p;

    %And chi
    FCTD.chi_tot = FCTD.chi ./ FCTD.r;

    %% Now do some house keeping - remove fields we don't want to output.
    FCTD=rmfield(FCTD,'Tz');
    FCTD=rmfield(FCTD,'Sz');
    FCTD=rmfield(FCTD,'Rrho');
    FCTD=rmfield(FCTD,'n2');
    FCTD=rmfield(FCTD,'alpha');
    FCTD=rmfield(FCTD,'beta');
    FCTD=rmfield(FCTD,'salinity1'); %Delete salinity calculated in this function
    FCTD=rmfield(FCTD,'density1'); %Delete density calculated in this function
    FCTD=rmfield(FCTD,'chi_pe_meas');
    FCTD=rmfield(FCTD,'r');


end



