% function [Pt_volt_f,Pt_Tg_k,chi,f,k,fc,kc,flag_tg_fc]=mod_efe_scan_chi(scan,fpo7_channel,Meta_Data,h_freq,FPO7noise)
function scan = mod_efe_scan_chi(scan,fpo7_channel,Meta_Data,h_freq,FPO7noise)

kmin=3;
% -------------------------------------------------------------------------
% Get constants and transfer functions
% NC - Changed w to absolute value so this also works for upcasts
w = abs(scan.w);

nfft=Meta_Data.PROCESS.nfft;
try
    Fs=Meta_Data.PROCESS.Fs_epsi;
catch
   Fs=Meta_Data.AFE.FS; 
end

switch fpo7_channel
    case 't1_volt'
        if isfield(Meta_Data,'AFE')
            %dTdV = Meta_Data.AFE.t1.cal;
            volts_to_C = Meta_Data.AFE.t1.volts_to_C; %[slope, intercept]
        else
            try
                %dTdV=Meta_Data.epsi.t1.dTdV;
            catch
                %dTdV=Meta_Data.epsi.t1.cal;
            end
        end
    case 't2_volt'
        if isfield(Meta_Data,'AFE')
            %dTdV = Meta_Data.AFE.t2.cal;
            volts_to_C = Meta_Data.AFE.t2.volts_to_C; %[slope, intercept]
        else
            try
            %dTdV=Meta_Data.epsi.t2.dTdV;
            catch
                %dTdV=Meta_Data.epsi.t2.cal;
            end
        end
    otherwise
        disp('wrong channel to compute chi, must be t1 or t2')
end

% If FPO7 noise is not specified, get it from Meta_Data
if nargin<5
    h_freq=get_filters_MADRE(Meta_Data,f);
    % get FPO7 channel average noise to compute chi
    switch Meta_Data.MAP.temperature
        case 'Tdiff'
            FPO7noise=load(fullfile(Meta_Data.paths.calibrations.fpo7,'FPO7_noise.mat'),'n0','n1','n2','n3');
        otherwise
            FPO7noise=load(fullfile(Meta_Data.paths.calibrations.fpo7,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
    end
end

h_freq.magsq=@(speed)(1./(1+((2*pi*(0.005*speed^(-0.32))).*f(:)).^2));
filter_TF=h_freq.FPO7(w);
dof=Meta_Data.PROCESS.dof;
% ---------------------------------------------------------------------
% Calculate spectra and chi

% Compute the frequency spectrum of timeseries in volts
[Pt_volt_f,f] = pwelch(detrend(scan.(fpo7_channel)),nfft,[],nfft,Fs,'psd');
k = f./w;

% Compute the frequency spectrum of timeseries in degrees C
[Pt_T_f,f] = pwelch(detrend(scan.(fpo7_channel)*volts_to_C(1) + volts_to_C(2)),nfft,[],nfft,Fs,'psd');

% Old way, using spectral dTdV - NC 10/06/25
% % Convert frequency spectrum of volt timeseries to frequency spectrum of
% % temperature in C
%Pt_T_f = (Pt_volt_f*(dTdV^2)) ./ filter_TF;

% Convert temperature frequency spectrum to temperature gradient wavenumber spectrum
Pt_Tg_k = ((2*pi*k).^2).*Pt_T_f.*w; %NC 9/2/21 - frequency spectrum should be MULTIPLIED by w, not divided

% Calculate chi
dk = mean(diff(k),'omitmissing');
fc_index = FPO7_cutoff(f,Pt_volt_f,FPO7noise);
fc = f(find(f<=f(fc_index),1,'last'));
kc = fc/w;
krange = find(k>=kmin & k<=k(fc_index));
chi = 6*scan.ktemp*dk.*sum(Pt_Tg_k(krange),'omitmissing');

Pxx=Pt_Tg_k(k>=kmin & k<=kc);
flag_tg_fc = fc_index<round(.95*length(f));

fitBatchelor = @(localepsi,localchi,local_nu,local_ktemp,k) ...
                (batchelor(localepsi,localchi,local_nu,local_ktemp,k)); %cpm 

chiSearch = chi.*[1e-1 1e1];
% get the correct FULL wavenumber range
SpecObs.k=k;
% get the cleaned spectrum corresponding to seg and the shear channel you want to look at.
SpecObs.P=Pt_Tg_k;
SpecObs.k=k;
SpecObs.kli=[kmin kc];
SpecObs.epsilon_final=scan.epsilon_final;
SpecObs.nu=scan.kvis;
SpecObs.ktemp=scan.ktemp;
SpecObs.dof=dof;

[chi_mle,~]=get_chi_mle(SpecObs,fitBatchelor,dof,chiSearch);


if (isnan(chi_mle) || isnan(SpecObs.epsilon_final))
    fom     = nan;
    fom_mle = nan;
else

    % compute fom
    %fom=compute_fom(scan,chi,k_cutoff,Pxx)
    fom     = compute_fom(SpecObs,chi,[kmin kc],Pxx);
    fom_mle = compute_fom(SpecObs,chi_mle,[kmin kc],Pxx);


    do_fom_fig = 0;
    if do_fom_fig
        figure(1)
        color_chi='g';
        if fom>1.15
            color_chi='r';
        end

        a=2;
        ax(a)=subplot(212);
        loglog(ax(a),k,Pt_Tg_k,'k-.')
        hold(ax(a),'on')
        loglog(ax(a),kin,Pxx,color_chi,'LineWidth',2)
        loglog(ax(a),k_model,P_model,'--','Color',[.5 .5 .5],'LineWidth',1)
        loglog(ax(a),kin,interp_P_model,'Color',[.5 .5 .5],'LineWidth',2)
        loglog(ax(a),k_model1,P_model1,'k--','LineWidth',1)
        loglog(ax(a),kin,interp_P_model1,'k','LineWidth',2)
        scatter(ax(a),kin(1),Pxx(1),20,'pr','filled')
        scatter(ax(a),kin(end),Pxx(end),20,'pr','filled')
        text(ax(a),2,1e-5,...
            ['\' sprintf('chi=%1.2e [˚C^2 s^{-1}]\n FOM= %1.2f ',...
            chi,fom)])

        hold(ax(a),'off')
        grid(ax(a),'on')
        ax(a).XLim=[1 1e3];
        ax(a).YLim=[1e-6 1e-1];
        pause(.1)

    end

end


% % Put new variables in the structure
% varList = {'Pt_volt_f','Pt_Tg_k','chi','fc','kc','flag_tg_fc'};
% for iVar=1:numel(varList)
%     scan.(varList{iVar}).(currChannel(1:2)) = eval(varList{iVar});
% end

scan.Pt_volt_f.(fpo7_channel(1:2))  = Pt_volt_f;
scan.Pt_Tg_k.(fpo7_channel(1:2))    = Pt_Tg_k;
scan.chi.(fpo7_channel(1:2))        = chi;
scan.chi_mle.(fpo7_channel(1:2))    = chi_mle;
scan.fc.(fpo7_channel(1:2))         = fc;
scan.kc.(fpo7_channel(1:2))         = kc;
scan.flag_tg_fc.(fpo7_channel(1:2)) = flag_tg_fc;
scan.kmin.(fpo7_channel(1:2))       = kmin;
scan.fom.(fpo7_channel(1:2))        = fom;
scan.fom_mle.(fpo7_channel(1:2))    = fom_mle;
end



function Psg=batchelor(epsilon,chi,nu,D,k)
% Usage: [k,Psg]=batchelor(epsilon,chi,nu,D,q);
%  inputs:
%    epsilon: turbulent dissipation rate, W/kg
%    chi: scalar dissipation rate, K^2/s or (c.u.)^2/s
%    nu: kinematic viscosity, m^2/s
%    D: scalar diffusivity, m^2/s
%    q: strain parameter, optional, q=3.7 is default
%  outputs:
%    k: wavenumber, cpm
%    Psg: power spectrum of scalar gradient, e.g. (K/m)^2/cpm if the
%     scalar is temperature and (c.u./m)^2/cpm if scalar is
%     salinity.
% Function: To evaluate the one-dimensional power spectrum of
%  scalar gradient using the theoretical form of Batchelor (1959).

%------------------------------------------------------------------------------
%
%        BATCHELOR              08-19-92               Ren-Chieh Lien
%
%        Batchelor temperature gradient spectrum
%
%        function batchelor(epsilon,chi,nu,D,q);
%
%        reference : 
%               Oakey, N. S., "Determination of the rate of dissipation of
%               turbulent energy from simultaneous temperature and velocity 
%               shear microstructure measurements", j.p.o., 12, 256-271, 1981.
%
%------------------------------------------------------------------------------

%alb change erf to erfc to change the slope of the roll-off
q = 3.7;
kb = (epsilon/nu/D^2)^(1/4);
a = sqrt(2*q)*2*pi*k/kb;
uppera = erfc(a/sqrt(2))*sqrt(pi/2);
g = 2*pi*a.*(exp(-a.^2/2) - a.*uppera);
Psg = g*sqrt(q/2)*(chi/kb/D);
Psg(Psg<=0)=0;
end

function [chi,misfit]=get_chi_mle(SpecObs,fitmyModel,dof,TDRSearch)
% find the healthy wavenumber range KMOIN and KMAX
ind=find(SpecObs.k>=SpecObs.kli(1) & SpecObs.k<=SpecObs.kli(2) & SpecObs.k>0);
% select only the Data you want to fit, i.e. on the healthy wavenumber
% range
SpecObs.k=SpecObs.k(ind);
SpecObs.P=SpecObs.P(ind);

%
% another anomymous function to fit nasmyth on the healthy wavenumber range
% feed tmpModel in mle_any_model
% TDR turbulence dissipation
%fitmyModel
% fitBatchelor = @(localepsi,localchi,k) ...
%                (batchelor(localepsi,localchi,scan.nu,scan.ktemp,k)); %cpm 

% tmpModel=@(TDR)(fitmyModel(SpecObs.epsilon_final,TDR,SpecObs.nu,SpecObs.ktemp,SpecObs.k));
tmpModel=@(TDR)(fitmyModel(SpecObs.epsilon_final,TDR,SpecObs.nu,SpecObs.ktemp,SpecObs.k));
% the magic happens here.
% SpecObs are the data, DOF is the degree of freedom, TDRSearch is the
% epsilon range to find a "solution", tmp model is the model we are trying
% to fit.
[chi,misfit.err]=mle_any_model(SpecObs,dof,TDRSearch, tmpModel, 0);

if isnan(chi)
    misfit.var=NaN;
    misfit.MAD=NaN;
else
    Pt=fitmyModel(SpecObs.epsilon_final,chi,SpecObs.nu,SpecObs.ktemp,SpecObs.k);
    [misfit.var, misfit.MAD]=ruddick_misfit(SpecObs.P,Pt);
end
if isempty(misfit.MAD)
    misfit.var=nan;
    misfit.MAD=nan;
end

end
function fom=compute_fom(SpecObs,chi,k_cutoff,Pxx)

    kmin=k_cutoff(1);
    kc=k_cutoff(2);
    P_model=batchelor(SpecObs.epsilon_final,chi,SpecObs.nu,SpecObs.ktemp,SpecObs.k);
    kin=SpecObs.k(SpecObs.k>=kmin & SpecObs.k<=kc);
    interp_P_model=interp1(SpecObs.k(~isnan(P_model)), ...
        P_model(~isnan(P_model)),...
        kin);
    adjust_chi=median(Pxx,'omitmissing')./ ...
          median(interp_P_model,'omitmissing');
    chi =chi.*adjust_chi;
    %Psg=batchelor(epsilon,chi,nu,D,k)
    P_model1=batchelor(SpecObs.epsilon_final,chi,SpecObs.nu,SpecObs.ktemp,SpecObs.k);
    interp_P_model1=interp1(SpecObs.k(~isnan(P_model1)), ...
        P_model1(~isnan(P_model1)),...
        kin);


    % high signal flag: The cut off frequency is very high.
    % this could mean that the whole scan is corrupt since the spectrum is way
    % above the noise floor.

    Pxx=SpecObs.P(SpecObs.k>=kmin & SpecObs.k<=kc);


    % ALB: Lueck 2024 says natural logarithm so I use log and not log10
    fom=log(Pxx(:)./interp_P_model1(:));
    mad_spec=mad(fom,1);
    sig_lnS=sqrt(5/4*(SpecObs.dof)^(-7/9));
    % Ns is the number of fourier coef in the kmin, kc range
    Ns=length(kin);
    Tm=0.8+sqrt(1.56/Ns);

    fom=mad_spec./sig_lnS./Tm;


end
