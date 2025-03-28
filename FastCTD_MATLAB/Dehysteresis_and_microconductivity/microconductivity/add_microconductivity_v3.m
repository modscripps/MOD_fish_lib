% MHA made add_microconductivity_v3 in Feb 2025
function [FCTD] = add_microconductivity_v3(FCTD)
%v3: MHA Feburary 2025. Skipped v2 to keep things consistent with the other
%routines that it calls.

% I do not see any use for lines 15 or 18 of the v2 code, as ucon and
% microtime are set in add_chi_microMHA_v3. So they are removed in this
% version.

% Note in v3 new fields are added to FCTD: w, chi, chi_tot and eps_chi.
% So the gridding code will need to be modified to grid these quantities.

% Load chi_param
chi_param=FCTD_DefaultChiParam;
chi_param.mfile = 'add_chi_microMHA_v3.m';
chi_param.fs=320;
chi_param.min_spd=0.01; %TFO RR2410 upcast have a slow last part of up cast.
chi_param.plotit = 0;

% Reshape the data into a long vector
%FCTD.ucon=reshape(FCTD.uConductivity',20*length(FCTD.time),1)/2^24;

% Make the higher resolution time vector for microconductivity
%FCTD.microtime=linspace(FCTD.time(1),FCTD.time(end),chi_param.fs./chi_param.fslow*length(FCTD.time))';

% Convert microconductivity to chi
FCTD = add_chi_microMHA_v3(FCTD,chi_param);
