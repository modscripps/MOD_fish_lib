function [epsilon,epsilon_non_adjust,kmin,kc,method]=eps1_mmp(k,Psheark,kvis,kmax)
% eps1_mmp
%   Usage: epsilon=eps1_mmp(k,Psheark,kvis,w);
%      k is a vector array with wavenumber in cpm
%      Psheark is a vector array with the shear spectrum
%      kvis is the kinematic viscocity, in m^2/s
%      w is the vehicle speed, in m/s
%      dk is the elementary waveSD! 1enumber determined by eps1_mmp
%      epsilon is the estimated dissipation rate, in W/kg
%      kc is the wavenumber at which integration is cutoff, in cpm
%   Function: To integrate airfoil shear versus k spectra to
%      obtain epsilon.  The algorithm is similar to that of
%      Wesson & Gregg, 1994, JGR, 99, -9877, but uses Panchev's
%      universal spectrum for reference and stops integration to 
%      avoid vibration peaks.  This limit, kmax, is determined by eps2_mmp.

%    Epsilon is determined iteratively, in three stages.
%    (1) integration between 2 and 10 cpm, eps1
%          The observed spectrum is interpolated onto a wavenumber grid

%       between 2 and 10 cpm, with 0.2 cpm steps.  The integral, shear 10,
%       is compared with an integral of Panchev's universal spectrum over
%       the same grid.  If log10(shear10)>-3, 2 to 10 cpm lies solely in
%       the inertial subrange, which does not depend on viscosity.  Epsilon
%       is obtained directly from a polynomial fit of epsilon to shear10 for
%       the Panchev spectrum.  If log10(shear10)<=-3, 2 to 10 cpm contains at
%       least some of the viscous rolloff and epsilon/7.5 nu is obtained from
%       a polynomial fit to minimize errors due to viscosity.
%
%    (2) integration to the wavenumber containing 90% variance of Panchev's
%       spectrum evaluated with eps1 and nu, eps2
%         The upper limit for integration is reduced if it exceeds kmax, the limit
%       determined for noise-free spectra by script eps2_mmp.m.  The lower bound
%       is at 1 cpmk.  If fewer than 4 spectral estimates are in the wavenumber band, 
%       no integration is done and epsilon is set to 1e-10, taken as the base level.  
%       The estimate is raised by script epsilon_correct.m if the signal has been 
%       reduced by probe attenuation.  
%
%    (3) repeat of (2) with wavenumber determined from eps2    

method=0; % we are not in the inertial subrange. 

% kmin=0;
% dKI=0.2;
% KI=(2:dKI:10); % wavenumber array for interpolation


% kmin =3 cpm for Apex
kmin=3;
% kmin =2 cpm for normal operations
%kmin=2;
% if k(1)==0
%     kmin=k(2);
% else
%     kmin=k(1);
% end


% dk=nanmean(diff(k));
% first estimation of epsilon of the sum of the shear variance is too high
% and falls into the inertial subrange (no roll off). I think it assumes that
% the shear variance directly follows a Panchev spectrum and predict
% epsilon directly without influence of the viscosity. I have no idea how
% these coefficients are computed.
% 
eps_fit_shear10=[8.6819e-04, -3.4473e-03, -1.3373e-03, 1.5248, -3.1607];

% If the shear variance include a part of the roll-off, we use these coef
% which (I think) give a panchev spectrum for a given shear variance value 
shtotal_fit_shear10=[6.9006e-04, -4.2461e-03, -7.0832e-04, 1.5275, 1.8564];


% first estimate, using Psheark between 2 & 10 cpm
% Interpolate Psheark onto 0.2 cpm grid & integrate
% Only this estimate is interpolated, as it is the only one input to
% a polynomial integrated with the same grid
krange=find(k>=kmin & k<10); 
% P_interpolated=interp1(k(krange),Psheark(krange),KI);
% ALB change to nansum since coherence correction can introduces nansclear 
% shear10=nansum(P_interpolated)*0.2;
if numel(krange)>2
shear10=trapz(k(krange),Psheark(krange));
%
% estimate epsilon using poly fits to log10(shear10)
logshear10=log10(shear10);
if logshear10>-3 % 2-10 cpm lies entirely in inertial subrange
	log10eps=polyval(eps_fit_shear10,logshear10);
	eps1=10^log10eps;
    method=1;% we ARE in the inertial subrange. 
else
	log10_sheartotal=polyval(shtotal_fit_shear10,logshear10);
	eps1=7.5*kvis*10^log10_sheartotal;
end

% second estimate: we use the first estimate of epsilon to find the
% kolmogrov scale and re-do the integration.
kc2 = 2*0.0816*( eps1  / kvis^3 )^(1/4);  % k for 90% variance of Panchev spectrum
if kc2>kmax
	kc2=kmax; % limit set by noise spectrum
end
krange=find(k>=kmin & k<=kc2);
% eps2=7.5*kvis*nansum(Psheark(krange))*dk/.9; % .9 we want to get 90% of the shear variance
if length(krange)>2
    eps2=7.5*kvis*trapz(k(krange),Psheark(krange)); %
else
    eps2=1e-12;
end

% third estimate: same as before.
kc=0.0816*( eps2 / kvis^3 )^(1/4);

if kc > kmax
	kc=kmax;
end
krange=find(k>=kmin & k<=kc);
if length(krange)>2
    eps3 = 7.5*kvis*trapz(k(krange),Psheark(krange)); 
else
    eps3=1e-12;
end

if eps3<1e-11 
	epsilon=1e-11;
else
    %ALB add a correction to fit best eps3 to the spectrum.
    % i think it shoulds only matter for very low epsilons
    [Ppan3,~]=nasmyth(eps3,kvis,k);
    idx_kmax2=find(k>=kc,1,'first');
    idx_kmax1=find(k>=kmin,1,'first');
    idx_range=idx_kmax1:idx_kmax2;
    % Pnasm_var=trapz(k(idx_range),Ppan3(idx_range));
    % obs_var=trapz(k(idx_range),Psheark(idx_range));
    adjust=median(Psheark(idx_range),'omitmissing')./ ...
           median(Ppan3(idx_range),'omitmissing');
    % ALB revisit this everytime you have doubt 
    % eps4=eps3*obs_var./Pnasm_var;
    % eps4=eps3*obs_var./Pnasm_var;
    eps4=eps3*adjust;
    % epsilon_option=[eps3 eps4 eps5];
    % [Ppan4,~]=nasmyth(eps4,kvis,k);
    % [Ppan5,~]=nasmyth(eps5,kvis,k);
    % best_fit=...
    %     [log10(Psheark(idx_range)./Ppan3(idx_range)) ...
    %      log10(Psheark(idx_range)./Ppan4(idx_range))  ...
    %      log10(Psheark(idx_range)./Ppan5(idx_range))];
    %  epsilon_best=epsilon_option(mad(best_fit)==min(mad(best_fit)));
    % 
    mf=epsilon2_correct(eps4, ...
                        kvis,kmin,kc);
	epsilon=mf*eps4;
	epsilon_non_adjust=mf*eps3;
    check_epsi=0;
    if check_epsi
        [Ppan6,~]=nasmyth(epsilon_best,kvis,k);
        [Ppan,~]=nasmyth(epsilon,kvis,k);
        figure(1)
        cla
        loglog(k,Psheark)
        hold on
        loglog(k,Ppan3,'r')
        scatter(k(find(k>kc,1,'first')),Psheark(find(k>kc,1,'first')),50,'pr','filled');
        grid on
        loglog(k,Ppan4,'k')
        loglog(k,Ppan5,'m')
        loglog(k,Ppan6,'c')
        loglog(k,Ppan,'k.-')
        legend('obs','eps3','eps4','eps5','epsi-good','epsi-final')

        pause(.5)
    end
end

%selecting the closest k for gc;
kc=k(find(k>kc,1,'first'));
else
    epsilon = nan;
    kmin    = nan;
    kc      = nan;
    method  = nan;
end

end







% function [epsilon,kc]=eps1_mmp(k,Psheark,kvis,kmax)
% eps1_mmp
%   Usage: epsilon=eps1_mmp(k,Psheark,kvis,w);
%      k is a vector array with wavenumber in cpm
%      Psheark is a vector array with the shear spectrum
%      kvis is the kinematic viscocity, in m^2/s
%      w is the vehicle speed, in m/s
%      dk is the elementary waveSD! 1enumber determined by eps1_mmp
%      epsilon is the estimated dissipation rate, in W/kg
%      kc is the wavenumber at which integration is cutoff, in cpm
%   Function: To integrate airfoil shear versus k spectra to
%      obtain epsilon.  The algorithm is similar to that of
%      Wesson & Gregg, 1994, JGR, 99, -9877, but uses Panchev's
%      universal spectrum for reference and stops integration to 
%      avoid vibration peaks.  This limit, kmax, is determined by eps2_mmp.

%    Epsilon is determined iteratively, in three stages.
%    (1) integration between 2 and 10 cpm, eps1
%          The observed spectrum is interpolated onto a wavenumber grid

%       between 2 and 10 cpm, with 0.2 cpm steps.  The integral, shear 10,
%       is compared with an integral of Panchev's universal spectrum over
%       the same grid.  If log10(shear10)>-3, 2 to 10 cpm lies solely in
%       the inertial subrange, which does not depend on viscosity.  Epsilon
%       is obtained directly from a polynomial fit of epsilon to shear10 for
%       the Panchev spectrum.  If log10(shear10)<=-3, 2 to 10 cpm contains at
%       least some of the viscous rolloff and epsilon/7.5 nu is obtained from
%       a polynomial fit to minimize errors due to viscosity.
%
%    (2) integration to the wavenumber containing 90% variance of Panchev's
%       spectrum evaluated with eps1 and nu, eps2
%         The upper limit for integration is reduced if it exceeds kmax, the limit
%       determined for noise-free spectra by script eps2_mmp.m.  The lower bound
%       is at 1 cpmk.  If fewer than 4 spectral estimates are in the wavenumber band, 
%       no integration is done and epsilon is set to 1e-10, taken as the base level.  
%       The estimate is raised by script epsilon_correct.m if the signal has been 
%       reduced by probe attenuation.  
%
%    (3) repeat of (2) with wavenumber determined from eps2    
% 
% dKI=0.2;
% KI=(2:dKI:10); % wavenumber array for interpolation
% 
% dk=nanmean(diff(k));
% % first estimation of epsilon of the sum of the shear variance is too high
% % and falls into the inertial subrange (no roll off). I think it assumes that
% % the shear variance directly follows a Panchev spectrum and predict
% % epsilon directly without influence of the viscosity. I have no idea how
% % these coefficients are computed.
% % 
% eps_fit_shear10=[8.6819e-04, -3.4473e-03, -1.3373e-03, 1.5248, -3.1607];
% 
% % If the shear variance include a part of the roll-off, we use these coef
% % which (I think) give a panchev spectrum for a given shear variance value 
% shtotal_fit_shear10=[6.9006e-04, -4.2461e-03, -7.0832e-04, 1.5275, 1.8564];
% 
% 
% % first estimate, using Psheark between 2 & 10 cpm
% % Interpolate Psheark onto 0.2 cpm grid & integrate
% % Only this estimate is interpolated, as it is the only one input to
% % a polynomial integrated with the same grid
% krange=find(k>=2 & k<10); 
% P_interpolated=interp1(k(krange),Psheark(krange),KI);
% % ALB change to nansum since coherence correction can introduces nansclear 
% % shear10=nansum(P_interpolated)*0.2;
% shear10=nansum(P_interpolated)*dKI;
% %
% % estimate epsilon using poly fits to log10(shear10)
% logshear10=log10(shear10);
% if logshear10>-3 % 2-10 cpm lies entirely in inertial subrange
% 	log10eps=polyval(eps_fit_shear10,logshear10);
% 	eps1=10^log10eps;
% else
% 	log10_sheartotal=polyval(shtotal_fit_shear10,logshear10);
% 	eps1=7.5*kvis*10^log10_sheartotal;
% end
% 
% % second estimate: we use the first estimate of epsilon to find the
% % kolmogrov scale and re-do the integration.
% kc2 = 0.0816*( eps1  / kvis^3 )^(1/4);  % k for 90% variance of Panchev spectrum
% if kc2>kmax
% 	kc2=kmax; % limit set by noise spectrum
% end
% krange=find(k>=2 & k<=kc2);
% eps2=7.5*kvis*nansum(Psheark(krange))*dk/.9; % .9 we want to get 90% of the shear variance
% 
% % third estimate: same as before.
% kc=0.0816*( eps2 / kvis^3 )^(1/4);
% if kc > kmax
% 	kc=kmax;
% end
% krange=find(k>=2 & k<=kc);
% eps3 = 7.5*kvis*nansum(Psheark(krange))*dk; 
% 
% if eps3<1e-11 || length(krange) < 4
% 	epsilon=1e-11;
% else
%   mf=epsilon2_correct(eps3,kvis);
% 	epsilon=mf*eps3;
% end
