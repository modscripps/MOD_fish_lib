function [chi_refit,misfit]=mod_efe_scan_chi_mle(scan)

% add missing variance to chi

fitBatchelor = @(epsi,chi,k)(batchelor(epsi,chi,scan.nu,scan.ktemp,k)); %cpm 

% ---------------------------------------------------------------------
% Calculate spectra and chi

krange1 = find(scan.k>2 & scan.k<=scan.kc(1));
krange2 = find(scan.k>2 & scan.k<=scan.kc(2));

% Calculate observed chi between 2 and kc
dk = mean(diff(scan.k),'omitmissing');
chi1_obs = 6*scan.ktemp*dk.*sum(scan.Pt.t1(krange1),'omitmissing');
chi2_obs = 6*scan.ktemp*dk.*sum(scan.Pt.t2(krange2),'omitmissing');

% Calculate theoretical chi between 2 and kc 
%P_B = vmp_batchelor_spectra(kB,k,chi,kappa,q)
kB = (scan.epsilon/scan.nu/scan.ktemp^2)^(1/4);
Pbatch1=batchelor(scan.epsilon,...
                  scan.chi(1),...
                  scan.nu,...
                  scan.ktemp,...
                  scan.k);
Pbatch2=batchelor(scan.epsilon,...
                  scan.chi(2),...
                  scan.nu,...
                  scan.ktemp,...
                  scan.k);

chi1_batch = 6*scan.ktemp*dk.*sum(Pbatch1(krange1),'omitmissing');
chi2_batch = 6*scan.ktemp*dk.*sum(Pbatch2(krange2),'omitmissing');



% Calculate chi with missing variance 
Pbatch1_adjust=chi1_obs./chi1_batch .* Pbatch1;
Pbatch2_adjust=chi2_obs./chi2_batch .* Pbatch2;
chi_refit1 = 6*scan.ktemp*dk.*sum(Pbatch1_adjust,'omitmissing');
chi_refit2 = 6*scan.ktemp*dk.*sum(Pbatch2_adjust,'omitmissing');

check_adjust1=sum(Pbatch1_adjust(krange1),'omitmissing')-sum(scan.Pt.t1(krange1));
check_adjust2=sum(Pbatch2_adjust(krange2),'omitmissing')-sum(scan.Pt.t2(krange2));

if(abs(check_adjust1)>1e-15)
disp("adjust1 wrong")
end
if(abs(check_adjust2)>1e-15)
disp("adjust2 wrong")
end


if isnan(chi_refit1)
    misfit.t1.var=NaN; 
    misfit.t1.MAD=NaN;
else
    Pt=fitBatchelor(scan.epsilon,chi_refit1,scan.k); 
    [misfit.t1.var, misfit.t1.MAD]=ruddick_misfit(scan.Pt.t1(krange1),...
                                                             Pt(krange1));
end
if isempty(misfit.t1.MAD)
    misfit.t1.var=NaN; 
    misfit.t1.MAD=NaN;
end
if isnan(chi_refit2)
    misfit.t2.var=NaN; 
    misfit.t2.MAD=NaN;
else
    Pt=fitBatchelor(scan.epsilon,chi_refit2,scan.k); 
    [misfit.t2.var, misfit.t2.MAD]=ruddick_misfit(scan.Pt.t1(krange2),...
                                                  Pt(krange2));
end
if isempty(misfit.t2.MAD)
    misfit.t2.var=NaN; 
    misfit.t2.MAD=NaN;
end


chi_refit=[chi_refit1 chi_refit2];



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
Psg = sqrt(q/2)*(chi/kb/D)*g;
Psg(Psg<=0)=0;
end
