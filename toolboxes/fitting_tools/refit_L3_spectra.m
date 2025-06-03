function [epsi,misfit,NLin,Power]=refit_L3_spectra(L3,seg,ii_sh,kli,epsiSearch)
% Runs fit on individual spectra b/w user specified kli (min and max) cpm
%seg=13; % segment to process
%ii_sh: % shear probe to process
% kli: min and max kli for the seg in question
% Output:
%   NLin: nonlinear output with fields epsilon and MAD
if nargin<5
    epsiSearch=[1e-11 1e-2];
end


% anomymous fucnction to create a nasmyth spectrum with epsilon value and
% wavenumber range (KVISC comes from L3)
fitmyModel=@(epsi,k)(nasmyth(epsi, L3.KVISC(seg), k)); %cpm 

%%
% get DOF. why? 
if length(L3.DOF)>1
   L3.DOF=L3.DOF(seg); 
end
% get the correct FULL wavenumber range
SpecObs.k=L3.KCYC(seg,:)';
% get the cleaned spectrum corresponding to seg and the shear channel you want to look at. 
SpecObs.P=squeeze(L3.SH_SPEC_CLEAN(seg,:,ii_sh))';
% find the healthy wavenumber range KMOIN and KMAX
ind=find(SpecObs.k>=kli(1) & SpecObs.k<=kli(2) & SpecObs.k>0);
% select only the Data you want to fit, i.e. on the healthy wavenumber
% range
SpecObs.k=SpecObs.k(ind);
SpecObs.P=SpecObs.P(ind);

%%
tt=1; %Why tt=1?
% another anomymous function to fit nasmyth on the healthy wavenumber range
% feed tmpModel in mle_any_model
tmpModel=@(chiS)(fitmyModel(chiS,SpecObs.k));
% the magic happens here. 
% SpecObs are the data, DOF is the degree of freedom, epsiSearch is the
% epsilon range to find a "solution", tmp model is the model we are trying
% to fit.
[epsi,misfit.err]=mle_any_model(SpecObs,L3.DOF,epsiSearch, tmpModel, 0);
           
if isnan(epsi)
    misfit.var(tt)=NaN; 
    misfit.MAD(tt)=NaN;
else
    Pt=fitmyModel(epsi(tt),SpecObs.k); 
    [misfit.var(tt), misfit.MAD(tt)]=ruddick_misfit(SpecObs.P,Pt);
end
   
%% Try nonlinear least square of Nasmyth

[NLin.epsi,r,J,COVB]  = nlinfit(SpecObs.k,SpecObs.P,fitmyModel,epsi);
Pt=fitmyModel(NLin.epsi,SpecObs.k); 
[NLin.MAD]=ruddick_misfit(SpecObs.P,Pt);

%% Try now Inertial subrange usning Constant from Sreenivasan
% fitmyModel=@(epsi,k)(inertial_model(epsi,k));

%% Power fit. Assumes supplied k are in deed in the inertial subrange. Should check kmax<0.02*eta
[A, ~ ,Aci]=fitpowerlaw(SpecObs.k,SpecObs.P,1/3);
kolmC=((4./3)*18./55).*1.5; 
Power.epsi=(A./((2*pi).^(4/3).*kolmC)).^(3/2);
Pt=inertial_model(Power.epsi,SpecObs.k); 
[Power.MAD]=ruddick_misfit(SpecObs.P,Pt);
end


function Pk=inertial_model(epsi,k)
Ai=(4/3)*18/55;
C=1.5; % Kolmogorov sconstant
Pk=(2*pi).^(4/3).*Ai.*C.*epsi.^(2/3).*k.^(1/3);

end