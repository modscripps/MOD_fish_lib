function [varm, MAD2]=ruddick_misfit(Pw,Pt)
% Ruddick et al misfit criteria
% Required inputs:
%  Pw: Spectral Observations
%  Pt: The theoretical spectrum



MAD2=(1./length(Pt)).*nansum(abs((Pw./Pt)-nanmean(Pw./Pt)));
varm=var(Pw./Pt);