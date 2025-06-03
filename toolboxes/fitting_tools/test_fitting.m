% Test ftting algorithms using same wavenumbers as Rolf!

clear; close all;
vmppath='/Users/cynthiabluteau/Dropbox/2_Shared/ATOMIX_Datasets/ShearProbes/Data/'; %Lueck';%Ice_VMP250upriser_ArcticOcean';%Lueck';
fname='VMP250_TidalChannel_024';
mpath=fullfile(vmppath,'Lueck'); % Ice_VMP250upriser_ArcticOcean' could be Lueck
%%
VMP=load_netcdf_data(fullfile(mpath,[fname,'.nc']),{'L3_spectra','L4_dissipation'});

%%
L3=VMP.L3_spectra;
L3.KVISC=visc35(L3.TEMP);
kmax=VMP.L4_dissipation.KMAX;

[nSeg,nSh]=size(VMP.L4_dissipation.EPSI);
ind=[];
for jj=1:nSh
    tt=find(VMP.L4_dissipation.METHOD(:,jj)==1);
    ind=union(tt,ind);
end

%%
% Init
L4.EPSI=NaN([nSeg nSh]); L4.MAD=L4.EPSI;
L4.EPSI_NL=L4.EPSI; L4.MAD_NL=L4.EPSI;
L4.EPSI_POWER=L4.EPSI;L4.MAD_POWER=L4.EPSI;
%%
for ii=1:length(ind)
    seg=ind(ii);
    for jj=1:nSh
        [L4.EPSI(seg,jj),misfit,NLin,Power]=refit_L3_spectra(L3,seg,jj,[0 kmax(seg)],[1e-6 1e-3]);
        L4.MAD(seg,jj)=misfit.MAD;
        L4.EPSI_NL(seg,jj)=NLin.epsi;
        L4.MAD_NL(seg,jj)=NLin.MAD;
        L4.EPSI_POWER(seg,jj)=Power.epsi;
        L4.MAD_POWER(seg,jj)=Power.MAD;
    end
end

%%
for ii=1:2
    subplot(2,2,ii)
    if ii==1
        plot(VMP.L4_dissipation.MAD,L4.MAD,'.'); hold on
        plot(VMP.L4_dissipation.MAD,L4.MAD_NL,'+'); 
        %plot(VMP.L4_dissipation.MAD,L4.MAD_POWER,'d'); 
        plot([.1 .5],[.1 .5],'k-')
        title('MAD')
    else
        loglog(VMP.L4_dissipation.EPSI,L4.EPSI,'.'); hold on
        loglog(VMP.L4_dissipation.EPSI,L4.EPSI_NL,'+'); 
        loglog(VMP.L4_dissipation.EPSI,L4.EPSI_POWER,'d'); 
        ax=logscatterlines(gca,[1e-5 1e-4]);
        title('\epsilon [W kg^{-1}]')
    end
    
    xlabel('PI')
    ylabel('Tester')
    axis square
end

subplot(2,2,3:4)
plot(L4.MAD,VMP.L4_dissipation.EPSI./L4.EPSI,'.'); hold on
plot(L4.MAD,VMP.L4_dissipation.EPSI./L4.EPSI_NL,'+')
xlabel('MAD tester')
ylabel('\epsilon_{PI}/\epsilon_{tester}')