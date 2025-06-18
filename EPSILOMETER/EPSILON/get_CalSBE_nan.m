function SBEcal=get_CalSBE_nan()

    SBEcal.SN=nan;
    % Temperature Cal
    SBEcal.TempCal_date=nan;
    SBEcal.ta0=nan;
    SBEcal.ta1=nan;
    SBEcal.ta2=nan;
    SBEcal.ta3=nan;
    % line=fgetl(fid);
    % SBEcal.toffset=str2double(line(strfind(line,'=')+1:end));
    
    % Conductivity Cal
    SBEcal.CondCal_date=nan;
    SBEcal.g=nan;
    SBEcal.h=nan;
    SBEcal.i=nan;
    SBEcal.j=nan;
    SBEcal.tcor=nan;
    SBEcal.pcor=nan;
    
    %line=fgetl(fid);
    %SBEcal.cslope=str2double(line(strfind(line,'=')+1:end));
    
    % Pressure Cal
    SBEcal.PresCal_date=nan;
    SBEcal.pa0=nan;
    SBEcal.pa1=nan;
    SBEcal.pa2=nan;
    SBEcal.ptca0=nan;
    SBEcal.ptca1=nan;
    SBEcal.ptca2=nan;
    SBEcal.ptcb0=nan;
    SBEcal.ptcb1=nan;
    SBEcal.ptcb2=nan;
    SBEcal.ptempa0=nan;
    SBEcal.ptempa1=nan;
    SBEcal.ptempa2=nan;


end
