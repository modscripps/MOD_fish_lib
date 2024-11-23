function SBEcal=get_CalSBE_v2(str_dcal)

%  reads str and parse the SBEcal data 
    str_dcal_split=strsplit(str_dcal,'\n');

    i=1;
    SBEcal.SN=str_dcal_split{i}(end-3:end);
    % Temperature Cal
    i=2;
    SBEcal.TempCal_date=str_dcal_split{i}(end-8:end);
    i=3;
    SBEcal.ta0=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=4;
    SBEcal.ta1=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=5;
    SBEcal.ta2=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=6;
    SBEcal.ta3=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=7;
    SBEcal.toffset=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    % Conductivity Cal
    i=8;
    SBEcal.CondCal_date=str_dcal_split{i}(end-8:end);
    i=9;
    SBEcal.g=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=10;
    SBEcal.h=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=11;
    SBEcal.i=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=12;
    SBEcal.j=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=13;
    SBEcal.pcor=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=14;
    SBEcal.tcor=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=15;
    SBEcal.cslope=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    
    % Pressure Cal
    i=16;
    SBEcal.PresCal_date=str_dcal_split{i}(end-8:end);
    i=17;
    SBEcal.pa0=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=18;
    SBEcal.pa1=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=19;
    SBEcal.pa2=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=20;
    SBEcal.ptca0=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=21;
    SBEcal.ptca1=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=22;
    SBEcal.ptca2=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=23;
    SBEcal.ptcb0=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=24;
    SBEcal.ptcb1=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=25;
    SBEcal.ptcb2=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=26;
    SBEcal.ptempa0=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=27;
    SBEcal.ptempa1=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=28;
    SBEcal.ptempa2=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));
    i=29;
    SBEcal.poffset=str2double(str_dcal_split{i}(strfind(str_dcal_split{i},'=')+1:end));

end
