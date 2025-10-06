function Meta_Data=mod_epsi_linear_calibration_FP07(Profile,saveData)
%function Meta_Data=mod_epsi_temperature_spectra_v3(Meta_Data,Profile,Profile,title1,np,dsp,tscan)
%
% ALB Nov 17th 2024 
% changing to linear (Volt/Celcius) calibration of the FP07 
%
% ALB Feb 2019.
%
% input
% Meta_Data: All environment informations
% EPSI and CTD data form 1 profile defined by EPSI_create_profiles
% tscan = length in seconds of the segment used to compute the spectra
%
% ctd_df= sampling frequancy of the CTD instrument
% ctd_df=6; % RBR ctd_df=16; % SBE

% For the temperature correction, we want a good chunk of data. Try
% tscan = 50 or if the profile is too short, find the longest profile
% in the deployment and set tscan to about 0.5-0.8 times the profile
% length
disp('--- mod_epsi_linear_calibration_FP07.m ---')
Meta_Data=Profile.Meta_Data;
if ~isfield(Meta_Data.CTD,'cal')
    Meta_Data.CTD.cal.SN=000; 
end

% Options for saving and plotting data
if nargin<2
    saveData = 1;
end

   
T=filloutliers(Profile.ctd.T,'linear','movmean',100);
% check if t1 and t2 channels exist.
indt1=find(cellfun(@(x) strcmp(x,'t1'),Meta_Data.PROCESS.channels));
indt2=find(cellfun(@(x) strcmp(x,'t2'),Meta_Data.PROCESS.channels));

if ~isempty(indt1)
    Profile.epsi.t1_volt=fillmissing(Profile.epsi.t1_volt,'linear');
    Profile.epsi.t1_volt=filloutliers(Profile.epsi.t1_volt,'center','movmedian',1000);
    nanmask=~isnan(Profile.epsi.t1_volt);
    it1_volt=interp1(Profile.epsi.dnum(nanmask),Profile.epsi.t1_volt(nanmask),Profile.ctd.dnum);
    nanmask=~isnan(it1_volt);
    CalT1 = polyfit(it1_volt(nanmask),T(nanmask),1);
    %dTdV(1)=CalT1(2);
else
    CalT1=[nan nan];
end
if ~isempty(indt2)
    Profile.epsi.t2_volt=fillmissing(Profile.epsi.t2_volt,'linear');
    Profile.epsi.t2_volt=filloutliers(Profile.epsi.t2_volt,'center','movmedian',1000);
    nanmask=~isnan(Profile.epsi.t2_volt);
    it2_volt=interp1(Profile.epsi.dnum(nanmask),Profile.epsi.t2_volt(nanmask),Profile.ctd.dnum);
    nanmask=~isnan(it2_volt);
    CalT2 = polyfit(it2_volt(nanmask),T(nanmask),1);
    %dTdV(2)=CalT2(2);
else
    CalT2=[nan nan];
end


% Sensitivity of probe, nominal. Save dTdV in Meta Data it will be use in
% the batchprocess
if isfield(Meta_Data,'AFE')
    field_name = 'AFE';
elseif isfield(Meta_Data,'epsi')
    field_name = 'epsi';
end
%Meta_Data.(field_name).t1.cal=dTdV(1);
%Meta_Data.(field_name).t2.cal=dTdV(2);
Meta_Data.(field_name).t1.volts_to_C=CalT1; %[slope, intercept]
Meta_Data.(field_name).t2.volts_to_C=CalT2; %[slope, intercept]


% NC 4/4/25 - I need to sleep!
saveData = 0;
% Save dTdV in Meta_Data
if saveData
    save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');
    % TODO change this so we append the new dVdT.
    % Right now I am doing a quick an dirty open write save.
    %t1
    calibration_file =  ...
        sprintf('Calibration_%s.txt',Meta_Data.(field_name).t1.SN);
    calibration_path =  ...
        fullfile(Meta_Data.AFE.tempcal_path,Meta_Data.(field_name).t1.SN);

    if ~exist(calibration_path,'dir')
       mkdir(calibration_path);
    end
    fid=fopen(fullfile(calibration_path,calibration_file),"w");
    calibration_string=sprintf('%s\r\n %03s, %2.1f, %s\r\n',...
                        datestr(now),Meta_Data.(field_name).t1.SN, ...
                        dTdV(1), ...
                        Meta_Data.CTD.cal.SN);
    fwrite(fid,calibration_string);
    fclose(fid);
    % t2
    calibration_file =  ...
        sprintf('Calibration_%s.txt',Meta_Data.(field_name).t2.SN);
    calibration_path =  ...
        fullfile(Meta_Data.AFE.tempcal_path,Meta_Data.(field_name).t2.SN);

    if ~exist(calibration_path,'dir')
       mkdir(calibration_path);
    end
    fid=fopen(fullfile(calibration_path,calibration_file),"w");
    calibration_string=sprintf('%s\r\n %03s, %2.1f, %s\r\n',...
                        datestr(now),Meta_Data.(field_name).t2.SN, ...
                        dTdV(2), ...
                        Meta_Data.CTD.cal.SN);
    fwrite(fid,calibration_string);
    fclose(fid);
    
end


