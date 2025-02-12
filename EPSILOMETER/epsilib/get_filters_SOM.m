function H=get_filters_SOM(Meta_Data,f)
%function H=get_filters_MADRE(Meta_Data,f)
%% Define the electronic filter, ADC filters.
%  The electronics filter Helec = charge amp filter (Hca) + Gain (Hg)
%  The ADC filter is a simple sinc^4 and the gain of the ADC is 2. Other options are available. 
%
%  TODO: enter argument to be able to change the ADC filter without hard
%  coding
% Feb 2019 ALB

forig=f;
if contains(Meta_Data.vehicle_name(:).','sa_apex')
    Fs_epsi=320;
    [~,f]  =  pwelch(0*(1:Meta_Data.PROCESS.nfft),...
    Meta_Data.PROCESS.nfft,[], ...
    Meta_Data.PROCESS.nfft, ...
    Fs_epsi,'psd');
end


try
try    
switch Meta_Data.AFE.s1.ADCfilter
    case 'sinc4'
        Hs1filter=(sinc(f/(2*f(end)))).^4;
end
catch
switch Meta_Data.AFE.f1.ADCfilter
    case 'sinc4'
        Hs1filter=(sinc(f/(2*f(end)))).^4;
end
end
switch Meta_Data.AFE.t1.ADCfilter
    case 'sinc4'
        Ht1filter=(sinc(f/(2*f(end)))).^4;
end
switch Meta_Data.AFE.a1.ADCfilter
    case 'sinc4'
        Ha1filter=(sinc(f/(2*f(end)))).^4;
end
catch
switch Meta_Data.Firmware.ADCfilter
    case 'sinc4'
        Hs1filter=(sinc(f/(2*f(end)))).^4;
        Ht1filter=(sinc(f/(2*f(end)))).^4;
        Ha1filter=(sinc(f/(2*f(end)))).^4;
end
end

if contains(Meta_Data.vehicle_name(:).','sa_apex')
        Hs1filter=interp1(f,Hs1filter ,forig);
        Ht1filter=interp1(f,Ht1filter ,forig);
        Ha1filter=interp1(f,Ha1filter ,forig);
end
f=forig;


% shear channels
%charge amp filter

%ca_filter = load('FILTER/charge_coeffilt.mat');
% calibrator gain 2.5 dB ~ Vout/Vin=1.33 for a shear probe of 750pF

%ca_filter = load('FILTER/charge_coeffilt_09312019.mat'); %from network analysis
ca_filter = load('cap1nFres200Meg_5KohmInput.mat'); %from network analysis
epsi_ca   = interp1(ca_filter.freq,ca_filter.coef_filt ,f);
%gain_ca   = -(max(20*log10(epsi_ca))+2.5); %TODO set a coef to get a -2.5dB TF as a first approx. I might get fancier by getting probe Cap
%gain_ca   = -(max(20*log10(epsi_ca))-2.5); %try with a 0 dB gain.
gain_ca   = 1; %try with a 0 dB gain.
H.electshear= epsi_ca*gain_ca;% 
%H.electshear= epsi_ca*10.^(gain_ca/20);% 
H.gainshear=1;
H.adcshear=H.gainshear.* Hs1filter;
H.shear=(H.electshear(:) .* H.adcshear(:)).^2;

%% FPO7 channels

H.gainFPO7=1;
H.electFPO7 = H.gainFPO7.*Ht1filter;
%speed% convert to m/s
%tau=0.005 * speed^(-0.32); % thermistor time constant
H.magsq=@(speed)(1 ./ (1+((2*pi*(0.005 * speed^(-0.32))).*f(:)).^2)); % magnitude-squared no units
H.phase=@(speed)(-2*atan( 2*pi*f(:)*(0.005 * speed^(-0.32))));   % no units
try
switch Meta_Data.AFE.temp_circuit
    case {'Tdiff','tdiff'}
        Tdiff_filter = load('FILTER/Tdiff_filt.mat');
        Tdiff_H = interp1(Tdiff_filter.freq,Tdiff_filter.coef_filt ,f);
        H.Tdiff=Tdiff_H;
        H.FPO7=@(speed)(H.electFPO7(:).^2 .* H.magsq(speed) .* H.Tdiff.^2);
    otherwise
        H.FPO7=@(speed)(H.electFPO7(:).^2 .* H.magsq(speed));
end
catch
    switch Meta_Data.MAP.temperature
        case {'Tdiff','tdiff'}
            Tdiff_filter = load('FILTER/Tdiff_filt.mat');
            Tdiff_H = interp1(Tdiff_filter.freq,Tdiff_filter.coef_filt ,f);
            H.Tdiff=Tdiff_H;
            H.FPO7=@(speed)(H.electFPO7(:).^2 .* H.magsq(speed) .* H.Tdiff(:).^2);
        otherwise
            H.FPO7=@(speed)(H.electFPO7(:).^2 .* H.magsq(speed));
    end
end



%% Accel channels
H.gainAccel  = 1;
H.electAccel = (H.gainAccel(:).*Ha1filter(:)).^2;

end