function [FCTD] = add_microconductivity(FCTD)
    % Check for chi data
    if isfield(FCTD,'uConductivity') && ~isfield(FCTD,'chi')
        
        % Load chi_param
        chi_param=FCTD_DefaultChiParam;
        chi_param.fs=320;
        chi_param.min_spd=0.01; %TFO RR2410 upcast have a slow last part of up cast.
        
        % Reshape the data into a long vector
        FCTD.ucon=reshape(FCTD.uConductivity',20*length(FCTD.time),1)/2^24;
        
        % Make the higher resolution time vector for microconductivity
        FCTD.microtime=linspace(FCTD.time(1),FCTD.time(end),chi_param.fs./chi_param.fslow*length(FCTD.time))';
        
        % Convert microconductivity to chi
        FCTD = add_chi_microMHA_v2(FCTD,chi_param);
        
    end
    
    % Remove descriptive fields
    if isfield(FCTD,'chi_param')
        FCTD = rmfield(FCTD,'chi_param');
    end
    if isfield(FCTD,'chi_meta')
        FCTD = rmfield(FCTD,'chi_meta');
    end
   