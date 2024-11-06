function J = cost_T_serial(FCTD,TParams)
% This function calculates the cost function J for the temperature
% hysteresis correction. This version calculates the cost function over all
% the profiles of one serial number
%
% The function calculates dT/dz vs T for both the up and down directions.
% It then maps dT/dz from both directions onto the same temperature grid.
% Then it computes the mean squared error between dT/dz up and dT/dz down.
% The cost is the error for the particular profile
%
% Inputs: 
%   FCTD - the standard FCTD structure of profiles, without any
%       corrections
%   TParams - a parameter by which the difference between True 
%       and measured temperature is weighted
%   lagg - a parameter of the offset between the true and measured
%   temperature
%
% Ouptut:
%   J - a numerical value representing the goodness of the hysteresis
%   correction. The closer it is to zero, the better the effect of the
%   correction. 
%
%
%
% Alex Andriatis
% 2021-02-07

% Initialize - the FCTD structure must be a standard up down profile
% structure with at least variables time, pressure, and temperature

dirnames = {'up','down'}; 
nprof = length(FCTD.down.temperature);
FCTDdehyst=[];

% 1. Calculate modified up and down profiles using the parameters and
% calculate the temperature gradients

for n=1:2
    direction = dirnames{n};
    data = FCTD.(direction);
    for i=1:nprof
        P = data.pressure{i};
        T = data.temperature{i};
        
        That = deHyst_T(T,TParams);
        
        dz=diff(P);
        % Zero pressure difference messes up the calculation
        I=find(dz==0);
        dz(I)=[];
        dT=diff(That);
        dT(I)=[];
        dTdz = dT./dz;
        That(I+1)=[];
        switch n
            case 1
                Tmid = That(2:end); % Going up, use temp at top of step
            case 2
                Tmid = That(1:end-1); % Going down, use temp at top of step
        end
        FCTDdehyst.(direction).T{i} = Tmid;
        FCTDdehyst.(direction).dTdz{i} = dTdz;
    end
end

% 2. For each set of profiles, get the range of temperatures, choose the
% middle 5/6ths, and interpolate onto the finer grid. Then compute the
% standard deviation

Js = NaN(1,nprof);
for i=1:nprof
    Tup = FCTDdehyst.up.T{i};
    Tdown = FCTDdehyst.down.T{i};
    
    dTdzup = FCTDdehyst.up.dTdz{i};
    dTdzdown = FCTDdehyst.down.dTdz{i};
    
    tmax = min(max(Tup),max(Tdown));
    tmin = max(min(Tup),min(Tdown));
    
    if tmax<=tmin
        Js(i)=NaN;
        disp(['Bad set of profiles at index ' num2str(i)]);
        continue;
    end
    
    trange = tmax-tmin;
    tlim=[tmin+trange/12 tmax-trange/12]; % Use the central 5/6ths of the available temperature range.
    dTlim=5; % Sometimes weird stuff happens - want to eliminate spurious temperature steps over 5 degrees / meter
    
    I = find(Tup>tlim(1) & Tup<tlim(2) & abs(dTdzup)<dTlim); 
    Tup = Tup(I);
    dTdzup = dTdzup(I);
    [Tup,I] = unique(Tup);
    dTdzup = dTdzup(I);
    
    I = find(Tdown>tlim(1) & Tdown<tlim(2) & abs(dTdzdown)<dTlim);
    Tdown = Tdown(I);
    dTdzdown = dTdzdown(I);
    [Tdown,I]=unique(Tdown);
    dTdzdown = dTdzdown(I);
    
    if length(Tdown)<2||length(Tup)<2
        Js(i)=NaN;
        disp(['Bad set of parameters for profiles at index ' num2str(i)]);
        continue;
    end
    
    %tgrid = linspace(tlim(1),tlim(2),max(length(Tdown),length(Tup))*10);
    
    % Make a non-uniform grid
    tgrid = [Tup; Tdown];
    tgrid = unique(tgrid);
    
    dTdzup = interp1(Tup,dTdzup,tgrid);
    dTdzdown = interp1(Tdown,dTdzdown,tgrid);
    mse = mean((dTdzup-dTdzdown).^2,'omitnan');
    Js(i)=mse;
    
    toplot=0; % Diagnostic switch
    if toplot
    % Test plots
        figure(1);
        clf(1);
        hold on;
        plot(tgrid,dTdzup);
        plot(tgrid,dTdzdown);
        legend('Up','Down');
        xlabel('Temperature');
        ylabel('dTdz');
        title(sprintf('Profile %i of %i',i,nprof));
        pause(0.3);
    end
        
end
J = mean(Js,'omitnan');
end