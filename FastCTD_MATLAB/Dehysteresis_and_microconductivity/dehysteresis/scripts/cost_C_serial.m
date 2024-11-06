function J = cost_C_serial(FCTD,aC,bC,gamma,lag)
% This function calculates the cost function J for the conductivity
% correction. Uses parameters gamma, alpha, TimeC, and temperature
% parameters tau and lag
%
% The function calculates dC/dz vs C for both the up and down directions.
% It then maps dC/dz from both directions onto the same conductivity grid.
% Then it computes the mean squared error between dC/dz up and dC/dz down.
% The cost is the error for the particular profile
%
% Inputs: 
%   FCTD - the standard FCTD structure of profiles, without any
% corrections
%   x0 - the conductivity parameters
%   offset - an integer offset between the temperature and conductivity
%   cells
%   tparams - temperature parameters
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
nprof = length(FCTD.down.conductivity);
FCTDdehyst=[];

% 1. Calculate modified up and down profiles using the parameters and
% calculate the temperature gradients

for n=1:2
    direction = dirnames{n};
    data = FCTD.(direction);
    
    for i=1:nprof
        P = data.pressure{i};
        C = data.conductivity{i};
        That = data.temperature{i};
        
        Ccorrected = deHyst_C(That,C,aC,bC,gamma,lag);
        
        dz=diff(P);
        % Zero pressure difference messes up the calculation
        I=find(dz==0);
        dz(I)=[];
        dC=diff(Ccorrected);
        dC(I)=[];
        dCdz = dC./dz;
        Ccorrected(I+1)=[];
        switch n
            case 1
                Cmid = Ccorrected(2:end); % Going up, use temp at top of step
            case 2
                Cmid = Ccorrected(1:end-1); % Going down, use temp at top of step
        end

        FCTDdehyst.(direction).C{i} = Cmid;
        FCTDdehyst.(direction).dCdz{i} = dCdz;
    end
end

% 2. For each set of profiles, get the range of temperatures, choose the
% middle 5/6ths, and interpolate onto the finer grid. Then compute the
% standard deviation

Js = NaN(1,nprof);
for i=1:nprof
    Cup = FCTDdehyst.up.C{i};
    Cdown = FCTDdehyst.down.C{i};

    dCdzup = FCTDdehyst.up.dCdz{i};
    dCdzdown = FCTDdehyst.down.dCdz{i};

    tmax = min(max(Cup),max(Cdown));
    tmin = max(min(Cup),min(Cdown));

    if tmax<=tmin
        Js(i)=NaN;
        disp(['Bad set of profiles at index ' num2str(i)]);
        continue;
    end
    
    trange = tmax-tmin;
    tlim=[tmin+trange/12 tmax-trange/12]; % Use the central 5/6 of the available conductivity range.

    dClim=0.5; % Sometimes weird stuff happens - want to eliminate spurious conductivity steps over 0.5
    
    I = find(Cup>tlim(1) & Cup<tlim(2) & abs(dCdzup)<dClim);
    Cup = Cup(I);
    dCdzup = dCdzup(I);
    [Cup,I] = unique(Cup);
    dCdzup = dCdzup(I);

    I = find(Cdown>tlim(1) & Cdown<tlim(2) & abs(dCdzdown)<dClim);
    Cdown = Cdown(I);
    dCdzdown = dCdzdown(I);
    [Cdown,I]=unique(Cdown);
    dCdzdown = dCdzdown(I);


    if length(Cdown)<2||length(Cup)<2
        Js(i)=NaN;
        disp(['Bad set of parameters for profiles at index ' num2str(i)]);
        continue;
    end
    
    % Make a uniform temeprature grid the same length as the longest record
    %tgrid = linspace(tlim(1),tlim(2),max(length(Cdown),length(Cup))*10);
    
    % Make a non-uniform grid
    tgrid = [Cup; Cdown];
    tgrid = unique(tgrid);

    dCdzup = interp1(Cup,dCdzup,tgrid);
    dCdzdown = interp1(Cdown,dCdzdown,tgrid);
    mse = mean((dCdzup-dCdzdown).^2,'omitnan');
    Js(i)=mse;
    
    toplot=0;
    if toplot
    % Test plots
        figure(1);
        clf(1);
        hold on;
        plot(tgrid,dCdzup);
        plot(tgrid,dCdzdown);
        legend('Up','Down');
        xlabel('Conductivity');
        ylabel('dCdz');
        title(sprintf('Profile %i of %i',i,nprof));
        pause(0.3);
    end
        
end
J = mean(Js,'omitnan');
end
