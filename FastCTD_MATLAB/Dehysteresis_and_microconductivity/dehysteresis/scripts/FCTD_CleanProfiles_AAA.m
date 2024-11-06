function FCTDprofiles_cleaned=FCTD_CleanProfiles_AAA(FCTDprofiles,datapath)
%
% FCTD_CleanProfiles_AAA applies criteria to FastCTD profiles to remove bad
% profiles
%
%
% Alex Andriatis
% 2021-07-15

FCTDprofiles_cleaned=[];
%datapath='/Volumes/Andriatis_T7BL/Data/TFO/TFO2/FCTD/Processing/Profiles/FCTDprofiles_cleaned.mat';

if ~ischar(datapath)
    error('The path for saving the FCTD profiles is not a character array');
end

structnames = {'up','down','header'};

if ~all(isfield(FCTDprofiles,structnames))
    error('One or more required fields is missing from the FCTD structure');
end


nprofs = length(FCTDprofiles.header.serial);
if ~isfield(FCTDprofiles.header,'qc')
    FCTDprofiles.header.qc=zeros(1,nprofs);
end
FCTDprofiles_cleaned.header = FCTDprofiles.header;

vars = fieldnames(FCTDprofiles.up);
directions={'up','down'};

for t=1:nprofs
    disp(['Checking profile ' num2str(t) ' of ' num2str(nprofs)]);
    profile=[];
    for n=1:length(directions)
        direction=directions{n};
        for m=1:length(vars)
            var = vars{m};
            profile.(direction).(var)=FCTDprofiles.(direction).(var){t};
        end
    end
%   
    badprof=0;
    cleaned=0;
    while ~badprof && ~cleaned
        % 1. Limit on fall speed

        % I want to avoid big gaps in the timeseries. For a max fall speed of,
        % let's say, 10m/s, I expect the pressure to change by at most 10 m/s /
        % 16 samples/second = 0.625 m/sample   

        dt = abs(mean(diff(profile.up.time)*86400)); % Seconds per sample
        dPlim = 10*dt; % 10m/s * dt seconds per sample = meters per sample
        dP=diff(profile.up.pressure);
        Imid = round(length(dP)/2);
        Istart = find(dP(1:Imid)<-dPlim,1,'last');
        if isempty(Istart)
            Istart=1;
        else
            Istart = Istart+2;
        end
        Iend = find(dP(Imid+1:end)<-dPlim,1,'first')+Imid;
        if isempty(Iend)
            Iend = length(profile.up.pressure);
        else
            Iend = Iend-1;
        end
        I = Istart:Iend;
        if any(I>length(profile.up.pressure)) || any(I<1)
            badprof=1;
            continue;
        end
        if isempty(I)
            badprof=1;
            continue;
        end
        if length(I)<16
            badprof=1;
            continue;
        end
        Iup=I;

        dP=diff(profile.down.pressure);
        Imid = round(length(dP)/2);
        Istart = find(dP(1:Imid)>dPlim,1,'last');
        if isempty(Istart)
            Istart=1;
        else
            Istart = Istart+2;
        end
        Iend = find(dP(Imid+1:end)>dPlim,1,'first')+Imid;
        if isempty(Iend)
            Iend = length(profile.down.pressure);
        else
            Iend = Iend-1;
        end
        I = Istart:Iend;
        if any(I>length(profile.down.pressure)) || any(I<1)
            badprof=1;
            continue;
        end
        if isempty(I)
            badprof=1;
            continue;
        end
        if length(I)<16
            badprof=1;
            continue;
        end
        Idown=I;    

        for n=1:length(directions)
            direction=directions{n};
            switch n
                case 1
                    I=Iup;
                case 2
                    I=Idown;
            end
            for m=1:length(vars)
                var = vars{m};
                profile.(direction).(var)=profile.(direction).(var)(I,:);
            end
        end
        disp('Good fall speed');



        % 2. Cut out sections that are held at the surface. Only want to keep
        % the times when the probe is actually profiling. Set velocitiy to at
        % least 0.1m/s

        plength=length(profile.down.time);
        w=diff(profile.down.pressure)./(diff(profile.down.time)*86400); %m/s
        wlim=0.1;
        Imid = round(length(w)/2);
        Istart = find(w(1:Imid)<wlim,1,'last');
        if isempty(Istart)
            Istart=1;
        else
            Istart = Istart+2;
        end
        Iend = find(w(Imid+1:end)<wlim,1,'first')+Imid;
        if isempty(Iend)
            Iend = plength;
        else
            Iend = Iend-1;
        end
        I = Istart:Iend;
        if any(I>plength) || any(I<1)
            badprof=1;
            continue;
        end
        if isempty(I)
            badprof=1;
            continue;
        end
        if length(I)<16
            badprof=1;
            continue;
        end
        Idown=I;  

        plength=length(profile.up.time);
        w=diff(profile.up.pressure)./(diff(profile.up.time)*86400); %m/s
        wlim=-0.1;
        Imid = round(length(w)/2);
        Istart = find(w(1:Imid)>wlim,1,'last');
        if isempty(Istart)
            Istart=1;
        else
            Istart = Istart+2;
        end
        Iend = find(w(Imid+1:end)>wlim,1,'first')+Imid;
        if isempty(Iend)
            Iend = plength;
        else
            Iend = Iend-1;
        end
        I = Istart:Iend;
        if any(I>plength) || any(I<1)
            badprof=1;
            continue;
        end
        if isempty(I)
            badprof=1;
            continue;
        end
        if length(I)<16
            badprof=1;
            continue;
        end
        Iup=I;  

        for n=1:length(directions)
            direction=directions{n};
            switch n
                case 1
                    I=Iup;
                case 2
                    I=Idown;
            end
            for m=1:length(vars)
                var = vars{m};
                profile.(direction).(var)=profile.(direction).(var)(I,:);
            end
        end   
        disp('Good profiling');



        % 3. Cut out spikes in both dTdz and dCdz to eliminate false spiking



        % Sometimes the data goes wild and there's crazy spikes. Want to
        % eliminate:
        dTdzlim=5; % At most 5 degrees / meter gradients
        dCdzlim = 0.5; % At most 0.5 mS/cm / meter gradients
        % Same principle as the velocity

        for n=1:length(directions)
            direction=directions{n};
            
            despiked=0;
            tmp = profile.(direction);
        
            plength=length(tmp.time);
            while ~despiked && plength>=16
                dTdz = diff(tmp.temperature)./diff(tmp.pressure);
                dCdz = diff(tmp.conductivity)./diff(tmp.pressure);

                I=find((abs(dTdz)>dTdzlim) | (abs(dCdz)>dCdzlim),1);
                if isempty(I)
                    despiked=1;
                else
                    for m=1:length(vars)
                       var = vars{m};
                       tmp.(var)(I:I+1,:)=[];
                    end
                end
                plength=length(tmp.time);
            end
            if length(tmp.time)<16
                badprof=1;
                continue
            end
            profile.(direction)=tmp;
        end
        if badprof==1
            continue
        end
        disp('No spikes');

        % 4. Profiles must be at least 10 meters deep
        if (max(profile.up.pressure)-min(profile.up.pressure))<10 || (max(profile.down.pressure)-min(profile.down.pressure))<10
            badprof=1;
            continue;
        end
        disp('Deep enough');

        if max(min(profile.down.pressure),min(profile.up.pressure))>min(max(profile.down.pressure),max(profile.up.pressure))
            badprof=1;
            continue; % Pressures of profiles must overlap
        end
        disp('Good overlap');
     cleaned=1;
     end
    
    
    for n=1:length(directions)
        direction=directions{n};
        for m=1:length(vars)
            var = vars{m};
            FCTDprofiles_cleaned.(direction).(var){t}=profile.(direction).(var);
        end
    end
    if badprof==1
        FCTDprofiles_cleaned.header.qc(t)=1;
        warning('Bad profile');
    end
end

if length(FCTDprofiles.up.time)~=length(FCTDprofiles.down.time)
    error('Number of up and down profiles dont match')
end

disp('Saving...');
save(datapath,'-struct','FCTDprofiles_cleaned','-v7.3');
disp(['Saved cleaned FCTD profiles in ' datapath]);

