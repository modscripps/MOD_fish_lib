function J = cost_Offset_P(data,offset)
% This function calculates the cost function for the difference between up
% and down casts for a given variable in pressure space
%
% Inputs:
%   data: A structure containing pressure and a variable. Designed to work
%   with either temperature or conductivity interchangeably.
%
% Outputs:
%   J: A numerical value representing the match between the up and down
%   profile. The close it is to zero, the better they match
%
% Alex Andriatis
% 2021-02-07

dirnames = {'up','down'}; 
nprof = length(data.up.var);

for k=1:nprof
    % Apply the offset correction
    for n=1:2
        dirname = dirnames{n};
        Var = data.(dirname).var{k};
        if offset>0
            Var(1:end-offset)=Var(1+offset:end);
            Var(end-offset+1:end)=Var(end-offset);
        elseif offset<0
            Var(1-offset:end)=Var(1:end+offset);
            Var(1:-offset)=Var(1-offset);
        end
        data.(dirname).var{k} = Var;
    end
    
    Pup = data.up.pressure{k};
    Pdown = data.down.pressure{k};
    
    pmax = min(max(Pup),max(Pdown));
    pmin = max(min(Pup),min(Pdown));
    
    if pmax<pmin
         Jrecord(k)=NaN;
         continue;
    end

    prange = pmax-pmin;
    plim=[pmin+prange/12 pmax-prange/12]; % Use the central 5/6 of the available pressure range.

    for n=1:2
        dirname = dirnames{n};

        P = data.(dirname).pressure{k};
        I = find(P>plim(1) & P<plim(2));
        P = P(I);
        Var = data.(dirname).var{k}(I);

        [P,I] = unique(P);
        Var = Var(I);

        switch n
            case 1
                Pup = P;
                Varup = Var;
            case 2
                Pdown = P;
                Vardown = Var;
        end
    end

    if length(Pup)<10 || length(Pdown)<10
        Jrecord(k)=NaN;
        continue;
    end
    
    % An older version optimized the missmatch using an evenly-spaced
    % pressure vector. This ends up adding more weight to sparse
    % measurements near the top of the profile.
    
    % Instead, to make the code consistent with the sample-weighted
    % optimization that I'm now using for temperature and conductivity,
    % I'll use a non-uniform grid made out of the sample points of the
    % measurements.
    
    % pgrid = linspace(plim(1),plim(2),max(length(Pup),length(Pdown)));
    
    % Make a non-uniform grid
    pgrid = [Pup; Pdown];
    pgrid = unique(pgrid);
    
    
    Vardown = interp1(Pdown,Vardown,pgrid);
    Varup = interp1(Pup,Varup,pgrid);
    
    I = find(~isnan(Varup) & ~isnan(Vardown));
    Varup = Varup(I);
    Vardown = Vardown(I);
   
    mse = mean((Varup-Vardown).^2);
    Jrecord(k)=mse;
end
J = mean(Jrecord,'omitnan');
