function FCTDgrid = FCTD_grid_AAA(FCTD,datapath,varargin)
% This is a test code for FCTD_grid_AAA
%
% FCTD_grid_AAA is a re-write of Rob's
% J_2016_FCTD_sal_dens.m
%
% It now works with a profile structure saved as cell arrays.
%
% Inputs:
%   FCTD: A strucure of profiles which is the output of
%   FCTD_responsematch_AAA. 
%
%1. Convert pressure to depth using GSW
%2. Find the min and max depth of all profiles
%3. Calculate salinity and potential density for each profile
%4. Grid the data onto a uniform depth grid
%5. Save variables
%
% Alex Andriatis
% 2021-02-01
%

% Define structure and transfer auxilary data
FCTDgrid=[];

P=inputParser;
addRequired(P,'FCTD',@isstruct);

validationFcn = @(x) validateattributes(x,{'char','string'},{});
addRequired(P,'datapath',validationFcn);

defaultFigpath = '';
validationFcn = @(x) validateattributes(x,{'char','string'},{});
addParameter(P,'figpath',defaultFigpath,validationFcn);

defaultPlotevery = 0;
addParameter(P,'plotevery',defaultPlotevery,@isnumeric);

default_dz = 0.5;
addParameter(P,'dz',default_dz,@isnumeric);

parse(P,FCTD,datapath,varargin{:});
FCTD = P.Results.FCTD;
datapath = P.Results.datapath;
figpath = P.Results.figpath;
plotevery = P.Results.plotevery;
dz = P.Results.dz;

directions = {'up','down'};
nprof = length(FCTD.down.time);

%% 1. Convert pressure to depth, find min and max depth
zmin = 14000;
zmax = 0;

for j=1:length(directions)
    direction = directions{j};
    data = FCTD.(direction);
    for k=1:nprof
       P = data.pressure{k};
       
       if isfield(data,'latitude')
        lat = mean(data.latitude{k},'omitnan');
        if isnan(lat)
            lat=35;
        end
       else
        lat = 35;
       end
       
       z = gsw_z_from_p(P,lat);
       depth = -z;
       FCTD.(direction).depth{k}=depth;
       
       if max(depth)>zmax
           zmax = max(depth);
       end
       if min(depth)<zmin
           zmin = min(depth);
       end
    end
end

zmin(zmin<0)=0;
zmin = floor(zmin);
zmax = ceil(zmax);
if zmax-zmin<10
    error('Max profile length less than 10 meters');
end
dgrid = zmin:dz:zmax;
dgrid = column_AAA(dgrid);

%% 3. Calculate salinity and potential density for each profile

for j=1:length(directions)
    direction = directions{j};
    data = FCTD.(direction);
    for k=1:nprof
        
        P = data.pressure{k};
        T = data.temperature{k};
        C = data.conductivity{k};
        if isfield(data,'latitude')
            lat = mean(data.latitude{k},'omitnan');
            if isnan(lat)
                lat=35;
            end
        else
            lat = 35;
        end
        if isfield(data,'longitude')
            lon = mean(data.longitude{k},'omitnan');
            if isnan(lon)
                lon=0;
            end
        else
            lon = 0;
        end       
        SP = gsw_SP_from_C(C*10,T,P);        
        SA = gsw_SA_from_SP(SP,P,lon,lat);
        CT = gsw_CT_from_t(SA,T,P);
        rho = gsw_rho(SA,CT,P);
        density = gsw_rho(SA,CT,0);
        FCTD.(direction).salinity{k}=SP;
        FCTD.(direction).rho{k}=rho;
        FCTD.(direction).density{k}=density;
    end
end


%% 4. Grid data onto a uniform depth grid

%vars2grid={'time','pressure','temperature','conductivity','salinity','latitude','longitude'};

for j=1:length(directions)
    direction = directions{j};
    data = FCTD.(direction);
    vars2grid = fieldnames(data);
    vars2grid(strcmp(vars2grid,'depth'))=[];
    vars2grid(strcmp(vars2grid,'depth'))=[];
    for i=1:length(vars2grid)
        variable = vars2grid{i};
        if ~isfield(data,variable)
            continue
        end
        disp(['Gridding variable ' variable ' in the ' direction ' direction']);
        vsize = size(data.(variable){1},2);
        tmpgrid=NaN(length(dgrid),nprof,vsize);
        for k=1:nprof
            tmp = data.(variable){k};
            if size(tmp,1)<2
                continue
            end
            depth = data.depth{k};
            for l=1:vsize
                % Grid the data using hanning averaging and linear
                % interpolation
                % For faster code, simple binning works fine
                [~,tmpgrid(:,k,l)]=binning_1d_AAA(depth,tmp(:,l),dgrid,'BinCentering','center','InterpFill','interp');
            end
        end
        tmpgrid=squeeze(tmpgrid);
        FCTDgrid.(direction).(variable)=tmpgrid;
    end
    FCTDgrid.(direction).depth = dgrid;
    FCTDgrid.(direction).time_mean=mean(FCTDgrid.(direction).time,1,'omitnan');
end

%% 6. Save variables
info.time = 'Time, [Matlab Datenum]';
info.time_mean = 'Average Profile Time, [Matlab Datenum]';
info.pressure = 'Pressure, [dbar]';
info.temperature = 'In-situ Temperature, [deg. C]';
info.conductvitiy = 'Conductivity, [S/m]';
info.salinity = 'Practical Salinity, [g/kg]';
info.rho = 'In-situ density, [kg/m^3]';
info.density = 'Potential density referenced to the surface, [kg/m^3]';
info.longitude = 'Longitude, [Decimal degrees, -180:180]';
info.latitude = 'Longitude, [Decimal degrees, -90:90]';
info.depth = 'Depth, [m from surface, positive down]';

FCTDgrid.info = info;
FCTDgrid.header = FCTD.header;

save(datapath,'-struct','FCTDgrid','-v7.3');
disp(['Saved Gridded FCTD profiles in ' savepath]);