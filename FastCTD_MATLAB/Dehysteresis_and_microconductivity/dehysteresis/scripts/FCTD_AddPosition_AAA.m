function FCTD = FCTD_AddPosition_AAA(FCTD,GPS,datapath,varargin)
% Add position data to the profiles before gridding, useful for salinity
% calculation.
%
% If GPS structure is empty, gridding will use default values of 35 N and 0 E
%
%
% Alex Andriatis
% 2021-02-01
%

% Define structure and transfer auxilary data

P=inputParser;
addRequired(P,'FCTD',@isstruct);

addRequired(P,'GPS',@isstruct);

validationFcn = @(x) validateattributes(x,{'char','string'},{});
addRequired(P,'datapath',validationFcn);

defaultFigpath = '';
validationFcn = @(x) validateattributes(x,{'char','string'},{});
addParameter(P,'figpath',defaultFigpath,validationFcn);

defaultPlotevery = 0;
addParameter(P,'plotevery',defaultPlotevery,@isnumeric);

parse(P,FCTD,GPS,datapath,varargin{:});
FCTD = P.Results.FCTD;
GPS = P.Results.GPS;
datapath = P.Results.datapath;

directions={'up','down'};
nprof = length(FCTD.down.time);

if ~isempty(GPS)
    if ~isfield(GPS,'time') || ~isfield(GPS,'latitude') || ~isfield(GPS,'longitude')
        error('GPS structure needs a time, latitude, and longitude variable');
    end
    [~,GPS.longitude]=interp_sort_AAA(GPS.time,GPS.longitude);
    [GPS.time,GPS.latitude]=interp_sort_AAA(GPS.time,GPS.latitude);
    
    for i=1:length(directions)
        direction = directions{i};
        for k=1:nprof
            t = FCTD.(direction).time{k};
            FCTD.(direction).latitude{k}=interp1(GPS.time,GPS.latitude,t);
            FCTD.(direction).longitude{k}=interp1(GPS.time,GPS.longitude,t);
        end
    end
    save(datapath,'-struct','FCTD','-v7.3');
    disp(['Added GPS data to FCTD profiles in ' datapath]);
else
    warning('GPS structure empty, position data not added');
end
end
           
