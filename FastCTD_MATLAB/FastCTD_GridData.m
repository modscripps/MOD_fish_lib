function DataGrid = FastCTD_GridData(FCTD,varargin)
% Version of FastCTD_GridData used on TFO Seamounts. Adds uConductivity and
% fluorometer data (lines 255- )
%
%  DataGrid = FCTDMakeGrid(FCTD);
%
%  INPUTS:
%       FCTD - structure with timeseries of time, pressure, temperature,
%       etc
%
%   Temperature, Conductivity, Pressure, Salinity and Potential Density
%   will be placed on a grid of temperature versus depth
%
%  DataGrid = FCTDMakeGrid(FCTD,downcast,todo);
%   in addition to the variables described above, one can specify other
%   variable by defining TODO using a cell structure of strings
%
%  DataGrid is a structure
%
% Written by Jody Klymak
% Updated 2011 07 14 by San Nguyen
% used 20 point median filter to smooth the pressure

%26 Aug 2024: MHA has redone the short term response matching code.

if ~isfield(FCTD,'pressure')
    disp(FCTD);
    DataGrid=[];
    return;
end

% Get the field names that match the size of FCTD.time
vars2Grid = get_FCTD_fields(FCTD);

pickDownCast = true;
zInterval = 0.5;
zMin = 0;
zMax = 2000;

persistent argsNameToCheck;
if isempty(argsNameToCheck)
    argsNameToCheck = {'VarsToGrid','upcast','downcast','zMin','zMax','zInterval'};
end

index = 1;
n_items = nargin-1;
while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        error('MATLAB:FastCTD_GridData:wrongOption','Incorrect option specified: %s',varargin{index});
    end
    
    switch i
        case 1 % varsToGrid
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
             vars2Grid = varargin{index+1};

             % NC 11/9/24 - Commenting this out because it makes my
             % vars2Grid list error out.
             %
            % if iscellstr(vars2Grid)
            %     for i = 1:length(vars2Grid)
            %         if ~isfield(FCTD,vars2Grid{i})
            %             error('MATLAB:FastCTD_GridData:wrongVar2Grid','Wrong variable to grid: %s', vars2Grid{i});
            %         else
            %             switch lower(vars2Grid{i})
            %                 case vars2Grid_default
            %                     vars2Grid{i} = lower(vars2Grid{i});
            %                     continue;
            %                 otherwise
            %                     error('MATLAB:FastCTD_GridData:wrongVar2Grid','Wrong variable to grid: %s', vars2Grid{i});
            %             end
            %         end
            %     end
            % elseif ischar(vars2Grid)
            %     switch lower(vars2Grid)
            %         case vars2Grid_default
            %             vars2Grid = {lower(vars2Grid)};
            %             continue;
            %         otherwise
            %             error('MATLAB:FastCTD_GridData:wrongVar2Grid','Wrong variable to grid: %s', vars2Grid);
            %     end
            % else
            %     error('MATLAB:FastCTD_GridData:wrongVar2Grid','Variable to grid must be specified as a cell strings of variables: pressure, temperature, conductivity');
            % end
            
            index = index +2;
            n_items = n_items-2;
        case 2 % upcast
            pickDownCast = false;
            
            index = index + 1;
            n_items = n_items - 1;
        case 3 % downcast
            pickDownCast = true;
            
            index = index + 1;
            n_items = n_items - 1;
        case 4 % zMin
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            zMin = varargin{index+1};
            if ~isnumeric(zMin)
                error('MATLAB:FastCTD_GridData:zMinNumeric','zMin must be numeric');
            elseif length(zMin) > 1 || isempty(zMin)
                error('MATLAB:FastCTD_GridData:zMinScalar','zMin must be a scalar');
            end
            index = index +2;
            n_items = n_items-2;
        case 5 % zMax
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            zMax = varargin{index+1};
            if ~isnumeric(zMax)
                error('MATLAB:FastCTD_GridData:zMaxNumeric','zMin must be numeric');
            elseif length(zMax) > 1 || isempty(zMax)
                error('MATLAB:FastCTD_GridData:zMinScalar','zMax must be a scalar');
            end
            index = index +2;
            n_items = n_items-2;
        case 6 % zInterval
            if n_items == 1
                error('MATLAB:FastCTD_GridData:missingArgs','Missing input arguments');
            end
            zInterval = varargin{index+1};
            if ~isnumeric(zInterval)
                error('MATLAB:FastCTD_GridData:zMaxNumeric','zMin must be numeric');
            elseif length(zInterval) > 1 || isempty(zInterval)
                error('MATLAB:FastCTD_GridData:zMinScalar','zMax must be a scalar');
            elseif zInterval <= 0
                error('MATLAB:FastCTD_GridData:zMinZero','zInterval must be greater than zero');
            end
            index = index +2;
            n_items = n_items-2;
    end
end

if zMin > zMax
    zTemp = zMin;
    zMin = zMax;
    zMax = zTemp;
end
clear zTemp;

if zInterval > (zMax-zMin)
    error('MATLAB:FastCTD_GridData:zIntervalTooBig','zInterval must less than the range of z');
end

if pickDownCast
    FCTD = FastCTD_FindCasts(FCTD);
else
    FCTD = FastCTD_FindCasts(FCTD,'upcast');
end

if ~isfield(FCTD,'drop')
    DataGrid = [];
    return;
end

zMin = zMin - zInterval/2;
zMax = zMax + zInterval/2;
DataGrid.depth=(zMin:zInterval:zMax)';



% % time correction for the Conductivity cell
% df = 2.2;
% dpha = 1.05;
%
% dt = dpha/df/2/pi;
% t = FCTD.time-dt/24/3600;
% conductivity = FCTD.conductivity;
%
% good = find(FCTD.time>0);
% if isempty(good);
%     DataGrid = [];
%     return;
% end;
%
% [it,ind] = unique(FCTD.time(good));
% good = good(ind);
% if length(good) > 3
%     conductivity(good) = interp1(FCTD.time(good),FCTD.conductivity(good),t(good));
% end
% FCTD.conductivity = conductivity;

drops = unique(FCTD.drop);
drops = drops(drops>0);

num = 0;
for i=1:length(drops)
    ind = find(FCTD.drop==drops(i));
    if max(FCTD.pressure(ind))-min(FCTD.pressure(ind))>10
        num = num+1;
    end
end

% allocate space to grid data
DataGrid.time = NaN(1,num);
for i = 1:length(vars2Grid)
    if ~strcmp(vars2Grid{i},'depth')
        DataGrid.(vars2Grid{i}) = NaN(length(DataGrid.depth)-1,num);
    end
end

% Add longitude and latitude
if isfield(FCTD,'GPS')
    FCTD.longitude = FCTD.GPS.longitude;
    FCTD.latitude = FCTD.GPS.latitude;
end

% Loop through the drops and concatenate them
num = 0;
for i=1:length(drops)
    ind = find(FCTD.drop==drops(i));
    if max(FCTD.pressure(ind))-min(FCTD.pressure(ind))>10 
        for j=1:length(vars2Grid)
            %ALB add try/catch because it break in the fluor during TFO seamount 2023
            try
                myFCTD.(vars2Grid{j}) = FCTD.(vars2Grid{j})(ind);
            catch
                myFCTD.(vars2Grid{j})=nan.*ind;
            end
            % ALB end of hack
        end
        %         myFCTD.temperature = FCTD.temperature(ind);
        %         myFCTD.pressure = FCTD.pressure(ind);
        %         myFCTD.conductivity = FCTD.conductivity(ind);


        num = num+1;
        DataGrid.time(num) = mean(FCTD.time(ind),'omitmissing');
        for j=1:length(vars2Grid)
            if ~strcmp(vars2Grid{j},'depth') %don't try to grid depth. We already have it
            DataGrid.(vars2Grid{j})(:,num) = bindata1d(DataGrid.depth,...
                myFCTD.depth, ...
                myFCTD.(vars2Grid{j}));
            end
        end
    end
end

DataGrid.depth = midpoints(DataGrid.depth);

if ~isfield(DataGrid,'temperature') || ~isfield(DataGrid,'pressure') || ~isfield(DataGrid,'conductivity')
    return;
end

DataGrid.salinity_despike = sw_salt(DataGrid.conductivity*10/sw_c3515,DataGrid.temperature,DataGrid.pressure);
fc = 1./10;
fs = 1./0.5;
[b,a] = cheby2(4,20,fc/(fs/2));
for p=1:length(DataGrid.time)
    nanmask=~(isnan(DataGrid.salinity_despike(:,p)));
    % NC 3/19/25 - this breaks if there isn't enough data. "Data length
    % must be greater than 12, which is the maximum between three times the
    % filter order 4 and 1." I don't know if the filter order is always 12,
    % but I'll use that
    if sum(nanmask)>12
        filtScorr=DataGrid.salinity_despike(:,p);
        filtScorr(nanmask)=filtfilt(b,a,DataGrid.salinity_despike(nanmask,p));
        DataGrid.salinity(:,p)=filtScorr;
    end
end
DataGrid.density = sw_pden(DataGrid.salinity,DataGrid.temperature,DataGrid.pressure,0);

% gridding in time
mintime = min(DataGrid.time,[],'omitmissing');
maxtime = max(DataGrid.time,[],'omitmissing');

DataGrid.tGrid.time = mintime:2*median(diff(DataGrid.time),'omitmissing'):maxtime; % every minute
DataGrid.tGrid.depth = DataGrid.depth;


% allocate space to grid data
for i = 1:length(vars2Grid)
    DataGrid.tGrid.(vars2Grid{i}) = NaN(length(DataGrid.tGrid.depth),length(DataGrid.tGrid.time)-1);
end

for i = 1:length(DataGrid.tGrid.depth)
    for j = 1:length(vars2Grid)
        if ~strcmp(vars2Grid{j},'depth')
        DataGrid.tGrid.(vars2Grid{j})(i,:) = bindata1d(DataGrid.tGrid.time,...
            DataGrid.time, DataGrid.(vars2Grid{j})(i,:));
        end
    end
end

DataGrid.tGrid.salinity = sw_salt(DataGrid.tGrid.conductivity*10/sw_c3515,DataGrid.tGrid.temperature,DataGrid.tGrid.pressure);
DataGrid.tGrid.density = sw_pden(DataGrid.tGrid.salinity,DataGrid.tGrid.temperature,DataGrid.tGrid.pressure,0);

DataGrid.tGrid.time = midpoints(DataGrid.tGrid.time);

return;
end
