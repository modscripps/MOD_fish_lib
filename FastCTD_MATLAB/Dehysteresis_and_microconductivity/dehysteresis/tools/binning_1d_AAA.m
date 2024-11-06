function [x_center,y_binned]=binning_1d_AAA(x,y,bins,varargin)
% This fucntion is designed to bin data 
%
% Alex Andriatis
% 2021-08-13
%

P=inputParser;
addRequired(P,'x',@isnumeric);
addRequired(P,'y',@isnumeric);
if size(x)~=size(y)
    error('Input x and y must be the same size');
end
addRequired(P,'bins',@isnumeric);

defaultEdges = 'center';
checkString=@(s) any(strcmp(s,{'center','edges'}));
addParameter(P,'BinCentering',defaultEdges,checkString);

defaultInterpFill = 'none';
checkString=@(s) any(strcmp(s,{'none','interp'}));
addParameter(P,'InterpFill',defaultInterpFill,checkString);

parse(P,x,y,bins,varargin{:});
x = P.Results.x;
y = P.Results.y;
bins = P.Results.bins;
BinCentering = P.Results.BinCentering;
InterpFill = P.Results.InterpFill;

use_edges=0;
if exist('option','var')
    if strcmp(option,'edges')
        use_edges=1;
    end
end

% Make sure x,y and bins are all column vectors
x = column_AAA(x);
y = column_AAA(y);
bins = column_AAA(bins);

[x,I]=sort(x);
y = y(I);

if strcmp(BinCentering,'center')
    edges(1)=bins(1);
    edges(2:length(bins))=bins(1:end-1)+diff(bins)/2;
    edges(length(bins)+1)=bins(end);
    x_center = bins;
elseif strcmp(BinCentering,'edges')
    edges = bins;
    x_center = edges(1:end-1)+diff(edges)/2;
else
    error('Unrecognized bin option');
end
bin_ind = discretize(x,edges);
I=~isnan(bin_ind);
y_binned = accumarray(bin_ind(I),y(I),[length(bins),1],@nanmean,NaN);

if strcmp(InterpFill,'interp') && sum(~isnan(y_binned))>=2
    y_binned = interp1(x_center(~isnan(y_binned)),y_binned(~isnan(y_binned)),x_center);
end
end






