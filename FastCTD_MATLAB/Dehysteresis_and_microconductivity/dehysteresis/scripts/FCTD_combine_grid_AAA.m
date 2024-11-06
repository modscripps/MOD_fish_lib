function FCTDgrid = FCTD_combine_grid_AAA(FCTD,datapath)

% Define structure and transfer auxilary data
FCTDgrid=[];

P=inputParser;
addRequired(P,'FCTD',@isstruct);

validationFcn = @(x) validateattributes(x,{'char','string'},{});
addRequired(P,'datapath',validationFcn);

parse(P,FCTD,datapath);
FCTD = P.Results.FCTD;
datapath = P.Results.datapath;


vars = fieldnames(FCTD.up);
vars(strcmp(vars,'depth'))=[];
for l=1:length(vars)
    varname = vars{l};
    FCTDgrid.(varname)= cat(2,FCTD.up.(varname),FCTD.down.(varname));
end
tmean = mean(FCTDgrid.time,1,'omitnan');
[~,I] = sort(tmean);
for l=1:length(vars)
    varname = vars{l};
    FCTDgrid.(varname) = FCTDgrid.(varname)(:,I,:);
end
FCTDgrid.depth = FCTD.up.depth;
FCTDgrid.info=FCTD.info;

vars = fieldnames(FCTD.header);
for l=1:length(vars)
    varname = vars{l};
    tmp = cat(2,FCTD.header.(varname),FCTD.header.(varname));
    FCTDgrid.header.(varname) = tmp(I);
end

save(datapath,'-struct','FCTDgrid','-v7.3');
disp(['Saved Gridded FCTD profiles in ' datapath]);