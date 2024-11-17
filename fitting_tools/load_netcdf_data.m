function [AllData,AllAtts,VarName,newName] = load_netcdf_data(filename,nlist)
%function [Data,AllAtts, varName] = load_netcdf(filename,nlist)
%  Load the NetCDF  data in "filename" into a matlab struct
%  Can optionally request a subset of data.
% Handles "flat" NetCDF or one-level deep NetCDF groups (ATOMIX-SCOR format)
% Optional Input:
%   nlist: cell array of desired variables if known... 
%           Use  ncdisp('name_of_file.nc','/') to identify variables and
%           avoid loading too much data. 
%           Alternatively: ncdisp('name_of_file.nc','GroupName/') for
%           variables stored in GroupName
%   It now works with flat NetCDFs and hierarchical NetCDF (1 level deep)
%       e.g., nlist={'Level1/TIME','Level1/PRES','Level2'} will load
%       TIME & PRES variables in group Level1, and ALL variables in Level2.
%       To read "ALL" Global level data,  specify 'Global'.
%       IF you want ONE Global variable (e.g. TIME), then you can specify
%       'TIME' OR 'Global/TIME'
% Outputs:
%   AllData: structure of the data in NetCDF in the same hierarchal format.
%        If NetCDF is flat, then each fields will have the same names as
%           the short_name variables in the NetCDF
%       If NetCDF has grouped data, then each field represent a structure
%       with the group's name. The fields in each group's structure
%       represent the variable field1 (e.g. AllData.group1.field1)
%   AllAtts: structure with fields  "global" and group names, which details
%   all the attributes associate with the NetCDF file
%   VarName: structure (grouped NetCDF) or cell array (flat NetCDF) of
%       possible variable names that can be loaded.
%   newName: is cell array (nvars x 3columns) with the entire list
%       The 1st column is the group, 2nd column the var, 3rd is properly
%       merged tother group/var
% Created by CBluteau in Oct 2017 when I couldn't open a netcdf file in
% ncview/ncbrowse.
% Modifications
%   2021/09/28. CB accomodated requesting subset of variables from NetDCF groups.
%%%
% Initialise/defaults
if nargin<2
    nlist=[];
end

AllData=struct;
VarName=struct;
%%

[gpath,filename]=fileparts(filename);
filename=fullfile(gpath,[filename,'.nc']);

nInfo=ncinfo(filename);
ncid = netcdf.open(filename,'NOWRITE');
gInfo=nInfo.Groups;

[Att,GroupAtt] = load_netcdf_attributes(nInfo);

%%
AllAtts.Global=Att;

gid=ncid;
gName{1}='Global';
if ~isempty(gInfo)
    nG=length(gInfo);
    AllAtts=append_struc(AllAtts,GroupAtt); 
    for kk=1:nG
        gName{kk+1}=gInfo(kk).Name;
        gid(kk+1) = netcdf.inqNcid(ncid,gName{kk+1});
    end
end

%% Loop
nG=length(gid);

for kk=1:nG
    disp(['Trying to load group data: ',gName{kk}])
    varids = netcdf.inqVarIDs(gid(kk));
    
    if isempty(varids)
        continue;
    end
    
  
    for ii=1:length(varids)
        VarName.(gName{kk}){ii}=netcdf.inqVar(gid(kk),varids(ii));
    end
    
    if strcmp(gName{kk},'Global')
        if any(strcmpi(nlist,'Global'))
            group_list=[]; % load the entire top (Global) list of vars
        else
            [extra_list]=split_list(nlist,gName{kk}); % handles the format Global/varname
            group_list=[nlist extra_list];% assign nlist and those that are "Global" will be loaded. Redundant variables, but the other checks will deal with this "bug
        end
    else
        [group_list]=split_list(nlist,gName{kk});
    end
    
    
    if isempty(group_list)
        nVar= VarName.(gName{kk});
    else
        if isempty(group_list{1})==1 % no variables desired
            % VarName=rmfield(VarName,gName{kk});
            continue;
        else
            [nVar,varids]=getvarName(VarName.(gName{kk}),varids,group_list);
        end
    end
      
    for ii=1:length(varids)
        if strcmp(gName{kk},'Global')==0
            tmp = netcdf.getVar(gid(kk),varids(ii));
            AllData.(gName{kk}).(nVar{ii})=convert_char_cell(tmp);
           
        else
            tmp = netcdf.getVar(gid(kk),varids(ii));
            AllData.(nVar{ii})=convert_char_cell(tmp);
        end
    end
      
    clear Data;
end


%% Assign possible varnames to newNAmes and struct VarName
if  all([nG==1 strcmp(gName{1},'Global')])
    % %      AllData=AllData.(gName{1});
    VarName=VarName.(gName{1});
    newName=VarName;
else
    
    vF=fieldnames(VarName);
    cc=0;
    for kk=1:length(vF)
        tmpName=VarName.(vF{kk});
        nF=length(tmpName);
        for ii=1:nF
            cc=cc+1;
            newName{cc,1}=vF{kk};
            newName{cc,2}=tmpName{ii};
             switch vF{kk}
                 case{'Global'}
                     newName{cc,3}=tmpName{ii}; % didn't want the word Global
                 otherwise
                     newName{cc,3}=[vF{kk},'/',tmpName{ii}];
             end
        end
    end
end


netcdf.close(ncid)
end % end of main

%% subfcts
function [nVar,vIds]=getvarName(varName,varids,nlist)
cc=0;
for ii=1:length(nlist)
    ind=find(strcmp(nlist{ii},varName));
   
    if isempty(ind)
        warning(['Variable ', nlist{ii},' doesnt exist for this group']);
    else
        cc=cc+1;
        nVar{cc}=varName{ind};
        vIds(cc)=varids(ind);
    end
end

end



%% Process list of vars into groups & variables
function [Grouplist]=split_list(nlist,gName)
% Sorts out a new variable list for the specified group in gName
% If nlist=[] assumes the entire list of names is desired.
% If no vars in nlist are associated with  gName, then no data/var from
% that group will be loaded.

if isempty(nlist)
    Grouplist=[]; % All var
    return;
end


nVar=length(nlist);
cc=0;
tempty=0;
for ii=1:nVar
    tmp=strsplit(nlist{ii},'/'); % see if groups & var specified
    
    if strcmp(tmp{1},gName)==1
        cc=cc+1;
        if length(tmp)>1 % group/var specified
            Grouplist{cc}=tmp{2};
        else
            Grouplist{cc}=[];
            tempty=1;
        end
    else
        continue;
    end
end

if tempty
    Grouplist=[]; % requesting all variables in this group
end

if cc==0
    Grouplist=cell(1);% no variables are wanted
end
end

%%
function dat=convert_char_cell(tmp)
% Converting characters arrays as cells
 if ischar(tmp)
     dat=char2cell(tmp);
 else
     dat=tmp;
 end

end
