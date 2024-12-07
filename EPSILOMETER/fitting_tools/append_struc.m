function StruR=append_struc(StruR,StrutoAppend,vF)
% function grabs all the fields in StrutoAppend, and rewrites/appends them
% to structure StruR.
% Optional
%   vF: cell array of fields to append. If noneprovided or empty, all the
%   fields in StruToAppend will be added
if nargin<3
    vF=fieldnames(StrutoAppend);
end

if isempty(vF)
    vF=fieldnames(StrutoAppend);
end

%%
for ii=1:length(vF)
    if isfield(StruR,vF{ii})
        disp(['Field ',vF{ii},' exists, so skipping'])
    else
        if isfield(StrutoAppend,vF{ii})
            StruR.(vF{ii})=StrutoAppend.(vF{ii});
        else
            disp(['Field ',vF{ii},' doesnt exists'])
        end
    end
end
