function varhat = deHyst_Offset_P(var,offset)
% This function adjusts a variable vector by an offset
varhat = NaN(size(var));
if offset>0
    varhat(1:end-offset)=var(1+offset:end);
    varhat(end-offset+1:end)=varhat(end-offset);
elseif offset<0
    varhat(1-offset:end)=var(1:end+offset);
    varhat(1:-offset)=varhat(1-offset);
else
    varhat = var;
end
end