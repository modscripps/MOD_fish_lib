function x = column_AAA(x)
% column_AAA converts data to a column vector
if size(x,2)>1
    x=x(:);
end
end