function trend = trend_AAA(x,y,n)
% trend_AAA calculates the trend of a dataset y using an nth degree
% polynomial least-squares fit. The x-vector can be omitted if the data is
% evenly spaced
%
% Examples:
% trend = trend_AAA([],y,0) % Removes a mean from evenly-spaced data y
% trend = trend_AAA(x,y,1) % Removes a linear trend from data y
%
% Alex Andriatis
% 2021-08-23

T=length(y);
if isempty(x)
    x=1:T;
end
p = polyfit(x,y,n);
trend = polyval(p,x);
trend = reshape(trend,size(y));
end

