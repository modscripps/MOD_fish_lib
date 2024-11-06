function That = deHyst_T(T,TParams)
% This function calculates the true temperature, That, from the measured
% temperature, T, and some parameters, TParams, following Lueck & Picklo
% 1990
That = T;
for t=1+TParams(1):length(T)
    That(t) = (T(t)-TParams(2)*T(t-TParams(1)))/(1-TParams(2));
end
That(1:TParams(1))=That(1+TParams(1));
end