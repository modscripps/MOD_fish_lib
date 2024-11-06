function Chat = deHyst_C(T,C,aC,bC,gamma,lag)
% This function calculates the true conductivity, Chat, from the true temperature, measured conductivity, and some parameters, following Lueck & Picklo
% 1990
Chat = zeros(length(C),1);  
for t=max(1+lag,2):length(T)
    Chat(t) = -bC*Chat(t-1)+gamma*aC*(T(t)-T(t-lag))/lag;
end
Chat = C+Chat;
Chat(1:max(1+lag,2)-1)=Chat(max(1+lag,2));
end