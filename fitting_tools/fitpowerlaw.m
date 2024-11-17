function [A b Aci bci R2]=fitpowerlaw(xi,yi,bf,alp)
%[A b Aci bci R2]=fitpowerlaw(xi,yi,bf,alp)
% Use least-sqaures analysis to fit a power law y=Ax^b relationship without
% resorting to log-transformation and linear least-squares
%Required inputs:
%   xi and yi: dependent and independent "sample" data to perform powerlaw
%   bf= power relationship if forcing it (useful for some turbulence relationships to find the "best" A constant. 
%		Set to [] if want to find the power b from least-sqaures
% Optional inputs:
%   alp: confidence intervals for coefficients A and b 100*(1-alp), 90% by default
% Outputs:
%   A and b:  fitted coefficients to xi and yi data pairs
% Returns the confidence interval of the fitted parameters...
%   Aci= conf interval for A
%   bci: confidenc interval for the power b
%   R: sqrt of the R2 of the goodness of fit

if nargin<4
    alp=0.10;% (1-alp)*100=90% confidence intervals for variables A & b
end


ind=find(~isnan(yi) & ~isnan(xi));
A=NaN;b=NaN;
Aci=[NaN NaN];
bci=[];
R2=[];
if isempty(ind)
    return;
end
    
yi=yi(ind);
xi=xi(ind);

   %%

[p,S]=polyfit(log10(xi),log10(yi),1);
p(2)=10.^(p(2));

if isempty(bf)
    myfun=@(betpar,xi)(betpar(2).*xi.^betpar(1)); 
else
    b=bf; %desired power law
    myfun=@(betpar,xi)(betpar.*xi.^b); 
    p=p(2); 
end

%Opt.TolX=1e-8;Opt.TolFun=1e-8;
[bet,r,J,COVB]  = nlinfit(xi,yi,myfun,p);%,Opt);

%% Error stuff

ci = nlparci(bet,r,'covar',COVB,'alpha',alp); % 90% ci

if isempty(bf)
    A=bet(2);
    b=bet(1); 
    Aci=ci(2,:);
    bci=ci(1,:);

    sse = sum(r.^ 2);
    sst=sum((yi-mean(yi)).^2);
    R2=(1-(sse./sst)); 
else
    A=bet;
    Aci=ci(1,:);

    sse = sum(r.^ 2);
    sst=sum((yi-mean(yi)).^2);
    R2=(1-(sse./sst)); 
 end

