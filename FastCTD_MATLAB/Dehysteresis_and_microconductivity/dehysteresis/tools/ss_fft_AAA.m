function [f,a]=ss_fft_AAA(x,dt)
% ss_fft_AAA creates the single-sided fft of a timeseries x
%
% Alex Andriatis
% 2021-08-22

xft = fft(x);
p = length(xft);
if mod(p,2)==0
    k=0:p/2; %Use only half of the fourier transform since it's symetric about its middle value
else
    k=0:(p-1)/2;
end
T=p*dt; %Total length of record
f = k/T; %Frequency axis
if mod(p,2)==0
   a=xft(1:p/2+1);
else
   a=xft(1:(p+1)/2);
end
end