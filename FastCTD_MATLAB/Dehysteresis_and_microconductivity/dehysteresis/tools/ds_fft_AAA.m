function [f,a]=ds_fft_AAA(x,dt)
% ss_fft_AAA creates the double-sided fft of a timeseries x
%
% Alex Andriatis
% 2021-08-22

xft = fft(x);
p = length(xft);
if mod(p,2)==0
    k=0:p/2; %Use only half of the fourier transform since it's symetric about its middle value
    k=[k flip(k(2:end-1))];
else
    k=0:(p-1)/2;
    k=[k flip(k(2:end))];
end
T=p*dt; %Total length of record
f = k/T; %Frequency axis
a = xft;
end