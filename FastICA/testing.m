clc,clear;
N=512;  %??
Fs=10000;  %????10khz
n=0:N-1;
t=n/Fs;   %????
%??
x=cos(2*pi*1000*t);
y=cos(2*pi*1000*(t+0.0002));

x1=hilbert(x);
y1=hilbert(y);
[cc,n]=xcorr(x1,y1);
[c,lag]=max(cc);
d=lag-N;