clc
clear all
close all
% 原始数据
[S,F]=rdsamp('f:\Database\CUDB\cu05\cu05',[],(60*8+28)*250);
t=S((60*8+23)*250+1:end,1);
f=F((60*8+23)*250+1:end,1);
% fs = 1000;
% ts = 1/fs;
% t=0:ts:0.3;
% z=2*sin(2*pi*10*t) + 5*sin(2*pi*100*t);
figure
% plot(t, z)
plot(t,f)
title('原始信号')
% EMD
% imf=emd(z);
% emd_visu(z,t,imf)
imf=emd(f);
emd_visu(f,t,imf)
[A,f,tt]=hhspectrum(imf);
[im,tt]=toimage(A,f);
disp_hhs(im);