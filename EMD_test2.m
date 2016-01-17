clc
clear all
close all
% 原始数据
fs = 1000;
ts = 1/fs;
t=0:ts:0.3;
y=2*sin(2*pi*10*t);% + 5.*sin(2*pi*100*t);
figure
plot(t, y)
title('原始信号')
% 求Hilbert-Huang谱
[A,fh,th] = hhspectrum(y);
figure
subplot(211)
plot(th*ts, A)
title('瞬时幅值') % 就是包络
subplot(212)
plot(th*ts, fh*fs)
title('瞬时频率')
% 显示结果
[im,tt,ff] = toimage(A,fh,th);
disp_hhs(im,tt)
colormap(flipud(gray))
% 编程实现显示
figure
imagesc(tt*ts,[0,0.5*fs],im);
ylabel('frequency/Hz')
set(gca,'YDir','normal')
xlabel('time/s')
title('Hilbert-Huang spectrum')
