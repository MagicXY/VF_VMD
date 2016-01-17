clc
clear all
close all
% ԭʼ����
fs = 1000;
ts = 1/fs;
t=0:ts:0.3;
y=2*sin(2*pi*10*t);% + 5.*sin(2*pi*100*t);
figure
plot(t, y)
title('ԭʼ�ź�')
% ��Hilbert-Huang��
[A,fh,th] = hhspectrum(y);
figure
subplot(211)
plot(th*ts, A)
title('˲ʱ��ֵ') % ���ǰ���
subplot(212)
plot(th*ts, fh*fs)
title('˲ʱƵ��')
% ��ʾ���
[im,tt,ff] = toimage(A,fh,th);
disp_hhs(im,tt)
colormap(flipud(gray))
% ���ʵ����ʾ
figure
imagesc(tt*ts,[0,0.5*fs],im);
ylabel('frequency/Hz')
set(gca,'YDir','normal')
xlabel('time/s')
title('Hilbert-Huang spectrum')
