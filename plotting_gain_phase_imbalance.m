close all;
clear;
fm = 10;  % Frequency in Hz (single tone)
fs = 1000;  % Sampling frequency in Hz
t = 0:1/fs:1-1/fs;
fc=50;
%plotting phase and gain impalance
i=0;
phi=0;

for g = logspace(-3,3,10000)
    i=i+1;
    I = cos(2*pi*fc*t);
    Q = g * sin (2*pi*fc*t - phi);
    out = I + 1j * Q;
    w = fftshift(abs(fft(out)));
    IQ_Imbalance(i) = 20*log10(max( w(1:500))/(max( w(500:1000))));
end
i=0;
g = 1;
for phi = -0.2:0.001:0.2
    i=i+1;
I = cos(2*pi*fc*t);
Q = g * sin (2*pi*fc*t - phi);
out = I + 1j * Q ;
w=fftshift(abs(fft(out)));
IQ_Imbalancep(i) = 20*log10(max( w(1:500))/(max( w(500:1000))));
end
subplot(1,2,1)
plot(20*log10(logspace(-3,3,10000)),IQ_Imbalance)
title('IQ Imbalance gain')
xlabel('gain imbalance , dB');
ylabel('Image Rejection ,dB');
grid on
xlim([-3 3])
ylim([-50 -10])
subplot(1,2,2)
plot((-0.2:0.001:0.2)*180/pi,IQ_Imbalancep)
title('IQ Imbalance phase')
xlabel('Phase imbalance , dB');
ylabel('Image Rejection ,dB');
grid on
xlim([-10 10])
ylim([-50 -20])