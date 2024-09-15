close all;
clear;
% Parameters
fm = 10;  % Frequency in Hz (single tone)
fs = 1000;  % Sampling frequency in Hz
t = 0:1/fs:1-1/fs;  % Time vector from 0 to 1 second

% Complex signal in time domain
signal = exp(1j*2*pi*fm*t);

% Fourier Transform (Frequency Domain)
N = length(signal);  % Number of points
fft_signal = fft(signal);  % Compute FFT

% Frequency vector
f = (-N/2:N/2-1)*(fs/N);  % Frequency range from -fs/2 to fs/2

% Shift zero frequency component to center of spectrum
fft_signal_shifted = fftshift(fft_signal);

% Magnitude of FFT
magnitude_fft = abs(fft_signal_shifted);

% Plot the magnitude of FFT (Frequency Domain)
figure;
plot(f, magnitude_fft);
title('Single Tone Signal in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

fc=50;
modulated_signal=real(signal).*cos(2*pi*fc*t)-imag(signal).*sin(2*pi*fc*t);
N_modulated = length(modulated_signal);  % Number of points
fft_modulated = fft(modulated_signal);  % Compute FFT

% Frequency vector
f_modulated = (-N_modulated/2:N_modulated/2-1)*(fs/N_modulated);  % Frequency range from -fs/2 to fs/2

% Shift zero frequency component to center of spectrum
fft_modulated_shifted = fftshift(fft_modulated);

% Magnitude of FFT
magnitude_fft_modulated = abs(fft_modulated_shifted);

% Plot the magnitude of FFT (Frequency Domain)
figure;
plot(f_modulated, magnitude_fft_modulated);
title('Modulated Single Tone Signal in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

phi=0;
g=1.1266;
demodulated_signal_1 = modulated_signal.*cos(2*pi*fc*t);
demodulated_signal_2 = modulated_signal.*(-1*g*sin(2*pi*fc*t+phi));
demod = demodulated_signal_1+1j*demodulated_signal_2;
N_demod = length(demod);  % Number of points
fft_demod = fft(demod);  % Compute FFT

% Frequency vector
f_demod = (-N_demod/2:N_demod/2-1)*(fs/N_demod);  % Frequency range from -fs/2 to fs/2

% Shift zero frequency component to center of spectrum
fft_demod_shifted = fftshift(fft_demod);

% Magnitude of FFT
magnitude_fft_demod = abs(fft_demod_shifted);

% Plot the magnitude of FFT (Frequency Domain)
figure;
plot(f_demod, magnitude_fft_demod);
title('demod Single Tone Signal in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

cutoff_freq = 20;  % Cutoff frequency for the LPF (slightly above fm)
[b, a] = butter(5, cutoff_freq/(fs/2));  % 5th order Butterworth filter
% Apply the low-pass filter
I_filtered = filtfilt(b, a, demodulated_signal_1);  % Filtered in-phase component
Q_filtered = filtfilt(b, a, demodulated_signal_2);  % Filtered quadrature component
baseband_signal = I_filtered + 1j * Q_filtered;

N_baseband = length(baseband_signal);  % Number of points
fft_baseband = fft(baseband_signal);  % Compute FFT

% Frequency vector
f_baseband = (-N_baseband/2:N_baseband/2-1)*(fs/N_baseband);  % Frequency range from -fs/2 to fs/2

% Shift zero frequency component to center of spectrum
fft_baseband_shifted = fftshift(fft_baseband);

% Magnitude of FFT
magnitude_fft_baseband = abs(fft_baseband_shifted);

% Plot the magnitude of FFT (Frequency Domain)
figure;
plot(f_baseband, magnitude_fft_baseband);
title('baseband Single Tone Signal in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


inv_g=1/g;
rejected_image_signal= I_filtered + 1j*inv_g*Q_filtered;
N_rejected = length(rejected_image_signal);  % Number of points
fft_rejected = fft(rejected_image_signal);  % Compute FFT

% Frequency vector
f_rejected = (-N_rejected/2:N_rejected/2-1)*(fs/N_rejected);  % Frequency range from -fs/2 to fs/2

% Shift zero frequency component to center of spectrum
fft_rejected_shifted = fftshift(fft_rejected);

% Magnitude of FFT
magnitude_fft_rejected = abs(fft_rejected_shifted);

% Plot the magnitude of FFT (Frequency Domain)
figure;
plot(f_rejected, magnitude_fft_rejected);
title('Baseband Single Tone Signal with no image in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

