clear;
tt=edfread('muscle3.edf') 
tt 
fs=1000;  
signalA=tt.ECGBIT; 
nw=max(size(signalA));
for i=1:nw 
    x(i,:)=signalA{i,:}; 
end
signal=reshape(x',1,nw*1000);
save('ECG_12.mat','fs','signal'); 

signaln=(signal-mean(signal))./std(signal);
nsamples=length(signaln);
Ts=1/fs; 
t1=[0:nsamples-1]*Ts;

tt = edfread('refrence3.edf');
tt
fs = 1000;
signalA = tt.EMGBITREV; % Cambiar a EMG
nw = max(size(signalA));
for i = 1:nw
    x(i,:) = signalA{i,:};
end
signal = reshape(x',1,nw*1000);
save('muscle3.mat','fs','signal');

signaln2 = (signal - mean(signal)) ./ std(signal);
nsamples = length(signaln2);
Ts = 1/fs;
t2 = [0:nsamples-1] * Ts;

fc=50
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
signaln = filtfilt(b, a, signaln);
signaln2 = filtfilt(b, a, signaln2);

figure;
subplot(2,2,1)
nsamples = length(signaln);
Ts = 1/fs;
t = (0:nsamples-1) * Ts;
plot(t, signaln+signaln2 (1:nsamples));
xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('ECG of a 30 year-old Woman with Muscle Noise');
subplot(2,2,2)
plot(t, signaln2 (1:nsamples));
xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('Muscle Noise');

f = (0:nsamples-1)*(fs/nsamples);
signaln_fft = abs(fft(signaln+ signaln2(1:nsamples)));
subplot(2,2,3)
plot(f, signaln_fft)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Spectrum of ECG with Muscle Noise')
signaln2_fft = abs(fft(signaln2(1:nsamples)));
subplot(2,2,4)
plot(f, signaln2_fft)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Spectrum of Muscle Noise')

%Muscle 1
ecg_noise = signaln + signaln2(1:nsamples); % ECG signal corrupted by noise
noise = signaln2(1:nsamples); % Noise signal
order = 4; 
lambda = 0.5; 
N = length(ecg_noise);
Rxx = xcorr(ecg_noise, ecg_noise); 
Rx = Rxx(N:N+order);
Rxy = xcorr(noise, ecg_noise);
P = Rxy(N:N+order);
H = P ./ Rx;
G = conj(H) ./ (abs(H).^2 + lambda);
denoised_signal = filter(G, 1, ecg_noise);
corr_coef = corrcoef(denoised_signal, signaln);
coeficiente_correlacion = corr_coef(1, 2);
figure;
subplot(2,1,1);
plot(t, ecg_noise);
title('ECG of a 30 year-old Woman with Muscle Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, signaln);
title('Denoised ECG Signal (Wiener Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');





%ECG1
fc=50
fs=1000
duration=20
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, signaln);
figure
subplot(3,1,1);
plot(t, signaln);
title('Sitting Down ECG Signal of a 55-Year-Old Woman');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(3,1,2);
plot(t, filtered_ecg);
title('Denoised ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
order = 4;
cutoff_freq = 30
nyquist_freq = fs/2;
Wn = cutoff_freq / nyquist_freq;
[b, a] = butter(order, Wn, 'low');
denoised_ecg = filtfilt(b, a, filtered_ecg);
subplot(3,1,3);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth + Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
figure
ecg_fft = abs(fft(signaln));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,1)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 55-Year-Old Woman');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(filtered_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,2)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 55-Year-Old Woman After Notch Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(denoised_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,3)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 55-Year-Old Woman After Notch abd Butterworth Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%ECG2
fc=50
fs=1000
duration=18
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, signaln);
figure
subplot(3,1,1);
plot(t, signaln);
title('Sitting Down ECG Signal of a 55-Year-Old Woman');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(3,1,2);
plot(t, filtered_ecg);
title('Denoised ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
order = 4;
cutoff_freq = 25
nyquist_freq = fs/2;
Wn = cutoff_freq / nyquist_freq;
[b, a] = butter(order, Wn, 'low');
denoised_ecg = filtfilt(b, a, filtered_ecg);
subplot(3,1,3);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth + Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
figure
ecg_fft = abs(fft(signaln));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,1)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 55-Year-Old Woman');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(filtered_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,2)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 55-Year-Old Woman After Notch Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(denoised_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,3)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 55-Year-Old Woman After Notch abd Butterworth Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


%ECG3
fc=50
fs=1000
duration=14
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, signaln);
figure
subplot(3,1,1);
plot(t, signaln);
title('ECG Signal of a 55-Year-Old Woman During Forced Respiration');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(3,1,2);
plot(t, filtered_ecg);
title('Denoised ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
Fc = 2;  
Wc = Fc / (fs / 2);
orden = 4;
[b, a] = ellip(orden, 3, 40, Wc, 'high');
denoised_ecg = filtfilt(b, a, filtered_ecg);
subplot(3,1,3);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Elliptic + Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
figure
ecg_fft = abs(fft(signaln));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,1)
plot(f, ecg_fft);
title('Spectrum of ECG Signal of a 55-Year-Old Woman During Forced Respiration');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(filtered_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,2)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 55-Year-Old Woman During Forced Respiration After Notch Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(denoised_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,3)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 55-Year-Old Woman During Forced Respiration After Notch and Elliptic Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


%ECG4
fc=50
fs=1000
duration=14
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, signaln);
figure
subplot(3,1,1);
plot(t, signaln);
title('Sitting Down ECG Signal of a 50-Year-Old Man');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(3,1,2);
plot(t, filtered_ecg);
title('Denoised ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
order = 4;
cutoff_freq = 30
nyquist_freq = fs/2;
Wn = cutoff_freq / nyquist_freq;
[b, a] = butter(order, Wn, 'low');
denoised_ecg = filtfilt(b, a, filtered_ecg);
subplot(3,1,3);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth + Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
figure
ecg_fft = abs(fft(signaln));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,1)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 50-Year-Old Man');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(filtered_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,2)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 50-Year-Old Man After Notch Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(denoised_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,3)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 50-Year-Old Man After Notch abd Butterworth Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%ECG5
fc=50
fs=1000
duration=15
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, signaln);
figure
subplot(3,1,1);
plot(t, signaln);
title('Laying Down ECG Signal of a 50-Year-Old Man');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(3,1,2);
plot(t, filtered_ecg);
title('Denoised ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
order = 4;
cutoff_freq = 25
nyquist_freq = fs/2;
Wn = cutoff_freq / nyquist_freq;
[b, a] = butter(order, Wn, 'low');
denoised_ecg = filtfilt(b, a, filtered_ecg);
subplot(3,1,3);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth + Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
figure
ecg_fft = abs(fft(signaln));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,1)
plot(f, ecg_fft);
title('Spectrum of Laying Down ECG Signal of a 50-Year-Old Man');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(filtered_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,2)
plot(f, ecg_fft);
title('Spectrum of Laying Down ECG Signal of a 50-Year-Old Man After Notch Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(denoised_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,3)
plot(f, ecg_fft);
title('Spectrum of Laying Down ECG Signal of a 50-Year-Old Man After Notch abd Butterworth Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%ECG6
fc=50
fs=1000
duration=20
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, signaln);
figure
subplot(3,1,1);
plot(t, signaln);
title('ECG Signal of a 50-Year-Old Man During Forced Respiration');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(3,1,2);
plot(t, filtered_ecg);
title('Denoised ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
Fc = 3;  
Wc = Fc / (fs / 2);
orden = 4;
[b, a] = ellip(orden, 3, 40, Wc, 'high');
denoised_ecg = filtfilt(b, a, filtered_ecg);
subplot(3,1,3);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Elliptic + Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
figure
ecg_fft = abs(fft(signaln));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,1)
plot(f, ecg_fft);
title('Spectrum of ECG Signal of a 50-Year-Old Man During Forced Respiration');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(filtered_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,2)
plot(f, ecg_fft);
title('Spectrum of ECG Signal of a 50-Year-Old Man During Forced Respiration After Notch Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(denoised_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,3)
plot(f, ecg_fft);
title('Spectrum of ECG Signal of a 50-Year-Old Man During Forced Respiration After Notch and Elliptic Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%ECG7
fc=50
fs=1000
duration=25
t = 0:1/fs:(duration-1/fs);
bandwidth=4
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, signaln);
figure
subplot(3,1,1);
plot(t, signaln);
title('Sitting Down ECG Signal of a 20-Year-Old Man');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(3,1,2);
plot(t, filtered_ecg);
title('Denoised ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
order = 4;
cutoff_freq = 15
nyquist_freq = fs/2;
Wn = cutoff_freq / nyquist_freq;
[b, a] = butter(order, Wn, 'low');
denoised_ecg = filtfilt(b, a, filtered_ecg);
subplot(3,1,3);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth + Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
figure
ecg_fft = abs(fft(signaln));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,1)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 20-Year-Old Man');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(filtered_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,2)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 20-Year-Old Man After Notch Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(denoised_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,3)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 20-Year-Old Man After Notch and Butterworth Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%ECG8
fc=50
fs=1000
duration=14
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, signaln);
figure
subplot(3,1,1);
plot(t, signaln);
title('ECG Signal of a 70-Year-Old Man During Forced Respiration');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(3,1,2);
plot(t, filtered_ecg);
title('Denoised ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
Fc = 3;  
Wc = Fc / (fs / 2);
orden = 4;
[b, a] = ellip(orden, 3, 40, Wc, 'high');
denoised_ecg = filtfilt(b, a, filtered_ecg);
subplot(3,1,3);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Elliptic + Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
figure
ecg_fft = abs(fft(signaln));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,1)
plot(f, ecg_fft);
title('Spectrum of ECG Signal of a 70-Year-Old Man During Forced Respiration');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(filtered_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,2)
plot(f, ecg_fft);
title('Spectrum of ECG Signal of a 70-Year-Old Man During Forced Respiration After Notch Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(denoised_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,3)
plot(f, ecg_fft);
title('Spectrum of ECG Signal of a 70-Year-Old Man During Forced Respiration After Notch and Elliptic Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%ECG9 
fc=50
fs=1000
duration=20
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, signaln);
figure
subplot(3,1,1);
plot(t, signaln);
title('Sitting Down ECG Signal of a 70-Year-Old Man');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(3,1,2);
plot(t, filtered_ecg);
title('Denoised ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
order = 4;
cutoff_freq = 25
nyquist_freq = fs/2;
Wn = cutoff_freq / nyquist_freq;
[b, a] = butter(order, Wn, 'low');
denoised_ecg = filtfilt(b, a, filtered_ecg);
subplot(3,1,3);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth + Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
figure
ecg_fft = abs(fft(signaln));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,1)
plot(f, ecg_fft);
title('Spectrum of Sitting Down ECG Signal of a 70-Year-Old Man');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(filtered_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,2)
plot(f, ecg_fft);
title('Spectrum After Notch Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ecg_fft = abs(fft(denoised_ecg));
f = linspace(0, fs, length(ecg_fft));
subplot(3,1,3)
plot(f, ecg_fft);
title('Spectrum After Notch and Butterworth Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');