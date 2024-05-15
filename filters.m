


%Notch filter Pline+ RealECG
fc=50
fs=500
duration=10
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, Pline_realECG);
figure
subplot(2,1,1);
plot(t, Pline_realECG);
title('Powerline Interference (50 Hz) + Real ECG');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, filtered_ecg);
title('Filtered ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
filtered_intereference = filtfilt(b, a, interference);
 mean_filtered_intereference=mean(filtered_intereference)
 var_filtered_intereference=var(filtered_intereference)
 power_filtered_intereference=var_filtered_intereference+(mean_filtered_intereference^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_intereference);
num_muestras =0.3* fs
mse = mean((filtered_ecg(num_muestras:(end-num_muestras))-realECGFiltered(num_muestras:(end-num_muestras))).^2)
corr_coef = corrcoef(filtered_ecg, realECGFiltered);
coeficiente_correlacion = corr_coef(1, 2);

%Notch filter Pline+ SynECG
fc=50
fs=500
duration=10
bandwidth=4
t = 0:1/fs:(duration-1/fs);
[b, a] = iirnotch(fc/(fs/2), bandwidth /(fs/2));
filtered_ecg = filtfilt(b, a, Pline_synECG);
figure
subplot(2,1,1);
plot(t, Pline_synECG);
title('Powerline Interference (50 Hz) + Synthetic ECG');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, filtered_ecg);
title('Filtered ECG Signal (Notch Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
filtered_intereference = filtfilt(b, a, interference);
 mean_filtered_intereference=mean(filtered_intereference)
 var_filtered_intereference=var(filtered_intereference)
 power_filtered_intereference=var_filtered_intereference+(mean_filtered_intereference^2)
filtered_synECG = filtfilt(b, a, synECG);
mean_filtered_synECG=mean(filtered_synECG)
 var_filtered_synECG=var(filtered_synECG)
 power_filtered_synECG=var_filtered_synECG+(mean_filtered_synECG^2)
snroutput_dB = 10 * log10(power_filtered_synECG / power_filtered_intereference);
mse = mean((filtered_ecg-synECG).^2)
corr_coef = corrcoef(filtered_ecg, synECG);
coeficiente_correlacion = corr_coef(1, 2);

%Butterworth Muscle Noise EMG + Real ECG
fs = 500;
duration=10
t = 0:1/fs:(duration-1/fs); 
fft_result = fft(muscleNoiseEMG);
frequencies = (0:length(fft_result)-1) * fs/length(fft_result);
figure;
plot(frequencies, abs(fft_result),'r');
hold on;
fft_result2= fft(realECGFiltered);
plot(frequencies, abs(fft_result2),'b');
title('Frequency Spectrum of Real ECG and Muscle Noise');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Muscle Noise', ' Real ECG');
    % --- LOWPASS FILTER 
cutoff_freq = 40;         
filter_order = 4;          
[b, a] = butter(filter_order, cutoff_freq / (fs/2), 'low');
denoised_ecg = filtfilt(b, a, s3);
figure;
subplot(2,1,1);
plot(t, s3);
title('Real ECG with Muscle Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
filtered_muscle= filtfilt(b, a, muscleNoiseEMG);
 mean_filtered_muscle=mean(filtered_muscle)
 var_filtered_muscle=var(filtered_muscle)
 power_filtered_muscle=var_filtered_muscle+(mean_filtered_muscle^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_muscle)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)

    % --- HIGHPASS FILTER 
cutoff_freq = 5;         
filter_order = 4;          
[b, a] = butter(filter_order, cutoff_freq / (fs/2), 'high');
denoised_ecg = filtfilt(b, a, s3);
figure;
subplot(2,1,1);
plot(t, s3);
title('Real ECG with Muscle Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
filtered_muscle= filtfilt(b, a, muscleNoiseEMG);
 mean_filtered_muscle=mean(filtered_muscle)
 var_filtered_muscle=var(filtered_muscle)
 power_filtered_muscle=var_filtered_muscle+(mean_filtered_muscle^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_muscle)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)

%LMS filter Real ECG + Power Intereference
L = 2;
lms = dsp.LMSFilter(L,'Method','LMS');
[mumaxlms,mumaxmselms] = maxstep(lms, Pline_realECG);
lms.StepSize = mumaxlms/2;
ref_signal = power_line_amplitude* cos(2 * pi * 50 * t);
[~, elms, wlms] = lms(ref_signal', Pline_realECG');
figure;
subplot(2,1,1);
plot(t, Pline_realECG);
title('Real ECG with Power Line Interference');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, elms);
title('Denoised ECG Signal (Adaptive Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (elms- mean(elms))/std(elms)
mse = mean((denoised_ecg_normalize - realECGFiltered_normalize').^2)
corr_coef = corrcoef(elms, realECGFiltered);
coeficiente_correlacion = corr_coef(1, 2);

%LMS filter Syn ECG + Power Intereference
L = 10;
lms = dsp.LMSFilter(L,'Method','LMS');
[mumaxlms,mumaxmselms] = maxstep(lms, Pline_synECG);
lms.StepSize = mumaxlms/2;
ref_signal = 1* cos(2 * pi * 50 * t);
[~, elms, wlms] = lms(ref_signal', Pline_synECG');
figure;
subplot(2,1,1);
plot(t, Pline_synECG);
title('Synthetic ECG with Power Line Interference (50 Hz)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, elms);
title('Denoised ECG Signal (Adaptive Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
synECGFiltered_normalize= (synECG-mean_synECG)/std(synECG)
denoised_ecg_normalize= (elms- mean(elms))/std(elms)
mse = mean((denoised_ecg_normalize - synECGFiltered_normalize').^2)
corr_coef = corrcoef(elms, synECG);
coeficiente_correlacion = corr_coef(1, 2);

%Powerline Interference + Real ECG (Butterowrth)
fs = 500
power_line_freq = 50
ecg_band = [0.7 50]; 
[b, a] = butter(4, [(power_line_freq-0.5) (power_line_freq+0.5)]/(fs/2), 'stop');
filtered_signal = filtfilt(b, a, Pline_realECG);
[b_ecg, a_ecg] = butter(4, ecg_band/(fs/2), 'bandpass');
filtered_ecg_signal = filtfilt(b_ecg, a_ecg, filtered_signal);
figure;
subplot(2,1,1);
plot(t, Pline_realECG);
title('Real ECG with Power Line Interference');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, filtered_ecg_signal);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
num_muestras =1* fs
filtered_intereference = filtfilt(b_ecg, a_ecg, interference);
 var_filtered_intereference=var(filtered_intereference)
 power_filtered_intereference=var_filtered_intereference+(mean_filtered_intereference^2)
filtered_realECG = filtfilt(b_ecg, a_ecg, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_intereference);
 mean_filtered_intereference=mean(filtered_intereference)
 realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (filtered_intereference- mean(filtered_intereference))/std(filtered_intereference)
mse = mean((realECGFiltered_normalize (num_muestras:(end-num_muestras)) - denoised_ecg_normalize(num_muestras:(end-num_muestras))).^2)
corr_coef = corrcoef(filtered_ecg_signal, realECGFiltered);
coeficiente_correlacion = corr_coef(1, 2);

%Powerline Interference + Syn ECG (Butterowrth)
fs = 500
power_line_freq = 50
ecg_band = [0.7 40]; 
[b, a] = butter(4, [(power_line_freq-0.5) (power_line_freq+0.5)]/(fs/2), 'stop');
filtered_signal = filtfilt(b, a, Pline_synECG);
[b_ecg, a_ecg] = butter(4, ecg_band/(fs/2), 'bandpass');
filtered_ecg_signal = filtfilt(b_ecg, a_ecg, filtered_signal);
figure;
subplot(2,1,1);
plot(t, Pline_synECG);
title('Synthetic ECG with Power Line Interference (50 Hz)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, filtered_ecg_signal);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
filtered_intereference = filtfilt(b_ecg, a_ecg, interference);
 var_filtered_intereference=var(filtered_intereference)
 power_filtered_intereference=var_filtered_intereference+(mean_filtered_intereference^2)
filtered_synECG = filtfilt(b_ecg, a_ecg, synECG);
mean_filtered_synECG=mean(filtered_synECG)
 var_filtered_synECG=var(filtered_synECG)
 power_filtered_synECG=var_filtered_synECG+(mean_filtered_synECG^2)
snroutput_dB = 10 * log10(power_filtered_synECG / power_filtered_intereference);
 mean_filtered_intereference=mean(filtered_intereference)
 synECGFiltered_normalize= (synECG-mean_filtered_synECG)/std(synECG)
num_muestras=1*fs;
 denoised_ecg_normalize= (filtered_intereference- mean(filtered_intereference))/std(filtered_intereference)
mse = mean((synECGFiltered_normalize (num_muestras:(end-num_muestras)) - denoised_ecg_normalize(num_muestras:(end-num_muestras))).^2)
corr_coef = corrcoef(filtered_ecg_signal, synECG);
coeficiente_correlacion = corr_coef(1, 2);


%Real ECG+ Baseline Wander Butterworth
fs = 500;
duration=10
t = 0:1/fs:(duration-1/fs); 
fft_result = fft(BaseLineNoise);
fft_result2 = fft(realECGFiltered);
frequencies = (0:length(fft_result)-1) * fs/length(fft_result);
figure;
plot(frequencies, abs(fft_result),'r');
hold on;
plot(frequencies, abs(fft_result2),'b');
title('Frequency Spectrum of Baseline Noise and Real ECG');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Baseline Noise', ' Real ECG');
fs=500
cutoff_frequency = 4; 
filter_order =8; 
[b, a] = butter(filter_order, cutoff_frequency / (fs/2), 'high');
denoised_ecg = filtfilt(b, a, BaseLine_RealECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_RealECG);
title('Real ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
filtered_baseline= filtfilt(b, a, BaseLineNoise);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_baseline)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, realECGFiltered);
coeficiente_correlacion = corr_coef(1, 2);

%Syn ECG+ Baseline Wander Butterworth
fs = 500;
duration=10
t = 0:1/fs:(duration-1/fs); 
fft_result = fft(0.005*BaseLineNoise);
fft_result2 = fft(synECG);
frequencies = (0:length(fft_result)-1) * fs/length(fft_result);
figure;
plot(frequencies, abs(fft_result),'r');
hold on;
plot(frequencies, abs(fft_result2),'b');
title('Frequency Spectrum of Baseline Noise and Synthetic ECG');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Baseline Noise', ' Synthetic ECG');
fs=500
cutoff_frequency = 4; 
filter_order =7; 
[b, a] = butter(filter_order, cutoff_frequency / (fs/2), 'high');
denoised_ecg = filtfilt(b, a, BaseLine_synECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_synECG);
title('Synthetic ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
filtered_baseline= filtfilt(b, a, 0.005*BaseLineNoise);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_synECG = filtfilt(b, a, synECG);
mean_filtered_synECG=mean(filtered_synECG)
 var_filtered_synECG=var(filtered_synECG)
 power_filtered_synECG=var_filtered_synECG+(mean_filtered_synECG^2)
snroutput_dB = 10 * log10(power_filtered_synECG / power_filtered_baseline)
synECGFiltered_normalize= (synECG-mean_filtered_synECG)/std(synECG)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((synECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, synECG);
coeficiente_correlacion = corr_coef(1, 2);

%Real ECG+ Baseline Wander II Butterworth
fs = 500;
duration=10
t = 0:1/fs:(duration-1/fs); 
fft_result = fft(BaseLineNoise_extracted);
fft_result2 = fft(realECGFiltered);
frequencies = (0:length(fft_result)-1) * fs/length(fft_result);
figure;
plot(frequencies, abs(fft_result),'r');
hold on;
plot(frequencies, abs(fft_result2),'b');
title('Frequency Spectrum of Baseline Noise and Real ECG');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Baseline Noise', ' Real ECG');
fs=500
cutoff_frequency = 4; 
filter_order =8; 
[b, a] = butter(filter_order, cutoff_frequency / (fs/2), 'high');
BaseLine_RealECG = realECGFiltered + BaseLineNoise_extracted
figure;
plot(t, BaseLine_RealECG);
title('Real ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
denoised_ecg = filtfilt(b, a, BaseLine_RealECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_RealECG);
title('Real ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
grid on;
filtered_baseline= filtfilt(b, a, BaseLineNoise_extracted);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_baseline)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, realECGFiltered);
coeficiente_correlacion = corr_coef(1, 2);

%Syn ECG+ Baseline Wander II Butterworth
fs = 500;
duration=10
t = 0:1/fs:(duration-1/fs); 
fft_result = fft(0.005*BaseLineNoise_extracted);
fft_result2 = fft(synECG);
frequencies = (0:length(fft_result)-1) * fs/length(fft_result);
figure;
plot(frequencies, abs(fft_result),'r');
hold on;
plot(frequencies, abs(fft_result2),'b');
title('Frequency Spectrum of Baseline Noise and Synthetic ECG');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Baseline Noise', ' Synthetic ECG');
fs=500
cutoff_frequency =3; 
filter_order =6; 
[b, a] = butter(filter_order, cutoff_frequency / (fs/2), 'high');
BaseLine_RealECG = synECG + 0.005* BaseLineNoise_extracted
figure;
plot(t, BaseLine_RealECG);
title('Synthetic ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
denoised_ecg = filtfilt(b, a, BaseLine_RealECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_RealECG);
title('Synthetic ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
grid on;
filtered_baseline= filtfilt(b, a, 0.005* BaseLineNoise_extracted);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_realECG = filtfilt(b, a, synECG);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_baseline)
realECGFiltered_normalize= (synECG-mean_filtered_realECG)/std(synECG)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, synECG);
coeficiente_correlacion = corr_coef(1, 2);

%Real ECG + Baseline Wander (Elliptic)
fs = 500;  
Rp = 0.1;  
Rs = 60;    
Fp = 1.5;     
Fs = 3;     
Fp_normalized = Fp / (fs/2);
Fs_normalized = Fs / (fs/2);
[n, Wn] = ellipord(Fp_normalized, Fs_normalized, Rp, Rs);
[b, a] = ellip(n, Rp, Rs, Wn,'high');
denoised_ecg = filtfilt(b, a, BaseLine_RealECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_RealECG);
title('Real ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Elliptic Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
filtered_baseline= filtfilt(b, a, BaseLineNoise);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_baseline)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, realECGFiltered);
coeficiente_correlacion = corr_coef(1, 2);

%Syn ECG + Baseline Wander (Elliptic)
fs = 500;  
Rp = 0.1;  
Rs = 60;    
Fp = 2;     
Fs = 3;     
Fp_normalized = Fp / (fs/2);
Fs_normalized = Fs / (fs/2);
[n, Wn] = ellipord(Fp_normalized, Fs_normalized, Rp, Rs);
[b, a] = ellip(n, Rp, Rs, Wn,'high');
denoised_ecg = filtfilt(b, a, BaseLine_synECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_synECG);
title('Synthetic ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Elliptic Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
filtered_baseline= filtfilt(b, a, 0.005 *BaseLineNoise);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_synECG = filtfilt(b, a, synECG);
mean_filtered_synECG=mean(filtered_synECG)
 var_filtered_synECG=var(filtered_synECG)
 power_filtered_synECG=var_filtered_synECG+(mean_filtered_synECG^2)
snroutput_dB = 10 * log10(power_filtered_synECG / power_filtered_baseline)
synECGFiltered_normalize= (synECG-mean_filtered_synECG)/std(synECG)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((synECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, synECG);
coeficiente_correlacion = corr_coef(1, 2);

%Real ECG + Baseline Wander II (Elliptic)
fs = 500;  
Rp = 0.1;  
Rs = 60;    
Fp = 1.5;     
Fs = 3;     
Fp_normalized = Fp / (fs/2);
Fs_normalized = Fs / (fs/2);
[n, Wn] = ellipord(Fp_normalized, Fs_normalized, Rp, Rs);
[b, a] = ellip(n, Rp, Rs, Wn,'high');
denoised_ecg = filtfilt(b, a, BaseLine_RealECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_RealECG);
title('Real ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Elliptic Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
filtered_baseline= filtfilt(b, a, BaseLineNoise_extracted);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_baseline)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, realECGFiltered);
coeficiente_correlacion = corr_coef(1, 2);

%Syn ECG + Baseline Wander II (Elliptic)
fs = 500;  
Rp = 0.1;  
Rs = 60;    
Fp = 1.5;     
Fs = 3;     
Fp_normalized = Fp / (fs/2);
Fs_normalized = Fs / (fs/2);
[n, Wn] = ellipord(Fp_normalized, Fs_normalized, Rp, Rs);
[b, a] = ellip(n, Rp, Rs, Wn,'high');
denoised_ecg = filtfilt(b, a, BaseLine_RealECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_RealECG);
title('Synthetic ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Elliptic Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
filtered_baseline= filtfilt(b, a, 0.005*BaseLineNoise_extracted);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_realECG = filtfilt(b, a, synECG);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_baseline)
realECGFiltered_normalize= (synECG-mean_filtered_realECG)/std(synECG)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, synECG);
coeficiente_correlacion = corr_coef(1, 2);

%Baseline Wander + Real ECG (Cheb. 2)
fs = 500;   
Rp = 0.1;  
Rs = 50;    
Fp = 0.5/(fs/2);   
Fs = 1/(fs/2);   
[n, Wn] = cheb2ord(Fp, Fs, Rp, Rs, 's');
[b, a] = cheby2(n, Rs, Wn, 'high');
denoised_ecg = filtfilt(b, a, BaseLine_RealECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_RealECG);
title('Real ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Chebyshev II Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
filtered_baseline= filtfilt(b, a, BaseLineNoise);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_baseline)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, realECGFiltered);
coeficiente_correlacion = corr_coef(1, 2);

%Baseline Wander + Syn ECG (Cheb. 2)
fs = 500;   
Rp = 0.1;  
Rs = 50;    
Fp = 0.5/(fs/2);   
Fs = 1/(fs/2);   
[n, Wn] = cheb2ord(Fp, Fs, Rp, Rs, 's');
[b, a] = cheby2(n, Rs, Wn, 'high');
denoised_ecg = filtfilt(b, a, BaseLine_synECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_synECG);
title('Synthetic ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Chebyshev II Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
filtered_baseline= filtfilt(b, a, 0.005*BaseLineNoise);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_synECG = filtfilt(b, a, synECG);
mean_filtered_synECG=mean(filtered_synECG)
 var_filtered_synECG=var(filtered_synECG)
 power_filtered_synECG=var_filtered_synECG+(mean_filtered_synECG^2)
snroutput_dB = 10 * log10(power_filtered_synECG / power_filtered_baseline)
synECGFiltered_normalize= (realECGFiltered-mean_filtered_synECG)/std(synECG)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((synECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, synECG);
coeficiente_correlacion = corr_coef(1, 2);

%Baseline Wander II + Real ECG (Cheb. 2)
fs = 500;   
Rp = 0.1;  
Rs = 50;    
Fp = 0.5/(fs/2);   
Fs = 1/(fs/2);   
[n, Wn] = cheb2ord(Fp, Fs, Rp, Rs, 's');
[b, a] = cheby2(n, Rs, Wn, 'high');
denoised_ecg = filtfilt(b, a, BaseLine_RealECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_RealECG);
title('Real ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Chebyshev II Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
filtered_baseline= filtfilt(b, a, BaseLineNoise_extracted);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_baseline)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, realECGFiltered);
coeficiente_correlacion = corr_coef(1, 2);

%Baseline Wander II + Syn ECG (Cheb. 2)
fs = 500;   
Rp = 0.1;  
Rs = 50;    
Fp = 0.5/(fs/2);   
Fs = 1/(fs/2);   
[n, Wn] = cheb2ord(Fp, Fs, Rp, Rs, 's');
[b, a] = cheby2(n, Rs, Wn, 'high');
denoised_ecg = filtfilt(b, a, BaseLine_RealECG);
figure;
subplot(2,1,1);
plot(t, BaseLine_RealECG);
title('Synthetic ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Chebyshev II Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
filtered_baseline= filtfilt(b, a, 0.005*BaseLineNoise_extracted);
 mean_filtered_baseline=mean(filtered_baseline)
 var_filtered_baseline=var(filtered_baseline)
 power_filtered_baseline=var_filtered_baseline+(mean_filtered_baseline^2)
filtered_realECG = filtfilt(b, a, synECG);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_baseline)
realECGFiltered_normalize= (synECG-mean_filtered_realECG)/std(synECG)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, synECG);
coeficiente_correlacion = corr_coef(1, 2);

%LMS filter Real ECG + Baseline Wander
L = 2;
lms = dsp.LMSFilter(L,'Method','LMS');
[mumaxlms,mumaxmselms] = maxstep(lms, BaseLine_RealECG);
lms.StepSize = mumaxlms/2;
ref_signal = BaseLineNoise_extracted
[~, elms, wlms] = lms(ref_signal', BaseLine_RealECG');
figure;
subplot(2,1,1);
plot(t, BaseLine_RealECG);
title('Real ECG with Baseline Wander Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, elms);
title('Denoised ECG Signal (Adaptive Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (elms- mean(elms))/std(elms)
mse = mean((denoised_ecg_normalize - realECGFiltered_normalize').^2)
corr_coef = corrcoef(BaseLineNoise, ref_signal)
coeficiente_correlacion = corr_coef(1, 2)

%Butterworth Muscle Noise + Real ECG
fs = 500;
duration=10
t = 0:1/fs:(duration-1/fs); 
fft_result = fft(muscleNoise);
frequencies = (0:length(fft_result)-1) * fs/length(fft_result);
figure;
plot(frequencies, abs(fft_result),'r');
hold on;
fft_result2= fft(realECGFiltered);
plot(frequencies, abs(fft_result2),'b');
title('Frequency Spectrum of Real ECG and Muscle Noise');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Muscle Noise', ' Real ECG');
    % --- LOWPASS FILTER 
cutoff_freq = 50;         
filter_order = 4;          
[b, a] = butter(filter_order, cutoff_freq / (fs/2), 'low');
denoised_ecg = filtfilt(b, a, s4);
figure;
subplot(2,1,1);
plot(t, s4);
title('Real ECG with Muscle Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
filtered_muscle= filtfilt(b, a, muscleNoise);
 mean_filtered_muscle=mean(filtered_muscle)
 var_filtered_muscle=var(filtered_muscle)
 power_filtered_muscle=var_filtered_muscle+(mean_filtered_muscle^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_muscle)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)

    % --- HIGHPASS FILTER 
cutoff_freq = 5;         
filter_order = 4;          
[b, a] = butter(filter_order, cutoff_freq / (fs/2), 'high');
denoised_ecg = filtfilt(b, a, s4);
figure;
subplot(2,1,1);
plot(t, s4);
title('Real ECG with Muscle Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
filtered_muscle= filtfilt(b, a, muscleNoise);
 mean_filtered_muscle=mean(filtered_muscle)
 var_filtered_muscle=var(filtered_muscle)
 power_filtered_muscle=var_filtered_muscle+(mean_filtered_muscle^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_muscle)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)


%Chev. II Muscle Noise + Real ECG
   % --- LOWPASS
fs = 500;   
Rp = 1;  
Rs = 50;    
Fp = 50/(fs/2);   
Fs = 150/(fs/2);   
[n, Wn] = cheb2ord(Fp, Fs, Rp, Rs, 's');
[b, a] = cheby2(n, Rs, Wn, 'high');
denoised_ecg = filtfilt(b, a, s4);
figure;
subplot(2,1,1);
plot(t, s4);
title('Real ECG with Muscle Noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Chebyshev II Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
filtered_muscle= filtfilt(b, a, muscleNoise);
 mean_filtered_muscle=mean(filtered_muscle)
 var_filtered_muscle=var(filtered_muscle)
 power_filtered_muscle=var_filtered_muscle+(mean_filtered_muscle^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_muscle)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)

%Butterworth Real ECG + Gaussian Noise
fs=500
duration=10
t = 0:1/fs:(duration-1/fs);
order = 4;
cutoff_freq = 60
nyquist_freq = fs/2;
Wn = cutoff_freq / nyquist_freq;
[b, a] = butter(order, Wn, 'low');
denoised_ecg = filtfilt(b, a, AWGN_Realsignal{7});
figure;
subplot(2,1,1);
plot(t, AWGN_Realsignal{7});
title('Real ECG with AWGN, SNR: 10 dB');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t, denoised_ecg);
title('Denoised ECG Signal (Butterworth Filter)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
grid on;
filtered_AWGN= filtfilt(b, a, white_noise{7});
 mean_filtered_AWGN=mean(filtered_AWGN)
 var_filtered_AWGN=var(filtered_AWGN)
 power_filtered_AWGN=var_filtered_AWGN+(mean_filtered_AWGN^2)
filtered_realECG = filtfilt(b, a, realECGFiltered);
mean_filtered_realECG=mean(filtered_realECG)
 var_filtered_realECG=var(filtered_realECG)
 power_filtered_realECG=var_filtered_realECG+(mean_filtered_realECG^2)
 snroutput_dB = 10 * log10(power_filtered_realECG /power_filtered_AWGN)
realECGFiltered_normalize= (realECGFiltered-mean_filtered_realECG)/std(realECGFiltered)
denoised_ecg_normalize= (denoised_ecg- mean(denoised_ecg))/std(denoised_ecg)
mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2)
corr_coef = corrcoef(denoised_ecg, realECGFiltered);
coeficiente_correlacion = corr_coef(1, 2);
% Buttwerworth AWGN + SynECG
fs = 500;
duration = 10;
t = 0:1/fs:(duration-1/fs);
order = 4;
cutoff_freq = 50;
nyquist_freq = fs/2;
snroutput_dB= [];
mse = [];
coeficiente_correlacion = [];
SNR=[-20:5:30]
for i = 1:length(SNR)

    Wn = cutoff_freq / nyquist_freq;
    [b, a] = butter(order, Wn, 'low');
    denoised_ecg = filtfilt(b, a, AWGN_SYNsignal{i});
    filtered_noise = filtfilt(b, a, white_noise{i});
    mean_filtered_noise = mean(filtered_noise);
    var_filtered_noise = var(filtered_noise);
    power_filtered_noise = var_filtered_noise + (mean_filtered_noise^2);
    filtered_ecg = filtfilt(b, a, synECG);
    mean_filtered_ecg = mean(filtered_ecg);
    var_filtered_ecg = var(filtered_ecg);
    power_filtered_ecg = var_filtered_ecg + (mean_filtered_ecg^2);
    snroutput_dB (i)= 10 * log10(power_filtered_ecg / power_filtered_noise);
    
    realECGFiltered_normalize = (synECG - mean_synECG) / std(synECG);
    denoised_ecg_normalize = (denoised_ecg - mean(denoised_ecg)) / std(denoised_ecg);
    mse (i) = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2);
    
    corr_coef = corrcoef(denoised_ecg, synECG);
    coeficiente_correlacion (i) = corr_coef(1, 2);
    figure;
    subplot(2,1,1);
    plot(t, AWGN_SYNsignal{i});
    title(['Synthetic ECG with AWGN, SNR: ' num2str(SNR(i)) ' dB']);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;

    subplot(2,1,2);
    plot(t, denoised_ecg);
    title('Denoised ECG Signal (Butterworth Filter)');
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;
end

%Chev. II Real ECG + Gaussian Noise
fs = 500;
duration = 10;
t = 0:1/fs:(duration-1/fs);
Rp = 1;
Rs = 50;
Fp = 50/(fs/2);
Fs = 60/(fs/2);
snroutput_dB_all = [];
mse_all = [];
coeficiente_correlacion_all = [];
for snr_index = 1:length(SNR)
    [n, Wn] = cheb2ord(Fp, Fs, Rp, Rs, 's');
    [b, a] = cheby2(n, Rs, Wn, 'low');
    denoised_ecg = filtfilt(b, a, AWGN_Realsignal{snr_index});
    filtered_AWGN = filtfilt(b, a, white_noise{snr_index});
    mean_filtered_AWGN = mean(filtered_AWGN);
    var_filtered_AWGN = var(filtered_AWGN);
    power_filtered_AWGN = var_filtered_AWGN + (mean_filtered_AWGN^2);

    filtered_realECG = filtfilt(b, a, realECGFiltered);
    mean_filtered_realECG = mean(filtered_realECG);
    var_filtered_realECG = var(filtered_realECG);
    power_filtered_realECG = var_filtered_realECG + (mean_filtered_realECG^2);

    snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_AWGN);
    snroutput_dB_all(snr_index) = snroutput_dB;

    realECGFiltered_normalize = (realECGFiltered - mean_filtered_realECG) / std(realECGFiltered);
    denoised_ecg_normalize = (denoised_ecg - mean(denoised_ecg)) / std(denoised_ecg);
    mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2);
    mse_all(snr_index) = mse;

    corr_coef = corrcoef(denoised_ecg, realECGFiltered);
    coeficiente_correlacion = corr_coef(1, 2);
    coeficiente_correlacion_all(snr_index) = coeficiente_correlacion;
end

%Chev. II Syn ECG + Gaussian Noise
fs = 500;
duration = 10;
t = 0:1/fs:(duration-1/fs);
Rp = 1;
Rs = 50;
Fp = 20/(fs/2);
Fs = 30/(fs/2);
snroutput_dB_all = [];
mse_all = [];
coeficiente_correlacion_all = [];

for snr_index = 1:length(SNR)
    [n, Wn] = cheb2ord(Fp, Fs, Rp, Rs, 's');
    [b, a] = cheby2(n, Rs, Wn, 'low');
    
    % Aplicar el filtro de Chebyshev II a la señal AWGN correspondiente
    denoised_ecg = filtfilt(b, a, AWGN_SYNsignal{snr_index});
    
    % Calcular SNR
    filtered_AWGN = filtfilt(b, a, white_noise{snr_index});
    mean_filtered_AWGN = mean(filtered_AWGN);
    var_filtered_AWGN = var(filtered_AWGN);
    power_filtered_AWGN = var_filtered_AWGN + (mean_filtered_AWGN^2);
    
    filtered_realECG = filtfilt(b, a, synECG);
    mean_filtered_realECG = mean(filtered_realECG);
    var_filtered_realECG = var(filtered_realECG);
    power_filtered_realECG = var_filtered_realECG + (mean_filtered_realECG^2);
    
    snroutput_dB = 10 * log10(power_filtered_realECG / power_filtered_AWGN);
    snroutput_dB_all(snr_index) = snroutput_dB;

    % Normalizar las señales para calcular el MSE
    realECGFiltered_normalize = (synECG - mean(filtered_realECG)) / std(filtered_realECG);
    denoised_ecg_normalize = (denoised_ecg - mean(denoised_ecg)) / std(denoised_ecg);
    
    % Calcular el MSE
    mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2);
    mse_all(snr_index) = mse;

    % Calcular el coeficiente de correlación
    corr_coef = corrcoef(denoised_ecg, synECG);
    coeficiente_correlacion = corr_coef(1, 2);
    coeficiente_correlacion_all(snr_index) = coeficiente_correlacion;

    % Visualización de las señales denoised
    figure;
    subplot(2,1,1);
    plot(t, AWGN_SYNsignal{snr_index});
    title(['Synthetic ECG with AWGN, SNR: ' num2str(SNR(snr_index)) ' dB']);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;
    
    subplot(2,1,2);
    plot(t, denoised_ecg);
    title('Denoised ECG (Chebyshev II)');
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;
end




%Wiener Filter AWGN + Real ECG 
order = 4; 
lambda = 0.5; 
denoised_signals_wiener = cell(size(SNR));
coeficientes_correlacion = zeros(size(SNR));

for i = 1:length(SNR)
    N = length(AWGN_Realsignal{i});
    Rxx = xcorr(AWGN_Realsignal{i}, AWGN_Realsignal{i}); 
    Rx = Rxx(N:N+order);
    Rxy = xcorr(AWGN_Realsignal{i}, realECGFiltered);
    P = Rxy(N:N+order);
    H = P ./ Rx;
    G = conj(H) ./ (abs(H).^2 + lambda);
    denoised_signals_wiener{i} = filter(G, 1, AWGN_Realsignal{i});
    corr_coef = corrcoef(denoised_signals_wiener{i}, realECGFiltered);
    coeficiente_correlacion = corr_coef(1, 2);
    coeficientes_correlacion(i) = coeficiente_correlacion;
end
for i = 1:length(SNR)
    figure;
    subplot(2,1,1);
    plot(t, AWGN_Realsignal{i});
    title(['Real ECG with AWGN, SNR: ' num2str(SNR(i)) ' dB']);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;

    subplot(2,1,2);
    plot(t, denoised_signals_wiener{i});
    title(['Denoised ECG Signal (Wiener Filter), SNR: ' num2str(SNR(i)) ' dB']);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;
end

%Wiener Filter AWGN + Syn ECG 
order = 4; 
lambda = 0.5; 
denoised_signals_wiener = cell(size(SNR));
coeficientes_correlacion = zeros(size(SNR));

for i = 1:length(SNR)
    N = length(AWGN_SYNsignal{i});
    Rxx = xcorr(AWGN_SYNsignal{i}, AWGN_SYNsignal{i}); 
    Rx = Rxx(N:N+order);
    Rxy = xcorr(AWGN_SYNsignal{i}, synECG);
    P = Rxy(N:N+order);
    H = P ./ Rx;
    G = conj(H) ./ (abs(H).^2 + lambda);
    denoised_signals_wiener{i} = filter(G, 1, AWGN_SYNsignal{i});
    corr_coef = corrcoef(denoised_signals_wiener{i}, synECG);
    coeficiente_correlacion = corr_coef(1, 2);
    coeficientes_correlacion(i) = coeficiente_correlacion;
end
for i = 1:length(SNR)
    figure;
    subplot(2,1,1);
    plot(t, AWGN_SYNsignal{i});
    title(['Synthetic ECG with AWGN, SNR: ' num2str(SNR(i)) ' dB']);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;

    subplot(2,1,2);
    plot(t, denoised_signals_wiener{i});
    title(['Denoised ECG Signal (Wiener Filter)']);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;
end

%LMS Real ECG + AWGN 
L = 50;
lms = dsp.LMSFilter(L,'Method','LMS');
mse_all = [];
coeficientes_correlacion = [];

 realECGFiltered_normalize = (realECGFiltered - mean(realECGFiltered)) / std(realECGFiltered);
for i = 1:length(SNR)
    [mumaxlms,~] = maxstep(lms, AWGN_Realsignal{i});
    lms.StepSize = mumaxlms / 2;
    ref_signal = white_noise{i};
    [~, elms, ~] = lms(ref_signal', AWGN_Realsignal{i}');
    
    corr_coef = corrcoef(elms, realECGFiltered);
    coeficientes_correlacion(i) = corr_coef(1, 2);
    denoised_ecg_normalize = (elms - mean(elms)) / std(elms);
    
    mse = mean((realECGFiltered_normalize - denoised_ecg_normalize).^2);
    mse_all(i) = mse;
end

%LMS Syn ECG + AWGN 
L = 60;
lms = dsp.LMSFilter(L,'Method','LMS');
mse = [];
coeficientes_correlacion = [];

for i = 1:length(SNR)
    [mumaxlms,~] = maxstep(lms, AWGN_SYNsignal{i});
    lms.StepSize = mumaxlms / 2;
    ref_signal = white_noise{i};
    [~, elms, ~] = lms(ref_signal', AWGN_SYNsignal{i}');
    
    corr_coef = corrcoef(elms, synECG);
    coeficientes_correlacion(i) = corr_coef(1, 2);

    figure;
    subplot(2,1,1);
    plot(t, AWGN_SYNsignal{i});
    title(['Synthetic ECG with AWGN, SNR: ' num2str(SNR(i)) ' dB']);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;

    subplot(2,1,2);
    plot(t, elms);
    title(['Denoised ECG Signal (LMS Filter)']);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;
end
