 run("ECGsyn.m")
 synECG=ecg;
 mean_synECG=mean(synECG)
 var_synECG=var(synECG)
 power_synECG=var_synECG+(mean_synECG^2)

load('ECGFiltered.mat')
realECGFiltered=val;
fsRealECG=500;
durationRealECG=10;
trealECG = 0:1/fsRealECG:(durationRealECG-1/fsRealECG);
plot(trealECG, realECGFiltered);
title('ECG Without Noise');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
 mean_realECGFiltered=mean(realECGFiltered)
 var_realECGFiltered=var(realECGFiltered)
 power_realECGFiltered=var_realECGFiltered+(mean_realECGFiltered^2)

signaln=(realECGFiltered-mean(realECGFiltered))./std(realECGFiltered);
nsamples=length(signaln);
Ts=1/fsRealECG; 
t1=[0:nsamples-1]*Ts;
plot(trealECG, signaln);
title('Filtered ECG from PhysioNet' );
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
 mean_realECGFiltered=mean(signaln)
 var_realECGFiltered=var(signaln)
 power_realECGFiltered=var_realECGFiltered+(mean_realECGFiltered^2)
fft1= fft(signaln);
fft1= fftshift(fft1)
N = length(signaln);
frequencies = (-N/2:N/2-1) * fsRealECG/N;
plot(frequencies, abs(fft1)/N,'b');
xlabel('Frequency (Hz)');
ylabel('Normalized magnitude');

% AWGN (for synECG)
SNR=[-20:5:30]
durationAWGN = 10; 
fsAWGN= 500;
t = 0:1/fsAWGN:(durationAWGN-1/fsAWGN);
powerAWGN = zeros(size(SNR));
amplitudeAWGN = zeros(size(SNR));
white_noise = cell(size(SNR));
AWGN_SYNsignal = cell(size(SNR));
for i=1: length(SNR)
powerAWGN(i)= power_synECG / (10^(SNR(i)/ 10));
amplitudeAWGN(i)=sqrt(powerAWGN(i))
white_noise{i} = amplitudeAWGN(i)* randn(1, durationAWGN * fsAWGN);
figure;
    plot(t, white_noise{i});
    title(['AWGN - SNR: ' num2str(SNR(i)) ' dB']);
    xlabel('Time (seconds)');
    ylabel('Amplitude (mV)');
    grid on;
end
for i=1: length(SNR)
AWGN_SYNsignal{i}=white_noise{i}+ synECG
figure;
    plot(t, AWGN_SYNsignal{i});
    title(['Synthetic ECG + AWGN (SNR: ' num2str(SNR(i)) ' dB)']);
    xlabel('Time (seconds)');
    ylabel('Amplitude (mV)');
    grid on;
end

% AWGN (for realECGFiltered)
SNR=[-20:5:30]
durationAWGN = 10; 
fsAWGN= 500;
t = 0:1/fsAWGN:(durationAWGN-1/fsAWGN);
powerAWGN = zeros(size(SNR));
amplitudeAWGN = zeros(size(SNR));
white_noise = cell(size(SNR));
AWGN_Realsignal = cell(size(SNR));
for i=1: length(SNR)
powerAWGN(i)= power_realECGFiltered / (10^(SNR(i)/ 10));
amplitudeAWGN(i)=sqrt(powerAWGN(i))
white_noise{i} = amplitudeAWGN(i)* randn(1, durationAWGN * fsAWGN);
figure;
    plot(t, white_noise{i});
    title(['AWGN - SNR: ' num2str(SNR(i)) ' dB']);
    xlabel('Time (seconds)');
    ylabel('Amplitude (mV)');
    grid on;
end
for i=1: length(SNR)
AWGN_Realsignal{i}=white_noise{i}+ signaln
figure;
    plot(t, AWGN_Realsignal{i});
    title(['Real ECG + AWGN (SNR: ' num2str(SNR(i)) ' dB)']);
    xlabel('Time (seconds)');
    ylabel('Amplitude (mV)');
    grid on;
end

fft1= fft(white_noise{6});
fft1= fftshift(fft1)
fft2= fft(signaln);
fft2= fftshift(fft2)
N1 = length(white_noise{6});
N2= length(signaln);
frequencies = (-N/2:N/2-1) * fsRealECG/N;
plot(frequencies, abs(fft1)/N1,'r');
hold on;
plot(frequencies, abs(fft2)/N2,'b');
xlabel('Frequency (Hz)');
ylabel('Normalized magnitude');
title ('Frequency Spectrum of AWGN (SNR:5 dB) + PhysioNet´s ECG')
legend ('AWGN (SNR:5 dB)', 'PhysioNet´s ECG')

% AWGN (for synECGFiltered)
SNR=[-20:5:10]
durationAWGN = 10; 
fsAWGN= 500;
t = 0:1/fsAWGN:(durationAWGN-1/fsAWGN);
powerAWGN = zeros(size(SNR));
amplitudeAWGN = zeros(size(SNR));
white_noise = cell(size(SNR));
AWGN_synsignal = cell(size(SNR));
for i=1: length(SNR)
powerAWGN(i)= power_synECG / (10^(SNR(i)/ 10));
amplitudeAWGN(i)=sqrt(powerAWGN(i))
white_noise{i} = amplitudeAWGN(i)* randn(1, durationAWGN * fsAWGN);
end
for i=1: length(SNR)
AWGN_synsignal{i}=white_noise{i}+ synECG
figure;
    plot(t, AWGN_synsignal{i});
    title(['Synthetic ECG + AWGN (SNR: ' num2str(SNR(i)) ' dB)']);
    xlabel('Time (seconds)');
    ylabel('Amplitude (mV)');
    grid on;
end
fft1= fft(white_noise{6});
fft1= fftshift(fft1)
fft2= fft(synECG);
fft2= fftshift(fft2)
N1 = length(white_noise{6});
N2= length(synECG);
frequencies = (-N/2:N/2-1) * fsRealECG/N;
plot(frequencies, abs(fft1)/N1,'r');
hold on;
plot(frequencies, abs(fft2)/N2,'b');
xlabel('Frequency (Hz)');
ylabel('Normalized magnitude');
title ('Frequency Spectrum of AWGN (SNR:5 dB) + Synthetic ECG')
legend ('AWGN (SNR:5 dB)', 'Synthetic ECG')

%Power Line Interference 10 sec
powerFs= 50; 
durationPowerLine = 10; 
peak = max(signaln) - min(signaln);
power_line_amplitude = 0.5 * peak;
interference = power_line_amplitude * sin(2 * pi * powerFs*trealECG);
mean_interference=mean(interference)
 var_interference=var(interference)
 power_interference=var_interference+(mean_interference^2)
t = 0:1/fsRealECG:(durationRealECG-1/fsRealECG);
 Pline_realECG= signaln + interference
 % Pline_synECG= synECG + interference
 figure;
plot(t, Pline_realECG);
title('Powerline Interference + PhysioNet ECG');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
figure;
plot(t, Pline_synECG);
title('Power Line Interference + Synthetic ECG');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
fft1= fft(interference);
fft1= fftshift(fft1)
fft2= fft(synECG);
fft2= fftshift(fft2)
N1 = length(interference);
N2= length(synECG);
frequencies = (-N/2:N/2-1) * fsRealECG/N;
plot(frequencies, abs(fft1)/N1,'r');
hold on;
plot(frequencies, abs(fft2)/N2,'b');
xlabel('Frequency (Hz)');
ylabel('Normalized magnitude');
title ('Frequency Spectrum of Powerline Interference + Synthetic ECG')
legend ('Powerline Interference ', 'Synthetic ECG')
figure;
fs=500
fft_RealECG = fft(Pline_realECG);
fft_RealECG= fftshift(fft_RealECG)
N = length(Pline_realECG);
frequencies = (-N/2:N/2-1) * fsRealECG/N;
plot(frequencies, abs(fft_RealECG)/N,'b');
title('Frequency Spectrum of PhysioNet´s  ECG and Powerline Interference');
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
snr_dB = 10 * log10(power_synECG / power_interference);

%Power Line Interference 10 sec (synECG)
powerFs= 50; 
durationPowerLine = 10;
power_line_amplitude = 1
interference = power_line_amplitude * sin(2 * pi * powerFs*trealECG);
mean_interference=mean(interference)
 var_interference=var(interference)
 power_interference=var_interference+(mean_interference^2)
t = 0:1/fsRealECG:(durationRealECG-1/fsRealECG);
  Pline_synECG= synECG + interference
figure;
plot(t, Pline_synECG);
title('Powerline Interference + Synthetic ECG');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
figure;
fs=500
fft_Pline = fft(interference);
fft_synECG = fft(synECG);
frequencies = (0:length(fft_Pline)-1) * fs/length(fft_Pline);
plot(frequencies, abs(fft_Pline),'r');
hold on;
plot(frequencies, abs(fft_synECG),'b');
title('Frequency Spectrum of Synthetic ECG and Powerline Interference');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Powerline Interference', ' Real ECG');
snr_dB = 10 * log10(power_synECG / power_interference);

%Baseline Wander Noise 10 sec
load('bwm.mat')
x=val
fsBaseline=360
durationBaseLine=10
fsBaseline_new=500
[P,Q] = rat(fsBaseline_new/fsBaseline);
BaseLineNoise = resample(x,P,Q);
BaseLineNoise=(BaseLineNoise-mean(BaseLineNoise))./std(BaseLineNoise);
t = 0:1/fsBaseline_new:(durationBaseLine-1/fsBaseline_new);
figure;
plot(t, BaseLineNoise);
title('Baseline Wander Noise');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
mean_BaseLineNoise=mean(BaseLineNoise)
 var_BaseLineNoise=var(BaseLineNoise)
 power_BaseLineNoise=var_BaseLineNoise+(mean_BaseLineNoise^2)
 BaseLine_RealECG= signaln + BaseLineNoise
 plot(t, BaseLine_RealECG);
title('Baseline Wander Noise + PhysioNet´s ECG');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
fft1= fft(BaseLineNoise);
fft1= fftshift(fft1)
fft2= fft(synECG);
fft2= fftshift(fft2)
N1 = length(BaseLineNoise);
N2= length(signaln);
frequencies = (-N/2:N/2-1) * fsRealECG/N;
plot(frequencies, abs(fft1)/N1,'r');
hold on;
plot(frequencies, abs(fft2)/N2,'b');
xlabel('Frequency (Hz)');
ylabel('Normalized magnitude');
title ('Frequency Spectrum of Baseline Wander Noise + Synthetic ECG')
legend ('Baseline Wander Noise ', 'Synthetic ECG')
 BaseLine_synECG= synECG + BaseLineNoise
 plot(t, BaseLine_synECG);
title('Baseline Wander Noise + Synthetic ECG');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
mean_BaseLineNoise=mean(0.005*BaseLineNoise)
 var_BaseLineNoise=var(0.005*BaseLineNoise)
 power_BaseLineNoise=var_BaseLineNoise+(mean_BaseLineNoise^2)
snr_dB = 10 * log10(power_synECG / power_BaseLineNoise);


 %Baseline Wander Noise 60 sec
load('bwm1min.mat')
x=val
fsBaseline=360
durationBaseLine=60
fsBaseline_new=500
[P,Q] = rat(fsBaseline_new/fsBaseline);
BaseLineNoise = resample(x,P,Q);
t = 0:1/fsBaseline_new:(durationBaseLine-1/fsBaseline_new);
figure;
plot(t, BaseLineNoise);
title('Baseline Wander Noise');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
mean_BaseLineNoise=mean(BaseLineNoise)
 var_BaseLineNoise=var(BaseLineNoise)
 power_BaseLineNoise=var_BaseLineNoise+(mean_BaseLineNoise^2)

t_new = 0:1/fsBaseline_new:(10-1/fsBaseline_new); % Time vector from 0 to 10 seconds
BaseLineNoise_extracted = BaseLineNoise(30*fsBaseline_new+1:40*fsBaseline_new); 
figure;
plot(t_new, BaseLineNoise_extracted);
title('Baseline Wander Noise');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
mean_BaseLineNoise_extracted = mean(0.005*BaseLineNoise_extracted);
var_BaseLineNoise_extracted = var(0.005*BaseLineNoise_extracted);
power_BaseLineNoise_extracted = var_BaseLineNoise_extracted + (mean_BaseLineNoise_extracted^2);
snr_dB = 10 * log10(power_synECG /(power_BaseLineNoise_extracted));
%Muscle Noise 10 sec
load('mam.mat')
x=val
fsMuscle=360
durationMuscle=10
fsMuscle_new=500
[P,Q] = rat(fsMuscle_new/fsMuscle);
muscleNoise = resample(x,P,Q);
muscleNoise=(muscleNoise-mean(muscleNoise))./std(muscleNoise)
t = 0:1/fsMuscle_new:(durationMuscle-1/fsMuscle_new);
figure;
plot(t, muscleNoise);
title('Muscle Noise');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
 mean_muscleNoise=mean(muscleNoise)
 var_muscleNoise=var(muscleNoise)
 power_muscleNoise=var_muscleNoise+(mean_muscleNoise^2)

%Muscle Noise 1 min
load('mam1min.mat')
x=val
fsMuscle=360
durationMuscle=60
fsMuscle_new=500
[P,Q] = rat(fsMuscle_new/fsMuscle);
muscleNoise = resample(x,P,Q);
t = 0:1/fsMuscle_new:(durationMuscle-1/fsMuscle_new);
figure;
plot(t, muscleNoise);
title('Muscle Noise');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;

mean_muscleNoise=mean(muscleNoise)
 var_muscleNoise=var(muscleNoise)
 power_muscleNoise=var_muscleNoise+(mean_muscleNoise^2)

 %Muscle Noise 10 sec EMG Database
load('emg_healthym (5).mat')
x=val
fsMuscleEMG=4000
durationMuscleEMG=10
fsMuscle_new=500
[P,Q] = rat(fsMuscle_new/fsMuscleEMG);
muscleNoiseEMG = resample(x,P,Q);
muscleNoiseEMG=(muscleNoiseEMG-mean(muscleNoiseEMG))./std(muscleNoiseEMG);
t = 0:1/fsMuscle_new:(durationMuscleEMG-1/fsMuscle_new);
figure;
plot(t, muscleNoiseEMG);
title('Muscle Noise');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
 mean_muscleNoise=mean(muscleNoiseEMG)
 var_muscleNoise=var(muscleNoiseEMG)
 power_muscleNoiseEMG=var_muscleNoise+(mean_muscleNoise^2)

%SynECG + muscle noise 
s1=0.5* muscleNoise+ synECG
figure;
plot(t, s1);
title('Synthetic ECG + Muscle Noise');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
snr_dB = 10 * log10(power_synECG / power_muscleNoise);
fft1= fft(muscleNoise);
fft1= fftshift(fft1)
fft2= fft(synECG);
fft2= fftshift(fft2)
N = length(s1);
frequencies = (-N/2:N/2-1) * fsRealECG/N;
plot(frequencies, abs(fft1)/N,'r');
hold on;
plot(frequencies, abs(fft2)/N,'b');
xlabel('Frequency (Hz)');
ylabel('Normalized magnitude');
title ('Frequency Spectrum of Muscle Noise + Synthetic ECG')
legend ('Muscle Noise', 'Synthetic ECG')


%SynECG + muscle noise EMG
s2=muscleNoiseEMG+synECG
figure;
plot(t, s2);
title('Synthetic ECG + Muscle Noise');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
snr_dB = 10 * log10(power_synECG / power_muscleNoiseEMG);

%RealECG + muscle noise EMG
s3=muscleNoiseEMG + realECGFiltered
figure;
plot(trealECG, s3);
title('Real ECG + Muscle Noise');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
snr_dB = 10 * log10(power_realECGFiltered / power_muscleNoiseEMG);

%RealECG + muscle noise
s4=muscleNoise+signaln
figure;
plot(trealECG, s4);
title('Muscle Noise + PhysioNet´s ECG');
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
grid on;
snr_dB = 10 * log10(power_realECGFiltered / power_muscleNoise);
fft1= fft(muscleNoise);
fft1= fftshift(fft1)
fft2= fft(signaln);
fft2= fftshift(fft2)
N = length(s4);
frequencies = (-N/2:N/2-1) * fsRealECG/N;
plot(frequencies, abs(fft1)/N,'r');
hold on;
plot(frequencies, abs(fft2)/N,'b');
xlabel('Frequency (Hz)');
ylabel('Normalized magnitude');
title ('Frequency Spectrum of Muscle Noise + PhysioNet´s ECG')
legend ('Muscle Noise', 'PhysioNet´s ECG')
