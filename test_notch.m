%% Read and display the signal of original signal
% Read the audio file (audio represents the signal plot)
[audio, fs] = audioread("audio09.wav");

%Play the audio file
%sound(audio, fs);
%pause(length(audio)/fs + 2);
%% Improved IIR Filter Design and Application
% Sampling frequency
fs = 11025; 

% Frequency to be notched out (240 Hz)
f_notch = 240; 
% Normalized notch frequency (240 Hz)
f_notch_norm = f_notch / (fs/2); 

% Design a notch filter using the IIR notch filter design function
wo = f_notch_norm; % Normalized frequency
bw = wo / 35; % Bandwidth, adjust as necessary
[b_notch, a_notch] = iirnotch(wo, bw);

% Apply the notch filter to the audio signal
audio_notched = filter(b_notch, a_notch, audio);

% Design a more aggressive low-pass IIR filter
Fc = 2000; % Lower cutoff frequency for more aggressive filtering
Fc_norm = Fc / (fs/2); % Normalized cutoff frequency
[b_lowpass, a_lowpass] = butter(6, Fc_norm, 'low'); % 6th order Butterworth filter

% Apply the low-pass filter to the notch-filtered audio signal
audio_iir = filter(b_lowpass, a_lowpass, audio_notched);

% Plot the frequency response of the notch filter
figure;
freqz(b_notch, a_notch, 1024, fs);
title('Frequency Response of Notch Filter');

% Plot the frequency response of the low-pass filter
figure;
freqz(b_lowpass, a_lowpass, 1024, fs);
title('Frequency Response of Low-Pass Filter');

% Plot the filtered audio signal (Notch + Low-Pass)
figure;
t = (1:length(audio_iir)) / fs;
plot(t, audio_iir);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Waveform of Notch and Low-Pass Filtered audio09.wav');
grid on;

% Compute the FFT of the filtered audio signal
Y_filtered = fft(audio_iir);

% Normalize the FFT output
Y_filtered = Y_filtered / length(audio_iir);

% Get the magnitude of the FFT
Y_filtered_magnitude = abs(Y_filtered);

% Use fftshift to center the FFT spectrum
Y_filtered_magnitude_shifted = fftshift(Y_filtered_magnitude);

% Create a frequency vector
f = linspace(-fs/2, fs/2, length(Y_filtered_magnitude));

% Plot the frequency spectrum of the filtered audio signal
figure;
plot(f, Y_filtered_magnitude_shifted);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Notch and Low-Pass Filtered audio09.wav (Shifted)');
grid on;

% Play the IIR filtered audio signal
%sound(audio_iir, fs);

%% Additional Evaluation Metrics

% Compute SNR of the original signal (SNR = signal-noise ratio)
snr_original = snr(audio);

% Compute SNR of the IIR filtered signal
snr_iir = snr(audio_iir);

% Calculate Mean Squared Error (MSE) for IIR filtered signal
mse_iir = sum((audio - audio_iir).^2) / length(audio);

% Display SNR and MSE values
fprintf('SNR of original signal: %.2f dB\n', snr_original);
fprintf('SNR of IIR filtered signal: %.2f dB\n', snr_iir);
fprintf('MSE of IIR filtered signal: %.6f\n', mse_iir);

% Use spectrogram to differentiate the signals plots
% Original audio spectrogram
figure;
spectrogram(audio, 256, [], [], fs, 'yaxis');
title('Spectrogram of Original Audio');

% IIR filtered audio spectrogram
figure;
spectrogram(audio_iir, 256, [], [], fs, 'yaxis');
title('Spectrogram of Improved IIR Filtered Audio');

audiowrite('audio09_IIR_notch.wav', audio_iir, fs);