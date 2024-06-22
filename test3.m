%% Read and display the signal of original signal
% Read the audio file (audio represents the signal plot)
[audio, fs] = audioread("audio09.wav");

%Play the audio file
sound(audio, fs);
pause(length(audio)/fs + 2);

% Display the audio information
disp('Sample Rate:');
disp(fs);

%% Plot the original signal graph
% Get the duration of the audio signal in seconds
duration = length(audio) / fs; 

% Create a time vector
t = linspace(0, duration, length(audio));

% Plot the original audio signal
figure;
plot(t, audio);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Waveform of audio09.wav');
grid on;

%% Plot frequency spectrum of the graph
% Compute the FFT of the audio signal
Y = fft(audio);

% Normalize the FFT output
Y = Y / length(audio);

% Get the magnitude of the FFT
Y_magnitude = abs(Y);

% Only keep the first half of the FFT output (positive frequencies)
Y_magnitude = Y_magnitude(1:floor(length(Y_magnitude)/2));

% Create a frequency vector
f = linspace(0, fs/2, length(Y_magnitude));

% Plot the frequency spectrum
figure;
plot(f, Y_magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of audio09.wav');
grid on;

%% Improved FIR Filter Design and Application
% Filter parameters
N = 180; % Increased Filter order
Fc = 0.3; % Cutoff frequency (normalized, 0 to 1, where 1 is the Nyquist frequency)

% Try different window functions
%b_fir = fir1(N, Fc, hamming(N+1)); % FIR filter coefficients using Hamming window
%b_fir = fir1(N, Fc, hann(N+1)); % Using Hanning window
%b_fir = fir1(N, Fc, blackman(N+1)); % Using Blackman window
b_fir = fir1(N, Fc, kaiser(N+1, 0.5)); % Using Kaiser window with beta parameter

% Plot the frequency response of the FIR filter using freqz
figure;
freqz(b_fir, 1, 1024, fs);
title('Frequency Response of Improved FIR Filter');

% Apply FIR filter to the audio signal
audio_fir = filter(b_fir, 1, audio);

% Plot the filtered audio signal (FIR)
figure;
plot(t, audio_fir);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Waveform of Improved FIR Filtered audio09.wav');
grid on;

% Compute the FFT of the FIR filtered audio signal
Y_fir = fft(audio_fir);

% Normalize the FFT output
Y_fir = Y_fir / length(audio_fir);

% Get the magnitude of the FFT
Y_fir_magnitude = abs(Y_fir);

% Use fftshift to center the FFT spectrum
Y_fir_magnitude_shifted = fftshift(Y_fir_magnitude);

% Create a frequency vector
f = linspace(-fs/2, fs/2, length(Y_fir_magnitude));

% Plot the frequency spectrum of the FIR filtered audio signal
figure;
plot(f, Y_fir_magnitude_shifted);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Improved FIR Filtered audio09.wav (Shifted)');
grid on;


%% Improved IIR Filter Design and Application
% Parameters for IIR filter
N = 2; % Filter order (adjusted to avoid excessive ringing)
Fc = 0.3; % Cutoff frequency (normalized, 0 to 1, where 1 is the Nyquist frequency)

% Design IIR Butterworth filter
[b_iir, a_iir] = butter(N, Fc);

% Plot the frequency response of the IIR filter using freqz
figure;
freqz(b_iir, a_iir, 1024, fs);
title('Frequency Response of Improved IIR Filter');

% Apply IIR filter to the audio signal
audio_iir = filter(b_iir, a_iir, audio);

% Plot the filtered audio signal (IIR)
figure;
plot(t, audio_iir);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Waveform of Improved IIR Filtered audio09.wav');
grid on;

% Compute the FFT of the IIR filtered audio signal
Y_iir = fft(audio_iir);

% Normalize the FFT output
Y_iir = Y_iir / length(audio_iir);

% Get the magnitude of the FFT
Y_iir_magnitude = abs(Y_iir);

% Use fftshift to center the FFT spectrum
Y_iir_magnitude_shifted = fftshift(Y_iir_magnitude);

% Create a frequency vector
f = linspace(-fs/2, fs/2, length(Y_iir_magnitude));

% Plot the frequency spectrum of the IIR filtered audio signal
figure;
plot(f, Y_iir_magnitude_shifted);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Improved IIR Filtered audio09.wav (Shifted)');
xlim([0 500]); % Zoom in on the frequency range from 0Hz to 500Hz
grid on;


%% Play the edited audio after FIR and IIR are applied
% Play the IIR filtered audio signal
sound(audio_iir, fs);

pause(length(audio_iir)/fs + 2);

% Play the FIR filtered audio signal
sound(audio_fir, fs);

%% Additional Evaluation Metrics

% Compute SNR of the original signal (SNR = signal-noise ratio)
snr_original = snr(audio);

% Compute SNR of the FIR filtered signal
snr_fir = snr(audio_fir);

% Compute SNR of the IIR filtered signal
snr_iir = snr(audio_iir);

% Calculate Mean Squared Error (MSE) for FIR filtered signal
mse_fir = sum((audio - audio_fir).^2) / length(audio);

% Calculate Mean Squared Error (MSE) for IIR filtered signal
mse_iir = sum((audio - audio_iir).^2) / length(audio);

% Display SNR and MSE values
fprintf('SNR of original signal: %.2f dB\n', snr_original);
fprintf('SNR of FIR filtered signal: %.2f dB\n', snr_fir);
fprintf('SNR of IIR filtered signal: %.2f dB\n', snr_iir);
fprintf('MSE of FIR filtered signal: %.6f\n', mse_fir);
fprintf('MSE of IIR filtered signal: %.6f\n', mse_iir);

% Use spectrogram to differentiate the signals plots
% Original audio spectrogram
figure;
spectrogram(audio, 256, [], [], fs, 'yaxis');
title('Spectrogram of Original Audio');

% FIR filtered audio spectrogram
figure;
spectrogram(audio_fir, 256, [], [], fs, 'yaxis');
title('Spectrogram of Improved FIR Filtered Audio');

% IIR filtered audio spectrogram
figure;
spectrogram(audio_iir, 256, [], [], fs, 'yaxis');
title('Spectrogram of Improved IIR Filtered Audio');


%% Save the edited audio
% From FIR design
audiowrite('audio09_FIR.wav', audio_fir, fs);

% From IIR design
audiowrite('audio09_IIR.wav', audio_iir, fs);
