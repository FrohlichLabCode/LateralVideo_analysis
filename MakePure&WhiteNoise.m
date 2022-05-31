
% for log spaceing: logspace(log10(startFreq),log10(endFreq),numFreq)

%% make pure tone
cf = 523.3;                 % carrier frequency (Hz)
fs = 96000;                 % sample frequency (Hz)
dur = 5;                    % duration (s)
n = fs * dur;                 % number of samples
s = (1:n) / fs;             % sound data preparation
s = sin(2 * pi * cf * s);   % sinusoidal modulation
sound(s, fs);               % sound presentation

audiowrite(['pure_' num2str(cf) 'Hz_' num2str(dur) 'sec.wav'],s,fs)

%% make white noise

fs = 96000;
dur = 5; % seconds
power = -13; % -13 is 80dB

wn = wgn(1,fs*dur,power);
sound(wn', fs);

audiowrite(['whitenoise_' num2str(dur) 'sec.wav'],wn,fs)