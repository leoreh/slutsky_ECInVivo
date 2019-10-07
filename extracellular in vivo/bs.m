ch = 16;
fs = 1250;
type = 'butter';
order = 3; 
winLength = 11;                 % [samples]    

% filter    
passband = [1 100];
sigfilt = filterLFP(double(lfp.data(:, ch)), 'type', type,...
    'passband', passband, 'order', order, 'graphics', false);

%%% Time-frequency analysis
idx = 30*60 * fs : 35*60 * fs;
% stft
figure; spectrogram(sigfilt(idx), 1000, [], [0 : 1 : 30], fs, 'yaxis')
% cwt
figure; cwt(sigfilt(idx), 'morse', fs, 'FrequencyLimits', [1 500])

%%% single IIS analysis
% single IIS
dur = fs * 0.8; % burst duration in samples
idx = 2582350 : 2582650;        % iis
idx = 323000 : 323600;          % burst
iis = lfp.data(idx, ch);
figure;
plot(idx / fs, iis)
axis tight

% filter    
passband = [30 500];
x = filterLFP(double(iis), 'type', type,...
    'passband', passband, 'order', order, 'graphics', false);

% fft
figure
l = length(x);
y = fft(x);
P2 = abs(y/l);
P1 = P2(1:l/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(l/2))/l;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


%%% detect BS
passband = [2 100];
x = filterLFP(double(lfp.data(:, 16)), 'order', 3, 'type', 'butter',...
    'passband', passband, 'graphics', false);

% x = double(lfp.data(:, 16));
t = [1 : length(x)] / fs;
winLength = 30;

% rectify
y = x .^ 2;
% moving average
y = movmean(y, winLength);
% standerdize
y = (y - mean(y(:))) / std(y(:));

% y = fastrms(x);
% y = (y - mean(y(:))) / std(y(:));

figure;
plot(t, x)
hold on
yyaxis right
p = plot(t, y);
p.Color(4) = 0.3;
% axis tight
xlim([1800 2200])

idx = 25*60 * fs : 50*60 * fs;
figure; spectrogram(x(idx), 10000, [], [0 : 10 : 200], fs, 'yaxis')
figure; cwt(x(idx), 'morse', fs, 'FrequencyLimits', [1 500])

bs = findBS('lfp', lfp, 'emgThr', 0, 'basepath', basepath, 'ch', 16);
