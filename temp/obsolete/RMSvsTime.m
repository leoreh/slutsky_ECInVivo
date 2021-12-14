% plots RMS, 50Hz and EMG as a function of time


% load data
basepath = 'E:\Data\Data Tanks\LH\DataTanks\LH4\LH4_13sep19';
[~, filename, ~] = fileparts(basepath);
filename = [filename '.dat'];
fs = 24414.125;
nCh = 16;
ch = 1 : 16;
dur = fs * 600;
dur = Inf;

dat = bz_LoadBinary(filename, 'frequency', fs, 'nChannels', nCh, 'channels', ch, 'duration', dur);


% calc RMS
win = 0.01;  % s
idx = 1 : fs * 1;
mRMS = dsp.MovingRMS(round(fs * win));
y = mRMS(single(dat(idx, :)));


figure
t = (1 : length(idx)) / fs;
plot(t, dat(idx, i))
hold on
plot(t, y(idx, i), 'Color', 'k', 'LineWidth', 3)

% Time-frequency analysis
y = stft(dat(idx));

spectrogram(single(dat(idx, 1)), fs / 10, [], 1 : 100, fs, 'yaxis')
