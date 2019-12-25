%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch = 5;
basepath = 'E:\Data\Dat\lh43';
[~, filename] = fileparts(basepath);
cd(basepath)
force = false;

if force
    load([filename '.lfp.mat'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inspect lfp 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = lfp.fs;
idx = 1 : 60 * fs;
x = lfp.data(idx, ch);
figure
plot(lfp.timestamps(idx), x)
xlim([idx(1), idx(1) + 10])

% for sleep-state analysis it is probably preferable to calculate the power
% spectrum (dB) rather than PSD (dB / Hz) because it is independent of
% other frequencies. however, when tested no differences were found.

win = 0.5 * fs;
freq = [0 : 100];

figure
[s, f, t, p] = spectrogram(x, win, round(win / 4), freq, fs, 'yaxis', 'power', 'onesided');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bimodal separation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
winmax = 0.8 * fs;      % window for moving max        
binw = 0.15 * fs;       % bin width for std calculation
dist = 0.2 * fs;        % distance for removing adjacent peaks
thr = 1;                % threshold of detection

% std
s = movstd(double(x), binw);


figure,plot(s)

% threshold for bs separation is determined according to the bimodel
% distribution of std values
thr = sepBimodel('x', s, 'lognorm', true, 'graphics', true);


histogram(s)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getBS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bs = getBS('sig', x, 'fs', lfp.fs, 'basepath', basepath,...
    'graphics', false, 'saveVar', false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(lfp.timestamps(idx), bs.binary)
hold on
yyaxis right
plot(lfp.timestamps(idx), lfp.data(idx, ch))
xlim([idx(1), idx(1) + 300])

dur = diff(bs.stamps, 1, 2);
ibi = bs.stamps(2 : end, 1) - bs.stamps(1 : end - 1, 2);

histogram(log10(dur))

figure, plot(log10(ibi), '*')
