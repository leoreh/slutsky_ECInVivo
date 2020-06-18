%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter LFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

order = 4;
passband = 'gamma';
switch passband
    case 'delta'
        passband = [0 6];
        order = 8;
    case 'theta'
        passband = [4 10];
    case 'spindles'
        passband = [9 17];
    case 'gamma'
        passband = [30 80];
    case 'ripples'
        passband = [100 250];
    otherwise
        error('Unknown frequency band')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find ripples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter
passband = [130 200];          
order = 3; 
type = 'butter';
filtered = filterLFP(double(lfp.data), 'type', type,...
    'passband', passband, 'order', order, 'graphics', false);

ripples = findRipples(lfp);
[maps,data,stats] = bz_RippleStatsLH(sigfilt(:, 5), lfp.timestamps, ripples);
bz_PlotRippleStats(maps, data, stats)

[wt,f] = cwt(double(lfp.data([1 : 100000, 5])), lfp.fs);
idx(1) = find(lfp.timestamps == ripples.timestamps(1, 1));
idx(2) = find(lfp.timestamps == ripples.timestamps(1, 2));
figure, plot(lfp.data(idx(1) : idx(2), 5))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate and plot spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frange = [1 128];
nf = 100;       
f = logspace(log10(frange(1)), log10(frange(2)), nf);
win = 1 * fs;
noverlap = 0.8 * fs;

% calculate
[spec, ~, t_FFT] = spectrogram(x, win, noverlap, f, fs);
spec = abs(spec)';

% plot
spectrogram(x, 'yaxis', win, noverlap, f, fs)
ax = gca;
ax.YScale = 'log';







