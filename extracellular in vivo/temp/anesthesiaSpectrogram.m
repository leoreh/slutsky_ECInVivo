% anesthesiaSpectrogram

dband = specBand('sig', sig, 'graphics', true, 'band', [1 4]);

% broad-band spectrogram
% instead of 10 s window w/ 1 s overlap, 0 overlap and smooth with 10 points
% this is so bins of bs and p are congruent
freq = linspace(1, 4, 100);
freq = logspace(0, 2, 100);
winsize = 10;       % win length [s]
win = hann(2 ^ nextpow2(winsize * fs));
[s, f, t, p] = spectrogram(sig, win, 0,...
    freq, fs, 'yaxis', 'psd');

% psd
p = 10 * log10(abs(p));

% smooth and z-score
smf = winsize;
z = zscore(movmean(p, smf, 2));
bs.bsr = movmean(bs.bsr, smf);

% integrate power over delta and sigma band
deltaf = [1 4];
[~, deltaidx] = min(abs(f - deltaf));
zdelta = sum(z(deltaidx(1) : deltaidx(2), :), 1);
sigmaf = [9 25];
[~, sigmaidx] = min(abs(f - sigmaf));
zsigma = sum(z(sigmaidx(1) : sigmaidx(2), :), 1);

zdelta = bz_NormToRange(zdelta, [0 1]);
zsigma = bz_NormToRange(zsigma, [0 1]);

% zdelta = movmean(zdelta, 15);

ff = figure;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
suptitle(basename)

% raw
s1 = subplot(3, 1, 1);
plot((1 : length(sig)) / fs / 60, sig)
hold on
axis tight
Y = ylim;
fill([bs.stamps fliplr(bs.stamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
    'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
ylabel('Voltage [mV]')
set(gca, 'TickLength', [0 0])
box off
title('Raw and BS')
if isfield(info, 'transitions')
    addLns('lns', info.transitions / fs / 60,...
        'lbs', {'1%', '1.5%', '2%'})
end

% spectrogram
s2 = subplot(3, 1, 2);
surf(t / 60, f, p, 'EdgeColor', 'none');
axis xy;
axis tight;
view(0,90);
origSize = get(gca, 'Position');
colormap(jet);
colorbar;
ylabel('Frequency [Hz]');
set(gca, 'YScale', 'log')
set(s2, 'Position', origSize);
set(gca, 'TickLength', [0 0])
box off
title('Wideband spectrogram')

% delta power
s3 = subplot(3, 1, 3);
hold on
plot(t / 60, zdelta, 'r')
plot(t / 60, zsigma, 'b')
plot(bs.cents / fs / 60, bs.bsr, 'k', 'LineWidth', 3)
legend({'[1-4 Hz]', '[9-25 Hz]', 'BSR'})
axis tight
ylim([0 1])
if isfield(info, 'transitions')
    addLns('lns', info.transitions / fs / 60,...
        'lbs', {'1%', '1.5%', '2%'})
end
xlabel('Time [min]');
ylabel('[a.u.]')
set(gca, 'TickLength', [0 0])
box off
title('Delta and Sigma power')

linkaxes([s1, s2, s3], 'x');

if saveFig
    figname = [basename '_anesthesia'];
    export_fig(figname, '-tif', '-transparent')
    % savePdf(figname, basepath, ff)
end
