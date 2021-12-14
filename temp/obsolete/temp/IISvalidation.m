extension = 'abf';
loadname = [basename '.' extension];

[raw, info] = abf2load(loadname);
fs_orig = 1 / (info.fADCSequenceInterval / 1000000);

minfactor = 29.994;
binsize = minfactor * fs_orig;

fs = 1250;
ds = resample(raw, fs, round(fs_orig));
fprintf('\n resampling from %.1f to %.1f\n\n', fs_orig, fs)

iis = getIIS('sig', -raw, 'fs', fs_orig, 'basepath', basepath,...
    'graphics', false, 'saveVar', false, 'binsize', binsize,...
    'marg', marg, 'basename', basename, 'thr', thr, 'smf', 1,...
    'saveFig', false, 'forceA', true, 'spkw', false, 'vis', true);

% x = data from clampfit
[c, ia, ic] = unique(x);
ii = 1 : max(c);
ii(c) = [];
counts = [c, accumarray(ic, 1)];
counts = [counts; [ii', zeros(length(ii), 1)]];
[~, ii] = sort(counts);
counts = counts(ii(:, 1), :);
histogram(round(iis.rate) - counts(:, 2), 10);

% compare with iis.rate from fs = 1250
figure
plot(iis.rate, 'b')
xlabel('Matlab')
axis tight
ax1 = gca;
ax1.XColor = 'b';
ax1.YColor = 'b';
ax1_pos = ax1.Position;
ax2 = axes('Position', ax1_pos,...
    'XAxisLocation', 'top',...
    'YAxisLocation', 'right',...
    'Color', 'none');
line(counts(:, 1), counts(:, 2), 'Parent', ax2, 'Color', 'k')
xlim([1 337])
xlabel('Clampfit')

idx = round(329 * fs_orig * minfactor) : round(330 * fs_orig * minfactor);
%  idx = length(lfpraw.data) - round(1 * fs_orig * minfactor) : length(lfpraw.data);
t = linspace(0, minfactor, length(idx));

figure
plot(t, raw(idx));
axis tight
hold on
plot(xlim, [thr(2) thr(2)], '--')


idx2 = round(329 * fs_orig * minfactor) : round(330 * fs_orig * minfactor);
t = linspace(idx2(1) / fs_orig / minfactor, idx2(end) / fs_orig / minfactor, length(idx2));

figure
plot(t, raw(idx2));
hold on
axis tight
plot(xlim, [-thr(2) -thr(2)], '--')
scatter(iis.peakPos(iis.peakPos > idx2(1) & iis.peakPos < idx2(end)) / fs_orig / minfactor,...
    -iis.peakPower(iis.peakPos > idx2(1) & iis.peakPos < idx2(end)), '*');
