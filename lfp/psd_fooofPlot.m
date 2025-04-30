function fh = psd_fooofPlot(psd_1of, varargin)

% plots the power spectral density and FOOOF model fits for specified bouts
% or the mean across all bouts.
%
% 02 feb 25 LH 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'bout_idx', [], @(x) isnumeric(x) || isempty(x));
addOptional(p, 'istate', 1, @(x) isnumeric(x) || isempty(x));

parse(p, varargin{:});
bout_idx = p.Results.bout_idx;
istate = p.Results.istate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get number of bouts if bout_idx is empty
nbouts = size(psd_1of.psd_orig, 2);
if isempty(bout_idx)
    bout_idx = 1 : nbouts;
end

% extract relevant data
freqs = psd_1of.freqs;
psd_orig = squeeze(psd_1of.psd_orig(istate, bout_idx, :));
psd_fit = squeeze(psd_1of.psd_fit(istate, bout_idx, :));
psd_ap = squeeze(psd_1of.psd_ap(istate, bout_idx, :));
residuals = squeeze(psd_1of.residuals(istate, bout_idx, :));
cf = squeeze(psd_1of.cf(istate, bout_idx, :));
amp = squeeze(psd_1of.amp(istate, bout_idx, :));
sd = squeeze(psd_1of.sd(istate, bout_idx, :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open figure
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [1, 2];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'FOOOF', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% plot spectra
axh = nexttile(th, 1, [1, 1]); cla; hold on
% plot mean across bouts using stdshade
ph(1) = plot_stdShade('dataMat', psd_orig', 'xVal', freqs, 'axh', axh, ...
    'clr', [0 0 0], 'alpha', 0.3);
ph(2) = plot_stdShade('dataMat', psd_fit', 'xVal', freqs, 'axh', axh, ...
    'clr', [1 0 0], 'alpha', 0.3);
ph(3) = plot_stdShade('dataMat', psd_ap', 'xVal', freqs, 'axh', axh, ...
    'clr', [0 0 1], 'alpha', 0.3);
legNames = {'Original', 'FOOOF fit', 'Aperiodic fit'};

% format
ylim([1 5])
set(axh, 'xscale', 'log')
xlabel('Frequency (Hz)')
ylabel('Power (log10)')
title(axh, 'FOOOF Model Fit')
legend(axh, ph, legNames)
grid on

% plot peaks
axh = nexttile(th, 2, [1, 1]); cla; hold on
% plot mean residuals using stdshade
ph(1) = plot_stdShade('dataMat', residuals', 'xVal', freqs, 'axh', axh, ...
    'clr', [0 0 0], 'alpha', 0.3);
legNames = cell(1, 1);
legNames{1} = 'Residuals';

% plot mean peaks
for ipeak = 1 : size(cf, 2)
    cf_mean = mean(cf(:, ipeak), 'omitnan');
    amp_mean = mean(amp(:, ipeak), 'omitnan');
    sd_mean = mean(sd(:, ipeak), 'omitnan');

    gaussian = amp_mean * exp(-((freqs - cf_mean) .^ 2 / (2 * sd_mean ^ 2)));
    ph(ipeak + 1) = plot(freqs, gaussian, '--', 'LineWidth', 1.5);
    legNames{ipeak + 1} = sprintf('Peak %.1f Hz', cf_mean);
end

set(axh, 'xscale', 'log')
xlabel('Frequency (Hz)')
ylabel('Power (log10)')
title(axh, 'Isolated Peaks')
grid on
legend(axh, ph, legNames)

end
% EOF