function fh = fooof_plotMdl(f1f, varargin)

% plots the power spectral density and FOOOF model fits for specified bouts
% or the mean across all bouts. note 
%
% 02 feb 25 LH 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'bout_idx', [], @(x) isnumeric(x) || isempty(x));

parse(p, varargin{:});
bout_idx = p.Results.bout_idx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get number of bouts if bout_idx is empty
nbouts = size(f1f.psd_orig, 1);
if isempty(bout_idx)
    bout_idx = 1 : nbouts;
end

% extract relevant data
freqs = f1f.freqs;
psd_orig = squeeze(f1f.psd_orig(bout_idx, :));
psd_fit = squeeze(f1f.psd_fit(bout_idx, :));
psd_ap = squeeze(f1f.psd_ap(bout_idx, :));
psd_res = log10(psd_orig) - log10(psd_ap);
squeeze(f1f.psd_res(bout_idx, :));
cf = squeeze(f1f.peaks.cf(bout_idx, :));
amp = squeeze(f1f.peaks.amp(bout_idx, :));
sd = squeeze(f1f.peaks.sd(bout_idx, :));

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
ylim(10.^[1 5])
set(axh, 'xscale', 'log', 'yscale', 'log')
xlabel('Frequency (Hz)')
ylabel('Power')
title(axh, 'FOOOF Model Fit')
legend(axh, ph, legNames)
grid on

% plot peaks
axh = nexttile(th, 2, [1, 1]); cla; hold on
set(axh, 'xscale', 'log', 'yscale', 'linear')
% plot mean residuals using stdshade
ph(1) = plot_stdShade('dataMat', psd_res', 'xVal', freqs, 'axh', axh, ...
    'clr', [0 0 0], 'alpha', 0.3);
legNames = cell(1, 1);
legNames{1} = 'Residuals';
set(axh, 'xscale', 'log', 'yscale', 'linear')

% plot mean peaks
for ipeak = 1 : size(cf, 2)
    cf_mean = mean(cf(:, ipeak), 'omitnan');
    amp_mean = mean(amp(:, ipeak), 'omitnan');
    sd_mean = mean(sd(:, ipeak), 'omitnan');

    gaussian = (amp_mean * exp(-((freqs - cf_mean) .^ 2 / (2 * sd_mean ^ 2))));
    ph(ipeak + 1) = plot(freqs, gaussian, '--', 'LineWidth', 1.5);
    legNames{ipeak + 1} = sprintf('Peak %.1f Hz', cf_mean);
end

xlabel('Frequency (Hz)')
ylabel('Power (log10)')
title(axh, 'Isolated Peaks')
grid on
legend(axh, ph, legNames)

end
% EOF