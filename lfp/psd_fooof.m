function psd_1of = psd_fooof(psd_data, freqs, varargin)

% Parameterizes neural power spectra using FOOOF (Fitting Oscillations &
% One-Over F), separating periodic (peaks) and aperiodic (1/f-like)
% components. Uses the MATLAB wrapper for Python FOOOF implementation.
%
% in windows cmd, make sure python includes relevant packages
% py -3.12 -m pip install numpy
% py -3.12 -m pip install fooof
% py -3.12 -m pip install matplotlib
%
% INPUT
%   psd_data    power spectra to analyze [bout x freq]. will analyze each
%               bout seperately. can also be a cell (e.g., per state) at
%               which case the output will be organized accordingly
%   freqs       vector of frequencies [1 x nfreq]
%   basepath    char. recording session path {pwd}
%   f_range     frequency range for fitting [1 100]
%   fooof_cfg   struct with optional fields for configuration. see below.
%   flg_pli     logical. interpret freqs of power line interference (true)
%   flg_plot    logical. plot fooof results per bout (true)
%   saveVar     logical / char. save variable {true}. if char then variable
%               will be named [saveVar].mat
%
% OUTPUT
%   psd_1of             struct array [nstates/nbouts x 1]. see end of
%                       function
%
% CALLS
%   fooof
%
% DEPENDENCIES
%   Python with packages installed:
%   py -3.12 -m pip install numpy fooof matplotlib
%
% SEE ALSO
%   https://fooof-tools.github.io/fooof/
%
% 09 Jan 24 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'fooof_cfg', []);
addOptional(p, 'basepath', pwd);
addOptional(p, 'f_range', [1 100]);
addOptional(p, 'flg_pli', true, @islogical);
addOptional(p, 'flg_plot', true, @islogical);
addOptional(p, 'saveVar', true);

parse(p, varargin{:})
fooof_cfg           = p.Results.fooof_cfg;
f_range             = p.Results.f_range;
basepath            = p.Results.basepath;
flg_pli             = p.Results.flg_pli;
flg_plot            = p.Results.flg_plot;
saveVar             = p.Results.saveVar;

% append info to output
psd_1of.info.fooof_cfg = fooof_cfg;
psd_1of.info.flg_pli = flg_pli;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename
[~, basename] = fileparts(basepath);
if ischar(saveVar)
    fname = [basepath, filesep, basename, '.' saveVar '.mat'];
else
    fname = fullfile(basepath, [basename, '.psd_1of.mat']);
end

% ensure row vectors
freqs = freqs(:)';

if iscell(psd_data)
    nstates = length(psd_data);
else
    nstates = 1;
    psd_data = {psd_data};
end

% initialize python and packages if needed
try         % check if modules already imported
    py.importlib.import_module('numpy');
    py.importlib.import_module('fooof');
    py.importlib.import_module('matplotlib');
catch ME
    clear classes
    clear py.*
    pyenv;   % initialize python environment
    pyversion
    % Verify FOOOF is installed in that Python environment
    py.help('fooof')
    py.importlib.import_module('numpy');
    py.importlib.import_module('fooof');
    py.importlib.import_module('matplotlib');
end

% model parameters
if isempty(fooof_cfg)
    % defaults
    % 'peak_width_limits', [0.5 12],...    % [min max] allowed peak width
    % 'max_n_peaks', inf,...               % maximum number of peaks to fit
    % 'min_peak_height', 0.0,...           % minimum peak height threshold
    % 'peak_threshold', 2,...              % relative height above background
    % 'aperiodic_mode', 'fixed',...        % 'fixed' or 'knee'
    % 'verbose', true,...                  % print output
    % 'return_model', true);               % return model results

    fooof_cfg = struct(...
        'peak_width_limits', [2 50],...
        'max_n_peaks', 6,...
        'min_peak_height', 0.2,...
        'peak_threshold', 1.5,...
        'aperiodic_mode', 'knee',...
        'verbose', true,...
        'return_model', true);
end

% test models with different peak numbers. adapted from
% https://www.biorxiv.org/content/10.1101/2024.08.01.606216v1.full.pdf
% (brainstorm). In sum, the idea is to create models with 0 to
% max_peaks and quantify the Bayesian Information Criterion for each
% model. The model with the lowest BIC is then subjected to Bayes factor
% inference against the aperiodic spectral model to adjudicate whether
% spectral peaks are likely. A Bayes factor < 1 is evidence in
% favour of periodic brain activity over the null hypothesis of no periodic
% activity. Note in their paper they implement this more efficiently, by
% adding peaks to the aperiodic component rather then remodeling the data.
% an example call to their function can be found at the end of this script.

% here, a model will be created for each combination of parameters below.
npeaks = [0 : fooof_cfg.max_n_peaks];
nmdls = length(npeaks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run fooof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
fooof_state = cell(1, nstates);

for istate = 1 : nstates

    % get number of bouts
    [nbouts, nfreq] = size(psd_data{istate});

    % validate psd input
    if nfreq ~= length(freqs)
        error('chech psd_data input')
    end

    % interpolate power line interference
    if flg_pli
        psd_mat = interp_pli(psd_data{istate}, freqs);
    end

    % initialize output
    fooof_data = cell(nbouts, 1);

    % loop through bouts
    for ibout = 1 : nbouts
        for imdl = 1 : nmdls

            % adjust model parameters
            fooof_cfg.max_n_peaks = npeaks(imdl);

            % run fooof
            fmdl(imdl) = fooof(freqs, psd_mat(ibout, :),...
                f_range, fooof_cfg, fooof_cfg.return_model);

            % calculate model fit
            ffit(imdl) = calc_fit(fmdl(imdl));

            % force distance between peaks
            thr_dist = 1;
            fmdl(imdl) = fooof_rmPeaks(fmdl(imdl), thr_dist);

            % calculate model fit
            ffit(imdl) = calc_fit(fmdl(imdl));
        end

        % select best model
        bic = vertcat(ffit.bic);
        [~, idx_best] = min(bic);

        % Bayes Factor between best model and apreiodic component
        bf = exp((bic(idx_best) - bic(1)) ./ 2);
        if bf > 1
            idx_best = 1;
        end

        % append
        fooof_data{ibout} = fmdl(idx_best);
        fooof_data{ibout}.fit = ffit(idx_best);

        % plot
        if flg_plot
            plot_fooof(fooof_data{ibout});
        end
    end

    fooof_state{istate} = cat(1, fooof_data{:});
end

% organize output [state x bout x frequency] or [state x bout x band]
psd_1of = fooof2params(fooof_state);

% save
if saveVar
    save(fname, 'psd_1of')
end

end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove overlapping peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fmdl, flg_rm] = fooof_rmPeaks(fmdl, thr_prox)
% Removes overlapping peaks from FOOOF model, keeping the higher amplitude peak.
% Overlap is determined by proximity threshold in standard deviations.

% Initialize output flag
flg_rm = false;

% Get peak parameters
peaks = fmdl.gaussian_params;
if ~any(peaks)
    return
end

% Sort peaks by frequency
[~, idx] = sort(peaks(:, 1));
peaks = peaks(idx, :);

% Calculate frequency-dependent threshold
freq_factor = log(peaks(:, 1) / min(peaks(:, 1)) + 1);  % increases with frequency
scaled_thresh = thr_prox * freq_factor;

% Calculate bounds for checking overlap with frequency-scaled threshold
bounds = [peaks(:, 1) - peaks(:, 3) .* scaled_thresh,...
    peaks(:, 1),...
    peaks(:, 1) + peaks(:, 3) .* scaled_thresh];

% Find peaks to remove
idx_rm = [];
for ipeak = 1 : size(bounds, 1) - 1

    % Get current and subsequent peak bounds
    b0 = bounds(ipeak, :);
    b1 = bounds(ipeak + 1:end, :);

    % Check if current peak extends into any subsequent peaks
    overlaps = b0(2) > b1(:, 1);

    if any(overlaps)
        % Get indices of overlapping peaks
        idx_overlaps = find(overlaps) + ipeak;

        % Compare amplitudes with each overlapping peak
        for idx_compare = idx_overlaps'
            if peaks(ipeak, 2) < peaks(idx_compare, 2)
                idx_rm = [idx_rm, ipeak];
            else
                idx_rm = [idx_rm, idx_compare];
            end
        end
    end
end
idx_rm = unique(idx_rm);

% Remove overlapping peaks
if ~isempty(idx_rm)
    flg_rm = true;
    peaks(idx_rm, :) = [];

    % Sort remaining peaks by amplitude (descending)
    [~, idx] = sort(peaks(:, 2), 'descend');
    peaks = peaks(idx, :);

    % Update model with new peaks
    fmdl.gaussian_params = peaks;

    % Convert gaussian parameters to peak parameters:
    % peak_params: [center frequency, power (area), bandwidth (width)]
    % gaussian_params: [center frequency, amplitude (height), sigma (std)]
    % Note: power is approximate area under gaussian curve
    fmdl.peak_params = [peaks(:, 1),...                     % center frequency
        peaks(:, 2) .* peaks(:, 3) * sqrt(2 * pi),...       % power (area)
        2 * peaks(:, 3)];                                   % bandwidth (FWHM)

    % Recalculate model components
    fmdl.fooofed_spectrum = peaks2spectrum(peaks, fmdl);
end

end

% create model fit (spectrum) from gaussian peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fooofed_spectrum] = peaks2spectrum(peaks, fmdl)

% Initialize periodic component
periodic = zeros(size(fmdl.freqs));

% Add each Gaussian peak
for ipeak = 1 : size(peaks, 1)
    cf = peaks(ipeak, 1);       % center frequency
    amp = peaks(ipeak, 2);      % amplitude
    sigma = peaks(ipeak, 3);    % width

    % Generate Gaussian
    gaussian = amp * exp(-((fmdl.freqs - cf) .^ 2 / (2 * sigma ^ 2)));
    periodic = periodic + gaussian;
end

% Combine components
fooofed_spectrum = fmdl.ap_fit + periodic;

end

% calculate model fit metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ffit = calc_fit(fmdl)

% Residuals
residuals = fmdl.power_spectrum - fmdl.fooofed_spectrum;

% Number of data points
n = length(fmdl.power_spectrum);

% Number of parameters (peaks + aperiodic)
k = length(fmdl.gaussian_params(:)) + size(fmdl.aperiodic_params, 2);

% Mean Squared Error of the model fit
ffit.mse = mean(residuals .^ 2);

% Log-likelihood of the model fit
ffit.loglik = -n / 2 * (log(2 * pi) + log(ffit.mse) + 1);

% Akaike Information Criterion
ffit.aic = 2 * k - 2 * ffit.loglik;

% Bayesian Information Criterion
ffit.bic = k * log(n) - 2 * ffit.loglik;

% r-squared
coefs = corrcoef(fmdl.power_spectrum, fmdl.fooofed_spectrum);
ffit.r_squared = coefs(1, 2);

% error (RMSE between original and model fit)
ffit.error = sqrt(mean((fmdl.power_spectrum - fmdl.fooofed_spectrum) .^ 2));

end

% organize output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psd_1of = fooof2params(fooof_cell)

nstates = length(fooof_cell);

for istate = 1 : nstates

    % convert to struct array
    pd = catfields(fooof_cell{istate}, 'addim', true);

    % reorganize peaks according to bands
    [pd.peak_params, pd.gaussian_params, pd.info] = peaks2bands(pd);
    pd.freqs = (pd.freqs(1, :, 1));

    % Separate parameters to their corresponding fields
    pd.cf = squeeze(pd.peak_params(:, 1, :));               % Center frequency [nbands x nbouts]
    pd.pow = squeeze(pd.peak_params(:, 2, :));              % Power [nbands x nbouts]
    pd.bw = squeeze(pd.peak_params(:, 3, :));               % Bandwidth [nbands x nbouts]
    pd.amp = squeeze(pd.gaussian_params(:, 2, :));          % Amplitude [nbands x nbouts]
    pd.sd = squeeze(pd.gaussian_params(:, 3, :));           % Standard deviation [nbands x nbouts]
    pd.ap_offset = squeeze(pd.aperiodic_params(1, 1, :))';
    pd.ap_exp = squeeze(pd.aperiodic_params(1, 2, :))';
    if size(pd.aperiodic_params, 2) == 2
        pd.ap_knee = [];
    else
        pd.ap_knee = squeeze(pd.aperiodic_params(1, 3, :))';
    end

    % rename fields
    pd.psd_orig = (pd.power_spectrum);
    pd.psd_fit = (pd.fooofed_spectrum);
    pd.psd_ap = (pd.ap_fit);

    % remove redundant fields
    pd = rmfield(pd, {'r_squared', 'error', 'power_spectrum', 'ap_fit',...
        'fooofed_spectrum', 'aperiodic_params', 'peak_params', 'gaussian_params'});

    p1of(istate) = pd;
end

psd_1of = catfields([p1of(:)], 'addim', true);
fit = psd_1of.fit;
fields = fieldnames(fit);
for ifield = 1:length(fields)
    val = fit.(fields{ifield});
    if ~isnumeric(val)
        continue
    end

    val = squeeze(val);
    dims = size(val);

    % If more than 1D, reorder so smallest dimension is first
    if length(dims) > 1
        [~, order] = sort(dims);
        val = permute(val, order);
    end

    % Assign back to struct
    fit.(fields{ifield}) = val;
end
psd_1of.fit = fit;
psd_1of.freqs = squeeze(psd_1of.freqs(1, :));
psd_1of.psd_orig = permute(squeeze(psd_1of.psd_orig), [3, 2, 1]);
psd_1of.psd_fit = permute(squeeze(psd_1of.psd_fit), [3, 2, 1]);
psd_1of.psd_ap = permute(squeeze(psd_1of.psd_ap), [3, 2, 1]);
psd_1of.cf = permute(psd_1of.cf, [3, 2, 1]);
psd_1of.pow = permute(psd_1of.pow, [3, 2, 1]);
psd_1of.bw = permute(psd_1of.bw, [3, 2, 1]);
psd_1of.sd = permute(psd_1of.sd, [3, 2, 1]);
psd_1of.amp = permute(psd_1of.amp, [3, 2, 1]);
psd_1of.ap_offset = squeeze(psd_1of.ap_offset);
psd_1of.ap_exp = squeeze(psd_1of.ap_exp);
psd_1of.ap_knee = squeeze(psd_1of.ap_knee);

% add porbability of detecting an oscillation in each band
[~, nbouts, nbands] = size(psd_1of.cf);
psd_1of.prob_osc = nan(nstates, nbands);
for istate = 1 : nstates   
    % for each band, count bouts with detected oscillation
    for iband = 1 : nbands
        psd_1of.prob_osc(istate, iband) = sum(~isnan(psd_1of.pow(istate, :, iband))) / nbouts;
    end
end

end

% restructure peaks according to pre-defined frequency bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pparams, gparams, info] = peaks2bands(fd)

% band params
bandNames = ["delta", "theta", "beta", "gamma"];
bandFreqs = [0, 5.5; 5.5, 12; 12, 30; 30, 100];
nbands = length(bandNames);

% initialize arrays
for iband = 1 : nbands
    pparams = nan(nbands, 3, size(fd.peak_params, 3));
    gparams = nan(nbands, 3, size(fd.gaussian_params, 3));
end

% Loop through bouts
for ibout = 1 : size(fd.peak_params, 3)
    peak_freqs = fd.peak_params(:, 1, ibout);
    peak_amps = fd.peak_params(:, 2, ibout);

    % For each band
    for iband = 1 : length(bandNames)
        % Find peaks in this band
        band_mask = peak_freqs >= bandFreqs(iband,1) & peak_freqs < bandFreqs(iband,2);

        if any(band_mask)

            % Get amplitudes of peaks in this band
            band_amps = peak_amps(band_mask);
            band_peaks = find(band_mask);

            % Find index of highest amplitude peak in this band
            [~, max_idx] = max(band_amps);
            peak_idx = band_peaks(max_idx);

            % Store parameters for highest amplitude peak
            pparams(iband, :, ibout) = fd.peak_params(peak_idx, :, ibout);
            gparams(iband, :, ibout) = fd.gaussian_params(peak_idx, :, ibout);
        end
    end
end

% append info
info.bandFreqs = bandFreqs;
info.bandNames = bandNames;

end

% interpolate psd to handle power line interference (PLI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psd_interp = interp_pli(psd_mat, freqs)

% PLI interpolation preserves ability to analyze full spectrum including gamma
% (30-100 Hz) while removing the confounding 50 Hz line noise artifact. This
% avoids edge effects from split frequency ranges and is more elegant than
% excluding peaks post-hoc. Particularly important when analyzing broadband
% phenomena like cross-frequency coupling. approach adapted from brainstorm.

freq_pli = [48 52];     % frequency range containing line noise [Hz]

% find indices for PLI
idx_pli = freqs >= freq_pli(1) & freqs <= freq_pli(2);
idx_clean = ~idx_pli;

% initialize interpolated psd
psd_interp = psd_mat;

% interpolate across PLI for each spectrum separately
for ispec = 1 : size(psd_mat, 1)
    % get values on either side of PLI
    x = freqs(idx_clean);
    y = psd_mat(ispec, idx_clean);

    % interpolate PLI region
    psd_interp(ispec, idx_pli) = interp1(x, y, freqs(idx_pli), 'pchip');
end

end

% plot psd and peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_fooof(fmdl)

% open figure
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [1, 2];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'FOOOF', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% original spectrum and fits
axh = nexttile(th, 1, [1, 1]); cla; hold on
hold on
plot(fmdl.freqs, fmdl.power_spectrum, 'k',...
    'LineWidth', 1.5, 'DisplayName', 'Original')
plot(fmdl.freqs, fmdl.fooofed_spectrum, 'r--',...
    'LineWidth', 1.5, 'DisplayName', 'FOOOF fit')
plot(fmdl.freqs, fmdl.ap_fit, 'b:',...
    'LineWidth', 1.5, 'DisplayName', 'Aperiodic fit')
xlabel('Frequency (Hz)')
ylabel('Power (log10)')
title(axh, sprintf('FOOOF Model Fit (R² = %.3f)', fmdl.r_squared))
legend('show')
grid on
set(axh', 'xscale', 'log')

% isolated peaks
axh = nexttile(th, 2, [1, 1]); cla; hold on

% plot residuals
residuals = fmdl.power_spectrum - fmdl.ap_fit;
plot(fmdl.freqs, residuals, 'k',...
    'LineWidth', 1, 'DisplayName', 'Residuals')

% plot detected peaks - no reshape needed
for ipeak = 1:size(fmdl.gaussian_params, 1)
    cf = fmdl.gaussian_params(ipeak, 1);       % center frequency
    sigma = fmdl.gaussian_params(ipeak, 3);    % width
    amp = fmdl.gaussian_params(ipeak, 2);      % amplitude

    % generate gaussian
    x = fmdl.freqs;
    gaussian = amp * exp(-((x - cf) .^ 2 / (2 * sigma ^ 2)));

    plot(x, gaussian, '--', 'LineWidth', 1.5,...
        'DisplayName', sprintf('Peak %.1f Hz', cf))
end

xlabel('Frequency (Hz)')
ylabel('Power (log10)')
title(axh, 'Isolated Peaks')
legend('show')
grid on
set(axh', 'xscale', 'log')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   fooof_results       struct array [nstates/nbouts x 1] with the following fields:
%
%   PERIODIC COMPONENTS (oscillations)
%   .peak_params        [n_peaks x 3] matrix of detected peaks with columns:
%                       1. center frequency (Hz)
%                       2. power (area under the Gaussian curve, not peak height)
%                       3. bandwidth (full width of the peak)
%   .gaussian_params    [n_peaks x 3] alternative parameterization of the same peaks:
%                       1. center frequency (Hz) - same as peak_params
%                       2. sigma (standard deviation of the Gaussian)
%                       3. amplitude (peak height)
%
%   APERIODIC COMPONENT (1/f background)
%   .aperiodic_params   [1 x 2] vector describing the 1/f background:
%                       1. offset (intercept in log-log space)
%                       2. exponent (slope in log-log space, steeper = more negative)
%                       Note: if aperiodic_mode='knee', returns [1 x 3] with
%                       additional knee parameter
%
%   MODEL FIT QUALITY
%   .error             goodness of fit (lower = better)
%                      measures difference between original and reconstructed spectrum
%   .r_squared         variance explained by model (higher = better)
%                      ranges from 0 to 1, with 1 being perfect fit
%
%   FULL MODEL COMPONENTS (if return_model=true)
%   .freqs             [1 x n_freqs] frequency vector
%   .power_spectrum    [1 x n_freqs] original power spectrum
%   .fooofed_spectrum  [1 x n_freqs] full model reconstruction
%                      combines aperiodic and periodic components
%   .ap_fit            [1 x n_freqs] isolated aperiodic component
%                      useful for calculating peak prominence above background
%
%   Note: Power values are in log10 units. To convert to linear power,
%   use: 10.^power_spectrum






% Another option for running fooof in matlab is through brainstorm
% opt.freq_range          = [1 100];
% opt.peak_width_limits   = [1 30];
% opt.max_peaks           = 6;
% opt.min_peak_height     = 1 / 10; % convert from dB to B
% opt.aperiodic_mode      = 'knee'; % 'knee'
% opt.peak_threshold      = 2;   % 2 std dev: parameter for interface simplification
% opt.border_threshold    = 1;   % 1 std dev: proximity to edge of spectrum, static in Python
% opt.return_spectrum     = 0;   % SPM/FT: set to 1
% opt.power_line          = '-5'; % otherwise 50, 60 Hz
% opt.proximity_threshold = 0.2;
% opt.optim_obj           = 'negloglike'; % negloglike
% opt.peak_type           = 'gaussian'; % 'cauchy', for interface simplification
% opt.guess_weight        = 'none'; % otherwise 'weak' or 'strong'
% opt.thresh_after        = true;   % Only when guess weight > 'None'
% opt.hOT                 = 1; % 0 if no optimization toolbox
%
% tf(:, 1, :) = log10(psd_data);
% [fs, fg] = ms_specparam(tf, freqs, opt, 1);
%
% % Select channel to visualize (e.g., channel 1)
% chan = 1;
%
% % Extract components
% original_spectrum = (squeeze(tf(1, chan, :)));  % Original spectrum
% modeled_spectrum = (fg(chan).fooofed_spectrum);  % Full modeled spectrum
% aperiodic_fit = (fg(chan).ap_fit);  % Aperiodic fit
% peak_fit = (fg(chan).peak_fit);  % Peaks only
%
% % Plot the results
% figure;
% hold on;
% plot(freqs, original_spectrum, 'k', 'LineWidth', 1.5, 'DisplayName', 'Original Spectrum');
% plot(fs, modeled_spectrum, 'r', 'LineWidth', 1.2, 'DisplayName', 'Modeled Spectrum');
% plot(fs, aperiodic_fit, 'b--', 'LineWidth', 1.2, 'DisplayName', 'Aperiodic Fit');
% plot(fs, peak_fit, 'g-.', 'LineWidth', 1.2, 'DisplayName', 'Peaks Only');
% legend show;
% r_squared = fg(chan).r_squared;
%     title(['Channel ', num2str(chan), ' - R^2 = ', num2str(r_squared)]);
% xlabel('Frequency (Hz)');
% ylabel('Power (log10)');
% title(['Spectral Fitting - Channel ', num2str(chan)]);
% grid on;
% hold off;
%
% % Extract peak parameters
% peak_params = fg(chan).peak_params;  % [Center, Height, Width]
% if isempty(peak_params)
%     disp('No peaks found for this channel.');
% else
%     % Overlay individual peaks on the plot
%     for i = 1:size(peak_params, 1)
%         peak_center = peak_params(i, 1);
%         peak_height = peak_params(i, 2);
%         peak_width = peak_params(i, 3);
%         % Generate the peak curve using Gaussian function
%         peak_curve = peak_height * exp(-(((freqs - peak_center) ./ peak_width) .^ 2) / 2);
%         hold on;
%         plot(freqs, peak_curve, '--', 'LineWidth', 1, 'DisplayName', ['Peak ', num2str(i)]);
%     end
%     legend show;
% end
