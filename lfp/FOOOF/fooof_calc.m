function f1f = fooof_calc(psd_data, freqs, varargin)
% fooof_calc - Parameterizes neural power spectra using FOOOF.
%
% This function leverages the FOOOF (Fitting Oscillations & One-Over F)
% algorithm to decompose neural power spectral densities (PSDs) into their
% fundamental components: aperiodic (1/f-like) activity and periodic
% oscillations (peaks). It utilizes a MATLAB wrapper for the Python
% implementation of FOOOF.
%
% Theoretical Considerations:
% Neural power spectra typically exhibit two main types of components:
% 1. Aperiodic Component: This is the 1/f-like characteristic, often
%    referred to as "background" activity. It reflects broadband,
%    non-oscillatory neural processes. FOOOF models this component using
%    the equation: P_ap(f) = L - log(k + f^X) for 'knee' mode, or
%    P_ap(f) = L - log(f^X) for 'fixed' mode (no knee).
%    - Offset (L): Broadband power shift of the spectrum.
%    - Knee (k): A parameter describing a bend in the aperiodic spectrum,
%      where the power law scaling changes. If k=0, it simplifies to the
%      'fixed' model.
%    - Exponent (X): The slope of the aperiodic component in log-log space,
%      thought to relate to the balance of excitation and inhibition in
%      the underlying neural population.
% 2. Periodic Components (Oscillations): These are deviations from the
%    aperiodic component, appearing as "bumps" or "peaks" in the spectrum.
%    They represent rhythmic neural activity in specific frequency bands.
%    FOOOF models these peaks as Gaussian functions, each characterized by:
%    - Center Frequency (CF): The frequency at which the peak is maximal.
%    - Power (PW): The power of the oscillation above the aperiodic component.
%    - Bandwidth (BW): The width of the peak.
%
% The function systematically evaluates models with a varying number of
% peaks (from 0 to a user-defined maximum). The Bayesian Information
% Criterion (BIC) is used to select the most parsimonious model that best
% fits the data, penalizing for model complexity. A Bayes Factor (BF) is
% then computed to compare the best oscillatory model against a purely
% aperiodic model (0 peaks). If BF > 1, the aperiodic model is favored,
% suggesting the absence of significant periodic activity.
%
% By separating these components, fooof_calc allows for a more nuanced
% analysis of neural data, enabling the study of oscillatory dynamics
% independent of broadband shifts in power and providing insights into
% the characteristics of the aperiodic background.
%
% INPUT
%   psd_data    [nBouts x nFreqs] double: Power spectral density data for
%               each bout (row) across frequencies (columns).
%   freqs       [1 x nFreqs] double: Vector of frequencies corresponding
%               to the columns of psd_data.
%   varargin    Optional input arguments (name-value pairs):
%     'basepath'  char: Path to the recording session. {pwd}
%     'f_range'   [1x2] double: Frequency range for fitting [min_freq, max_freq].
%                 Defaults to the full range of 'freqs'.
%     'fooof_cfg' struct: Configuration structure for the FOOOF algorithm.
%                 If empty or not provided, default settings are used.
%                 See "Model Parameters" section within the code for defaults.
%                 Key fields include:
%                   'peak_width_limits': [min max] width for peaks.
%                   'max_n_peaks': Maximum number of peaks to fit.
%                   'min_peak_height': Minimum height of a peak.
%                   'peak_threshold': Threshold for detecting peaks.
%                   'aperiodic_mode': 'fixed' or 'knee'.
%                   'verbose': Print FOOOF outputs.
%     'flg_plot'  logical: If true, plots FOOOF results for each bout. {true}
%     'saveVar'   logical | char: If true, saves the output 'f1f' struct.
%                 If a char, it's used as part of the filename. {true}
%
% OUTPUT
%   f1f         struct: Contains the FOOOF results and related information.
%     .psd_orig      [nBouts x nFreqs] - Original power spectrum (linear scale).
%     .psd_fit       [nBouts x nFreqs] - Full FOOOF model fit (aperiodic + periodic, linear scale).
%     .psd_ap        [nBouts x nFreqs] - Aperiodic component of the spectrum (linear scale).
%     .psd_res       [nBouts x nFreqs] - Residual spectrum (psd_orig - psd_ap, linear scale).
%     .psd_res_fit   [nBouts x nFreqs] - Residual of model fit (psd_fit - psd_ap, linear scale).
%     .ap            struct with aperiodic parameters:
%       .offset      [nBouts x 1] - Aperiodic offset.
%       .exp         [nBouts x 1] - Aperiodic exponent.
%       .knee        [nBouts x 1] - Aperiodic knee (NaN if not present or not modeled).
%     .peaks         struct with peak parameters organized by band:
%       .cf          [nBouts x nBands] - Center frequencies of selected peaks.
%       .pow         [nBouts x nBands] - Power of selected peaks.
%       .bw          [nBouts x nBands] - Bandwidth of selected peaks.
%       .amp         [nBouts x nBands] - Amplitude of Gaussian fit for peaks.
%       .sd          [nBouts x nBands] - Standard deviation of Gaussian fit for peaks.
%       .prob        [1 x nBands]    - Probability of detecting an oscillation in each band
%                                      (proportion of bouts with a peak in that band).
%       .info        struct: Contains band definitions (.bandNames, .bandFreqs).
%     .bandPow       [nBouts x nBands] - Sum of power in the *residual* spectrum (psd_res)
%                                      within predefined frequency bands.
%     .freqs         [1 x nFreqs_fit] - Frequency vector corresponding to the fitted range.
%     .fit           struct with goodness-of-fit metrics for the selected model:
%       .rmse        [nBouts x 1] - Root mean squared error.
%       .loglik      [nBouts x 1] - Log-likelihood.
%       .aic         [nBouts x 1] - Akaike Information Criterion.
%       .bic         [nBouts x 1] - Bayesian Information Criterion.
%       .r_squared   [nBouts x 1] - R-squared value of the fit.
%     .info          struct with analysis details:
%       .fooof_cfg   struct - FOOOF configuration used.
%       .runtime     datetime - Timestamp of when the analysis was run.
%       .cfg         (duplicate of .fooof_cfg, consider removing for clarity)
%
% CALLS
%   fooof (Python FOOOF via MATLAB wrapper)
%   peaks2bands (local function)
%   calc_fit (local function)
%   catfields (custom utility)
%   psd_fooofPlot (custom utility)
%
% DEPENDENCIES
%   MATLAB R2019b or later (for Python integration).
%   Python (e.g., 3.8-3.12) with the following packages installed:
%     numpy, fooof, matplotlib
%   To install Python packages: `pip install numpy fooof matplotlib` or
%   `py -m pip install numpy fooof matplotlib` (Windows example).
%
% SEE ALSO
%   Official FOOOF documentation: https://fooof-tools.github.io/fooof/
%
% VERSION HISTORY
%   09 Jan 24 LH - Initial version
%   [Current Date] [Your Initials] - Documentation enhanced.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 1: Input Parsing & Initialization
% This section handles the optional input arguments provided to the function,
% sets default values if arguments are not specified, and prepares basic
% parameters like the output filename and ensures the frequency vector is in
% the correct orientation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'fooof_cfg', []);
addOptional(p, 'basepath', pwd);
addOptional(p, 'f_range', []);
addOptional(p, 'flg_plot', true, @islogical);
addOptional(p, 'saveVar', true);

parse(p, varargin{:})
fooof_cfg           = p.Results.fooof_cfg;
f_range             = p.Results.f_range;
basepath            = p.Results.basepath;
flg_plot            = p.Results.flg_plot;
saveVar             = p.Results.saveVar;

if isempty(f_range)
    f_range = [min(freqs), max(freqs)];
end

% ensure row vectors
freqs = freqs(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 2: Filename Generation and Path Handling
% Determines the full path and filename for saving the output .mat file
% based on the 'basepath' and 'saveVar' inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, basename] = fileparts(basepath);
if ischar(saveVar)
    fname = [basepath, filesep, basename, '.' saveVar '.mat'];
else
    fname = fullfile(basepath, [basename, '.psd_1of.mat']);
end

% default band params
bands.names = ["delta", "theta", "gamma"]; 
bands.freqs = [2, 5; 6, 12; 30, 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 3: Python Environment and FOOOF Initialization
% This section ensures that the necessary Python environment is configured and
% that the required Python modules (numpy, fooof, matplotlib) are imported
% and accessible from MATLAB. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try         % check if modules already imported
    py.importlib.import_module('numpy');
    py.importlib.import_module('fooof');
    py.importlib.import_module('matplotlib');
catch ME
    clear classes % Clears MATLAB class definitions, which can sometimes help resolve Python integration issues.
    clear py.*    % Clears Python-related variables in the MATLAB workspace.
    pyenv;        % Initializes or re-initializes the Python environment MATLAB is using.
                  % Ensure this environment has fooof, numpy, matplotlib installed.
    pyversion     % Displays the Python version MATLAB is currently configured to use.
    % Verify FOOOF is installed in that Python environment
    py.help('fooof') % Attempts to access the help documentation for the fooof Python module.
                     % This serves as a check that the module is callable.
    py.importlib.import_module('numpy');
    py.importlib.import_module('fooof');
    py.importlib.import_module('matplotlib');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 4: FOOOF Model Configuration
% This section defines the parameters for the FOOOF model. If a custom
% configuration ('fooof_cfg') is not provided via input arguments,
% default settings are applied. These settings control aspects of the peak
% fitting process (e.g., peak width limits, maximum number of peaks) and
% the aperiodic fit (e.g., 'knee' vs 'fixed' mode).
% The selected configuration is stored in 'f1f.info.fooof_cfg' for record-keeping.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(fooof_cfg)
    % Default FOOOF settings if not specified by the user.
    % 'peak_width_limits': [min max] in Hz, defining the allowed range for the
    %                      bandwidth of detected peaks.
    % 'max_n_peaks': The upper limit on the number of peaks the model will attempt to fit.
    %                The model selection process (BIC) will determine the optimal number up to this limit.
    % 'min_peak_height': Minimum absolute height of a peak (in log10 power units) above the
    %                    aperiodic fit to be considered.
    % 'peak_threshold': A relative threshold (in standard deviations of the background noise)
    %                   for a peak's height relative to the aperiodic component.
    % 'aperiodic_mode': Specifies the model for the aperiodic component.
    %                   'knee': Fits the aperiodic component with an offset, knee, and exponent.
    %                   'fixed': Fits without a knee parameter (offset and exponent only).
    % 'verbose': If true, FOOOF will print fitting information to the console.
    % 'return_model': If true, the full FOOOF model object is returned by the Python call.

    fooof_cfg = struct(...
        'peak_width_limits', [1 12],...   % Wider peaks allowed by default
        'max_n_peaks', 6,...                % Max 6 peaks
        'min_peak_height', 0,...            % Min height in log10(power)
        'peak_threshold', 1,...             % Relative height threshold
        'aperiodic_mode', 'fixed',...        % Use knee fit for aperiodic
        'verbose', true,...
        'return_model', true);
end

% append info to output
f1f.info.fooof_cfg = fooof_cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 5: Model Selection Strategy - Iterative Fitting and BIC
% This section outlines the strategy for selecting the best FOOOF model.
% Instead of fitting a single model with 'max_n_peaks', the code iterates
% through models with a varying number of peaks, from 0 (aperiodic only)
% up to the 'max_n_peaks' defined in 'fooof_cfg'.
% For each number of peaks, a FOOOF model is fitted. The Bayesian
% Information Criterion (BIC) is calculated for each of these models.
% The model yielding the lowest BIC is initially considered the "best"
% as BIC penalizes model complexity, favoring simpler models that explain
% the data well. Subsequently, a Bayes Factor is used to compare this
% "best" oscillatory model against the purely aperiodic model (0 peaks)
% to ensure that the inclusion of peaks is justified.
%
% This approach is adapted from the method described in:
% https://www.biorxiv.org/content/10.1101/2024.08.01.606216v1.full.pdf
% (see supplementary materials for Brainstorm implementation).
% The rationale is to provide a more robust selection of periodic features.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here, a model will be created for each combination of parameters below.
npeaks = [0 : fooof_cfg.max_n_peaks];
nmdls = length(npeaks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 6: FOOOF Execution and Model Evaluation per Bout
% This is the main processing loop. It iterates through each 'bout' (epoch or
% segment of PSD data). For each bout:
% 1. It iterates through the different numbers of peaks specified in 'npeaks'
%    (from 0 to 'max_n_peaks').
% 2. For each number of peaks, it adjusts 'fooof_cfg.max_n_peaks' and calls
%    the 'fooof' Python function via the MATLAB wrapper. The FOOOF output
%    (spectra and parameters) is in log10 power units.
% 3. Goodness-of-fit metrics (RMSE, log-likelihood, AIC, BIC, R^2) are
%    calculated for the current model using the local function 'calc_fit'.
% 4. After all models (0 to 'max_n_peaks') are fitted for the current bout,
%    the model with the minimum BIC is identified as the provisionally best model.
% 5. A Bayes Factor (BF) is calculated: BF = exp((BIC_best_osc - BIC_aperiodic) / 2).
%    If BF > 1, it indicates that the aperiodic model (0 peaks) is more likely
%    given the data than the best model with oscillations. In this case, the
%    aperiodic model is chosen as the final model for the bout.
% 6. The selected FOOOF model results and its fit metrics for the current
%    bout are stored.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get number of bouts
[nBouts, nfreq] = size(psd_data); % Directly use psd_data

% validate psd input
if nfreq ~= length(freqs)
    error('chech psd_data input')
end

% initialize output
clear ffMdl
% loop through bouts
for ibout = 1 : nBouts
    for imdl = 1 : nmdls

        % adjust model parameters
        fooof_cfg.max_n_peaks = npeaks(imdl);

        % run fooof. note output (and subsequent spectrums) are in power
        % (log10)
        ffMdl(imdl) = fooof(freqs, psd_data(ibout, :),... 
            f_range, fooof_cfg, fooof_cfg.return_model);

        % calculate model fit
        mdl_fit(imdl) = calc_fit(ffMdl(imdl));
    end

    % select best model
    bic = vertcat(mdl_fit.bic);
    [~, idx_best] = min(bic);

    % Bayes Factor between best model and apreiodic component
    bf = exp((bic(idx_best) - bic(1)) ./ 2);
    if bf > 1
        idx_best = 1;
    end
    
    ffTmp = ffMdl(idx_best);
    ffTmp.fit = mdl_fit(idx_best);
    
    % append
    ffIn(ibout) = ffTmp;
end

% concatenate
ffIn = catfields(ffIn, 'addim');

% organize output 
f1f = fooof_org(ffIn, bands); 
mask_freq = freqs >= f_range(1) & freqs <= f_range(2);
f1f.freqs = freqs(mask_freq);

% plot
if flg_plot
    fh = fooof_plotMdl(f1f);
end

% organize info
f1f.info.runtime = datetime("now");
f1f.info.cfg = fooof_cfg;

% save
if saveVar
    save(fname, 'f1f')
end

end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc_fit: Calculate Model Fit Metrics
%
% This local function computes various goodness-of-fit statistics for a
% given FOOOF model output ('fmdl'). These metrics help evaluate how well
% the model (aperiodic + periodic components) describes the original power
% spectrum.
%
% INPUT
%   fmdl    struct: A single FOOOF model output structure, typically for one
%           bout and one specific model configuration. It must contain:
%             .power_spectrum      - The original power spectrum (log10 units).
%             .fooofed_spectrum    - The full model fit (log10 units).
%             .gaussian_params     - Parameters of the Gaussian peaks.
%             .aperiodic_params    - Parameters of the aperiodic fit.
%
% OUTPUT
%   ffit    struct: Contains the calculated fit metrics:
%     .rmse        Root Mean Squared Error between original and model spectrum.
%     .loglik      Log-likelihood of the model.
%     .aic         Akaike Information Criterion.
%     .bic         Bayesian Information Criterion.
%     .r_squared   Coefficient of determination (R^2) between original
%                  and model spectrum.
%
% Theory:
% - RMSE: Measures the average difference between observed and predicted values.
% - Log-likelihood: Quantifies the likelihood of observing the data given the model.
%   Higher is better.
% - AIC (Akaike Information Criterion): Balances model fit (likelihood) against
%   model complexity (number of parameters). Lower AIC is preferred.
%   AIC = 2k - 2*log(L), where k is number of parameters, L is likelihood.
% - BIC (Bayesian Information Criterion): Similar to AIC but with a stronger
%   penalty for model complexity, especially for larger datasets. Lower BIC
%   is preferred. BIC = k*log(n) - 2*log(L), where n is number of data points.
% - R-squared: Proportion of variance in the original spectrum explained by
%   the model. Ranges from 0 to 1 (or -inf to 1), higher is better.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ffit = calc_fit(fmdl)

% Number of data points
n = length(fmdl.power_spectrum);

% Number of parameters (peaks + aperiodic)
k = length(fmdl.gaussian_params(:)) + size(fmdl.aperiodic_params, 2);

% Root Mean Squared Error of the model fit (difference between original and model fit)
ffit.rmse = sqrt(mean((fmdl.power_spectrum - fmdl.fooofed_spectrum) .^ 2));

% Log-likelihood of the model fit
ffit.loglik = -n / 2 * (log(2 * pi) + log(ffit.rmse) + 1);

% Akaike Information Criterion
ffit.aic = 2 * k - 2 * ffit.loglik;

% Bayesian Information Criterion
ffit.bic = k * log(n) - 2 * ffit.loglik;

% r-squared
coefs = corrcoef(fmdl.power_spectrum, fmdl.fooofed_spectrum);
ffit.r_squared = coefs(1, 2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fooof_org: Organize FOOOF Output Structure
%
% This local function takes the concatenated raw output from multiple FOOOF
% model fits (typically after 'catfields') and reorganizes it into a more
% user-friendly and standardized structure ('ffOut').
% Key operations include:
% 1. Extracting and reshaping spectra (original, full fit, aperiodic component).
%    Converts spectra from log10 power (as returned by FOOOF) to linear power.
% 2. Calculating residual spectra (original - aperiodic, full fit - aperiodic).
% 3. Extracting aperiodic parameters (offset, knee, exponent), handling cases
%    with and without the knee parameter.
% 4. Calling 'peaks2bands' to process and assign peak parameters (center
%    frequency, power, bandwidth, Gaussian amplitude, Gaussian standard
%    deviation) to predefined frequency bands.
% 5. Calculating the probability of detecting an oscillation within each
%    predefined band across all processed bouts.
%
% INPUT
%   ffIn    struct: A structure containing concatenated fields from FOOOF
%           outputs across multiple bouts. Expected fields include:
%             .freqs             - Frequency vector (from FOOOF model).
%             .power_spectrum    - Original power spectra (log10 units).
%             .fooofed_spectrum  - Full model fits (log10 units).
%             .ap_fit            - Aperiodic component fits (log10 units).
%             .aperiodic_params  - Aperiodic parameters.
%             .peak_params       - Peak parameters from FOOOF.
%             .gaussian_params   - Gaussian parameters for peaks from FOOOF.
%
% OUTPUT
%   ffOut   struct: An organized structure containing:
%     .freqs         [1 x nFreqs]       - Frequency vector.
%     .psd_orig      [nBouts x nFreqs]  - Original PSD (linear).
%     .psd_fit       [nBouts x nFreqs]  - Full model fit (linear).
%     .psd_ap        [nBouts x nFreqs]  - Aperiodic component (linear).
%     .psd_res       [nBouts x nFreqs]  - Residual (original - aperiodic, linear).
%     .psd_res_fit   [nBouts x nFreqs]  - Residual of model fit (fit - aperiodic, linear).
%     .ap            struct:
%       .offset      [nBouts x 1]
%       .knee        [nBouts x 1] (NaN if not present)
%       .exp         [nBouts x 1]
%     .peaks         struct: (populated by 'peaks2bands')
%       .cf          [nBouts x nBands]
%       .pow         [nBouts x nBands]
%       .bw          [nBouts x nBands]
%       .amp         [nBouts x nBands]
%       .sd          [nBouts x nBands]
%       .prob        [1 x nBands]
%       .info        struct (band definitions from 'peaks2bands')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ffOut = fooof_org(ffIn, bands) 

ffOut.freqs = squeeze(ffIn.freqs(1, :, 1)); 
nBouts = length(ffIn.aperiodic_params(1, 1, :));
nFreqs = length(ffOut.freqs);

% Convert from power to linear
ffOut.psd_orig = reshape([10 .^ squeeze(ffIn.power_spectrum)]', [nBouts, nFreqs]);
ffOut.psd_fit = reshape([10 .^ squeeze(ffIn.fooofed_spectrum)]', [nBouts, nFreqs]);
ffOut.psd_ap = reshape([10 .^ squeeze(ffIn.ap_fit)]', [nBouts, nFreqs]);
ffOut.psd_res = ffOut.psd_orig - ffOut.psd_ap;
ffOut.psd_res_fit = ffOut.psd_fit - ffOut.psd_ap;

% Aperiodic parameters - can be [1x2] or [1x3]
ffOut.ap.offset = squeeze(ffIn.aperiodic_params(:, 1, :));
if size(ffIn.aperiodic_params, 2) == 2
    ffOut.ap.knee = nan(nBouts, 1); 
    ffOut.ap.exp = squeeze(ffIn.aperiodic_params(:, 2, :));
else
    ffOut.ap.knee = squeeze(ffIn.aperiodic_params(:, 2, :));
    ffOut.ap.exp = squeeze(ffIn.aperiodic_params(:, 3, :));
end

% Peak parameters need to be processed by peaks2bands
% peaks2bands will take the struct array and return [nbands x 3 x nBouts]
[pParams, gParams] = peaks2bands(ffIn, bands);

% peak_params_banded is [nbands x 3 x nBouts]
% We want cf, pow, bw as [nBouts x nbands]
ffOut.peaks.cf = squeeze(pParams(:, 1, :))';            % Center frequency 
ffOut.peaks.pow = squeeze(pParams(:, 2, :))';           % Power
ffOut.peaks.bw = squeeze(pParams(:, 3, :))';            % Bandwidth

% gaussian_params_banded is [nbands x 3 x nBouts]
% We want amp, sd as [nBouts x nbands]
ffOut.peaks.amp = squeeze(gParams(:, 2, :))';           % Amplitude
ffOut.peaks.sd = squeeze(gParams(:, 3, :))';            % Standard deviation

% add probability of detecting an oscillation in each band
nBands = length(bands.names); 
ffOut.peaks.prob = nan(1, nBands); 
for iband = 1 : nBands
    ffOut.peaks.prob(1, iband) = sum(~isnan(ffOut.peaks.pow(:, iband))) / nBouts;
end

% To quantify oscillatory activity within predefined frequency bands after
% accounting for the aperiodic component, two distinct metrics were derived
% from the original power spectral density (PSD(f)) and the fitted
% aperiodic component (AP(f)). First, to capture the total absolute power
% of periodic activity, Oscillatory Band Energy was calculated by summing
% the linear residuals within each band: sum(PSD(f)−AP(f)). This metric
% directly reflects the absolute energetic contribution of oscillations,
% which may inherently give more weight to activity at lower frequencies
% where baseline power is often higher. Second, to assess oscillatory
% strength relative to the varying aperiodic background and thus provide a
% measure less dominated by absolute power levels, a Logarithmic Band Index
% was computed as the sum of the differences between their respective
% log10-transformed values: sum(log10(PSD(f)/(AP(f)). This latter index
% emphasizes the consistency and relative magnitude of deviations from the
% aperiodic fit across the band.
resLog = log10(ffOut.psd_orig) - log10(ffOut.psd_ap);
ffOut.bands.pow = nan(nBouts, nBands);
ffOut.bands.energy = nan(nBouts, nBands);
for iband = 1 : nBands
    mask_band = ffOut.freqs >= bands.freqs(iband, 1) & ffOut.freqs < bands.freqs(iband, 2);
    ffOut.bands.energy(:, iband) = sum(ffOut.psd_res(:, mask_band), 2);
    ffOut.bands.pow(:, iband) = sum(resLog(:, mask_band), 2);
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peaks2bands: Restructure Peak Parameters into Frequency Bands
%
% This local function takes the raw peak parameters and Gaussian parameters
% from the FOOOF output (for all bouts, potentially after concatenation by
% 'catfields') and organizes them according to predefined frequency bands.
%
% For each bout and for each predefined frequency band:
% 1. It identifies all peaks from the FOOOF output whose center frequencies
%    fall within that band.
% 2. If multiple peaks are found in a band for a given bout, it selects the
%    peak with the highest power.
% 3. The parameters (Center Frequency, Power, Bandwidth from 'peak_params';
%    Amplitude, Standard Deviation from 'gaussian_params') of this selected
%    peak are stored for that band and bout. If no peak is found in a band,
%    NaN values are stored.
%
% INPUT
%   ffIn    struct: Contains concatenated FOOOF outputs. Relevant fields:
%             .peak_params      [max_n_peaks x 3 x nBouts] - CF, Power, BW
%             .gaussian_params  [max_n_peaks x 3 x nBouts] - CF, Amp, StdDev
%                                (Note: Gaussian CF is redundant with peak_params CF)
%   nBouts  scalar: The number of bouts, passed from 'fooof_org'.
%
% OUTPUT
%   pParams [nBands x 3 x nBouts]: Peak parameters (CF, Power, BW) for the
%            dominant peak in each band for each bout.
%   gParams [nBands x 3 x nBouts]: Gaussian parameters (CF, Amplitude, StdDev)
%            for the dominant peak in each band for each bout.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pParams, gParams] = peaks2bands(ffIn, bands) 

% params
nbands = length(bands.names);
nBouts = length(ffIn.aperiodic_params(1, 1, :));

% Initialize output arrays: [nbands x 3 params x nBouts]
pParams = nan(nbands, 3, nBouts);
gParams = nan(nbands, 3, nBouts);

% Loop through bouts (elements of the struct array fd)
for ibout = 1 : nBouts
    pBout = ffIn.peak_params(:, :, ibout);            % [nPeaks x 3]
    gBout = ffIn.gaussian_params(:, :, ibout);        % [nPeaks x 3]
    
    if isempty(pBout) || all(isnan(pBout(:)))
        continue; % Skip if no peaks for this bout
    end

    peak_freqs = pBout(:, 1);
    peak_powers = pBout(:, 2);          

    % For each band
    for iband = 1 : nbands
        
        % Find peaks in this band
        band_mask = peak_freqs >= bands.freqs(iband,1) &...
            peak_freqs < bands.freqs(iband,2);

        if any(band_mask)
            
            % Get power of peaks in this band
            powBand = peak_powers(band_mask);
            idxBand = find(band_mask);

            % Find index of highest power peak in this band (relative to selection)
            [~, idxMax] = max(powBand);
            
            % Get the original index from current_bout_peak_params
            idxPeak = idxBand(idxMax);

            % Store parameters for highest power peak
            pParams(iband, :, ibout) = pBout(idxPeak, :);
            gParams(iband, :, ibout) = gBout(idxPeak, :);
        end
    end
end
end


