function frFit = mea_frFit(frVec, pertOnset, varargin)
% MEAFRFIT Fits double-exponential recovery model to firing rate data.
%
% SUMMARY:
% This function fits a double-exponential model to capture firing rate
% recovery dynamics after a perturbation. The model captures three phases:
%   1. Baseline period: Robust linear fit to pre-perturbation data
%   2. Perturbation period: Linear decline from baseline to trough
%   3. Recovery period: Double-exponential model with primary recovery
%      and optional overshoot components
%
% INPUT (Required):
%   frVec      - Vector of firing rate values [Hz] for a single unit.
%   pertOnset   - Index of perturbation onset.
%
% INPUT (Name-Value):
%   flgPlot     - Flag to plot results {false}.
%   baseMdls    - Cell array of base model names to include {'gompertz', 'richards', 'peak'}.
%                 Options: 'gomp', 'rich', 'exp1', 'exp2'
%
% OUTPUT:
%   frFit       - Structure containing fit results:
%     .frTrough       - Firing rate at recovery onset (trough) [Hz].
%     .rcvOnset       - Index of recovery onset based on best fitted curve.
%     .mdlName        - Name of the best-fitting model ('expo', 'gomp', 'rich').
%     .pFit           - Fitted parameters for the best model.
%     .rsquare        - R-squared goodness of fit for the best model.
%     .exitflag       - Exit flag of the best-fitting model.
%     .overshootTime  - Time to peak overshoot [s] (if applicable).
%     .overShootFr    - Peak overshoot firing rate [Hz] (if applicable).
%     .fitCurve       - Complete fitted curve [baseline; perturbation; recovery].
%     .frBsl          - Baseline firing rate from robust linear fit [Hz].
%     .bslFit         - Robust linear fit parameters [intercept, slope].
%     .frRcv          - Steady-state firing rate from best model [Hz].
%     .mdlFits        - Detailed results for each model fit.
%     .gof            - Goodness-of-fit stats (AIC, BIC) for each model.
%     .goodFit        - Flag indicating whether the unit is considered good for fitting.
%
% DEPENDENCIES:
%   Optimization Toolbox (for lsqcurvefit)
%   Statistics and Machine Learning Toolbox (for robustfit)
%
% HISTORY:
%   Sep 2024 - Extracted from mea_frRecovery.m as standalone function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'frVec', @isnumeric);
addRequired(p, 'pertOnset', @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'baseMdls', {'gomp', 'rich', 'exp2'}, @iscell);

parse(p, frVec, pertOnset, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE OUTPUT STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frFit = struct('frTrough', NaN, 'rcvOnset', NaN, 'mdlName', '', 'pFit', [], ...
    'rsquare', NaN, 'fitCurve', [], 'frBsl', NaN, 'bslFit', NaN, ...
    'frRcv', NaN, 'mdlFits', [], 'mdlIdx', [], 'goodFit', false, 'pertDepth', [], ...
    'rcvGain', NaN, 'rcvErr', NaN, 'fitQlty', struct('fitR2', false, ...
    'bslStab', false, 'rcvKin', false, 'rcvFail', false, 'pertResp', false, ...
    'bslRate', false, 'rcvRate', false));

% Hard limit on number of bins with spikes
if sum(frVec > 0) < 10
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUESS INITIAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert indices to relative time
tIdx = (1:length(frVec))';

% Initial guess for recovery onset using continuous derivative method
dfr = diff(frVec);
dfrPost = dfr(pertOnset:end);
runsUp = conv(dfrPost > 0, ones(1, 5), 'valid');
longUp = find(runsUp >= 5, 1, 'first');

if ~isempty(longUp)
    rcvOnset = pertOnset + longUp - 1;
else
    % Fallback: first non-negative derivative
    rcvDelay = find(dfrPost >= 0, 1, 'first');
    rcvOnset = pertOnset + rcvDelay - 1;
end

% Convert to relative time for fitting
tPost = tIdx(rcvOnset:end) - rcvOnset;
frPost = frVec(rcvOnset:end);

% Recovery firing rate: Mean of last 20% of data. 
frRcv = mean(frPost(round(0.8 * length(frPost)):end), 'omitnan');
frTrough = frVec(rcvOnset);
% Add lower bound for stability
frRcv = max([frRcv, eps]);
frTrough = max([frTrough, eps]);

% Recovery kinetics: Based on time to 50% recovery, or last bin if no clear recovery
rcvFr = frTrough + 0.5 * (frRcv - frTrough);
rcvIdx = find(frPost >= rcvFr, 1, 'first');
if ~isempty(rcvIdx)
    % Clear recovery: use time to 50% recovery, with upper bound for
    % stability
    tHalf = max([1, tPost(rcvIdx)]);
    kRcv = log(2) / tHalf;
else
    % No clear recovery: use rate from last bin. This will be slow if the
    % unit hasn't recovered.
    lastBinFr = frPost(end);
    kRcv = (lastBinFr - frTrough) / tPost(end);
end

% Hard coded limit to neurons with FR profiles that don't fit any model
if frTrough / frRcv > 2000
    warning('This unit does not make sense, returning...')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND PEAKS OF RECOVERY FOR OVERSHOOT MODELING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use findpeaks to detect nPeaks significant peaks. Require 10% prominence
minProminence = 0.2 * abs(frRcv - frTrough); 

[peakFr, peakTime] = findpeaks(frPost, tPost, ...
    'NPeaks', 1, ...
    'SortStr', 'descend', ...
    'MinPeakHeight', frRcv, ...
    'MinPeakProminence', minProminence, ...
    'MinPeakWidth', round(0.05 * length(tPost)), ...
    'MinPeakDistance', round(0.05 * length(tPost)));

% Approximate overshoot amplitude relative to recovered FR

% Calculate initial parameters for the overshoot term A*t*exp(-k*t)
% Peak of this term is at t = 1/k, value is A/(k*exp(1))
if peakTime > 0
    kOver = 1 / peakTime;
    kFall = 1 / (max(tPost) - peakTime);
    overshootAmp = peakFr - frRcv;
    aOver = overshootAmp * kOver * exp(1);
else
    kOver = 1 / 1800; % Default: characteristic time of 30 min
    kFall = 1 / 1800;
    aOver = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASE MODELS DESCRIPTION:
% All models use the same parameter order: [frRcv, k, frTrough, ...]
%   frRcv: steady-state firing rate (parameter 1)
%   k: rate parameter (parameter 2)
%   frTrough: trough value (parameter 3)
%
% Available base models:
%   - 'exp1': Exponential recovery model
%     y(t) = frRcv - (frRcv - frTrough) * exp(-k*t)
%
%   - 'gomp': Gompertz growth model
%     y(t) = frRcv * exp(log(frTrough/frRcv) * exp(-k*t))
%
%   - 'rich': Richards growth model
%     y(t) = frRcv / (1 + ((frRcv/frTrough)^a - 1) * exp(-k*t))^(1/a)
%     Parameters: [..., a] where a is shape parameter
%
%   - 'exp2': Peak model with rising and falling phases
%     Rising phase (t ≤ t_peak): frTrough + (peakFr - frTrough) * exp(k_rise * (t - t_peak))
%     Falling phase (t > t_peak): frRcv + (peakFr - frRcv) * exp(-k_fall * (t - t_peak))
%     Parameters: [frRcv, k_rise, frTrough, t_peak, peakFr, k_fall]

% Create model struct array for base models
mdls = struct('name', {}, 'func', {}, 'p0', {}, 'lb', {}, 'ub', {});

% Loop through specified base models
for iMdl = 1:length(p.Results.baseMdls)
    mdlName = p.Results.baseMdls{iMdl};
    
    switch lower(mdlName)
        case 'gomp'
            mdls(iMdl).name = "gomp";
            mdls(iMdl).func = @(p, tData) p(1) .* exp(log(max(p(3),eps) ./ p(1)) .* exp(-p(2) * tData));
            mdls(iMdl).p0 = [frRcv, kRcv, frTrough];
            mdls(iMdl).lb = [eps, 1e-6, eps];
            mdls(iMdl).ub = [3*max(frVec), 1, max(frVec)];
            
        case 'rich'
            mdls(iMdl).name = "rich";
            mdls(iMdl).func = @(p, tData) p(1) ./ (1 + ((p(1)./max(p(3),eps)).^p(4) - 1) .* exp(-p(2)*tData)).^(1./p(4));
            mdls(iMdl).p0 = [frRcv, kRcv, frTrough, 1];
            mdls(iMdl).lb = [eps, 1e-6, eps, 1e-3];
            mdls(iMdl).ub = [3*max(frVec), 1, max(frVec), 100];
            
        case 'exp1'
            mdls(iMdl).name = "exp1";
            mdls(iMdl).func = @(p, tData) p(1) - (p(1) - p(3)) .* exp(-p(2) * tData);
            mdls(iMdl).p0 = [frRcv, kRcv, frTrough];
            mdls(iMdl).lb = [eps, 1e-6, eps];
            mdls(iMdl).ub = [3*max(frVec), 1, max(frVec)];
            
        case 'exp2'
            mdls(iMdl).name = "exp2";
            mdls(iMdl).func = @(p, tData) ...
                (tData <= p(4)) .* ... % Rising phase (exponential approach to peak)
                (p(3) + (p(5)-p(3)) .* exp(p(2) * (tData - p(4)))) + ...
                (tData > p(4)) .* ...  % Falling phase (exponential decay from peak)
                (p(1) + (p(5)-p(1)) .* exp(-p(6) * (tData - p(4))));
            mdls(iMdl).p0 = [frRcv, kOver, frTrough, peakTime, peakFr, kFall];
            mdls(iMdl).lb = [eps, 1e-3, eps, eps, eps, 1e-3];
            mdls(iMdl).ub = [3*max(frVec), 20, max(frVec), max(tPost), 3*max(frVec), 20];
            
        otherwise
            warning('Unknown model type: %s. Skipping.', mdlName);
            continue;
    end
end

% If peaks were detected, create augmented models with overshoot terms
if peakTime > 0
    mdls = [mdls, mdls];
    for iMdl = length(mdls) / 2 + 1 : length(mdls)
        
        % Create new model with overshoot
        mdls(iMdl).name = mdls(iMdl).name + "+PK";
        
        % Store the base model function and number of parameters
        baseFunc = mdls(iMdl).func;
        numBaseParams = length(mdls(iMdl).p0);
        
        % Create simple overshoot model with 1 peak
        mdls(iMdl).func = @(p, tData) baseFunc(p(1:numBaseParams), tData) + ...
            p(numBaseParams + 1) * tData .* exp(-p(numBaseParams + 2) * tData);
        
        % Append parameters for overshoot terms
        mdls(iMdl).p0 = [mdls(iMdl).p0, aOver, kOver];
        mdls(iMdl).lb = [mdls(iMdl).lb, 0, 1e-6];
        mdls(iMdl).ub = [mdls(iMdl).ub, 3*max(frVec), 1];
    end
    
    % Remove exp2 model peak version from mdls struct
    idxExp2 = find(strcmpi([mdls.name], "exp2+PK"));
    if ~isempty(idxExp2)
        mdls(idxExp2) = [];
    end
else
    % Remove exp2 if no peak detected
    idxExp2 = find(strcmpi([mdls.name], "exp2"));
    if ~isempty(idxExp2)
        mdls(idxExp2) = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORM MODEL FITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nModels = length(mdls);

% Fit each model using the helper function
for iMdl = 1:nModels
    mdlFit(iMdl) = fit_mdl(frPost, tPost, mdls(iMdl));
end
frFit.mdlFits = catfields(mdlFit, 1, true);

% Check if any fit was successful
if all(frFit.mdlFits.exitflag <= 0)
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL SELECTION & GOODNESS-OF-FIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, frFit.mdlIdx] = min(frFit.mdlFits.BIC);
frFit.mdlName = mdls(frFit.mdlIdx).name;
frFit.pFit = frFit.mdlFits.pFit(frFit.mdlIdx, :);
frFit.rsquare = frFit.mdlFits.rsquare(frFit.mdlIdx);

% Update frTrough and frRcv based on the best model's fit
frFit.frTrough = frFit.pFit(3);
frFit.frRcv = frFit.pFit(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE BASELINE WITH ROBUST LINEAR FIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define baseline window (pre-perturbation onset) with margin
mrgnBins = 10;
bslWin = 1:max(1, pertOnset - mrgnBins);

% Fit robust linear regression to baseline period
[bslFit, ~] = robustfit(tIdx(bslWin), frVec(bslWin));

% Generate baseline curve using the linear fit
bslCurve = bslFit(1) + bslFit(2) * tIdx(1:(pertOnset-1));

% Use fitted line for calculating baseline FR
frFit.frBsl = mean(bslCurve);
frFit.bslFit = bslFit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE COMPLETE FITTED CURVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perturbation period: exponential decay
frTrough = max([frFit.frTrough, eps]);
tPert = tIdx(pertOnset:(rcvOnset-1));
pertCurve = linspace(bslCurve(end), frTrough, length(tPert))';

% Recovery period: fitted model
rcvCurve = frFit.mdlFits.rcvCurve(frFit.mdlIdx, :);

% Extend recovery curve to full data length if needed
if length(rcvCurve) < (length(frVec) - rcvOnset + 1)
    rcvCurve = mdls(frFit.mdlIdx).func(frFit.pFit, tIdx(rcvOnset:end) - rcvOnset);
end

% Combine segments
frFit.fitCurve = [bslCurve; pertCurve; rcvCurve'];

% Calculate recovery onset based on the best fitted curve
% Find the minimum of the fitted curve (trough)
% rcvCurve = frFit.fitCurve(pertOnset:end);
% [~, minIdx] = min(rcvCurve);
% frFit.rcvOnset = pertOnset + minIdx - 1;
frFit.rcvOnset = rcvOnset;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDENTIFY BAD UNITS BASED ON FIT QUALITY AND STABILITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section implements quality control to identify units unsuitable for recovery analysis.
% The criteria ensure reliable and interpretable results by checking fit quality, baseline
% stability, perturbation response, recovery kinetics, and minimum recovery rates.

% Initialize quality flags
fitQlty = struct();

% FIT QUALITY
% R-squared threshold ensures the model explains a meaningful portion of variance
thrRsquare = 0.1;
fitQlty.fitR2 = frFit.rsquare >= thrRsquare;
fitQlty.fitR2(isnan(frFit.rsquare)) = false; % Mark NaN R-squares as bad

% BASELINE STABILITY   
% Baseline slope should be minimal to ensure stable pre-perturbation activity
% Use normalized slope (% change per bin) to be independent of baseline FR
thrBslSlope = 0.2; 
if frFit.frBsl > 0
    bslSlope = abs(frFit.bslFit(2)) / frFit.frBsl; 
    fitQlty.bslStab = bslSlope <= thrBslSlope;
else
    fitQlty.bslStab = false;
end

% RECOVERY KINETICS
% Recovery time constant should be physiologically plausible
% Too fast (< 1 bin) suggests unreliable fitting or noise
kRcv = frFit.pFit(2);
tauRcv = 1 / kRcv; % Recovery time constant in bins
minTau = 1; 
fitQlty.rcvKin = tauRcv >= minTau;

% RECOVERY FAILURE
thrRcv = -4.6; % log2(0.01) for minimum ~1% recovery
frFit.rcvGain = log2(frFit.frRcv / frFit.frBsl);
frFit.rcvErr = abs(frFit.rcvGain);
fitQlty.rcvFail = frFit.rcvGain >= thrRcv;

% PERTURBATION RESPONSE
% Units should show meaningful perturbation effects (>X% reduction from baseline)
thrPert = -1; % log2(0.5) for ~50% reduction threshold
frFit.pertDepth = log2(frFit.frTrough / frFit.frBsl);
fitQlty.pertResp = frFit.pertDepth <= thrPert;

% BASELINE AND RECOVERY RATE
% Units should have sufficient activity for reliable analysis
thrFr = 0.001; % Hz
fitQlty.bslRate = frFit.frBsl >= thrFr;
fitQlty.rcvRate = frFit.frRcv >= thrFr;

% COMBINE 
qualityFields = fieldnames(fitQlty);
frFit.goodFit = true;
for iFld = 1:length(qualityFields)
    frFit.goodFit = frFit.goodFit && fitQlty.(qualityFields{iFld});
end

% OVVERIDE AND SELECT SPECIFIC CHECKS
frFit.goodFit = fitQlty.fitR2 & fitQlty.rcvKin & fitQlty.pertResp & ...
    fitQlty.bslRate;

% Store detailed quality information for debugging
frFit.fitQlty = fitQlty;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RESULTS (OPTIONAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p.Results.flgPlot
    plot_frFit_results(pertOnset, frFit, frVec, mdls);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: FIT_MDL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mdlFit = fit_mdl(fr, t, mdl)
% FIT_MDL Helper function to fit a model to firing rate data
%
% INPUT:
%   fr      - Firing rate data vector
%   t       - Time vector
%   mdl     - Model struct containing model function, initial parameters, and bounds
%
% OUTPUT:
%   mdlFit  - Structure containing fit results:
%     .pFit       - Fitted parameters
%     .resnorm    - Residual norm
%     .exitflag   - Exit flag from optimization
%     .rsquare    - R-squared goodness of fit
%     .AIC        - Akaike Information Criterion
%     .BIC        - Bayesian Information Criterion
%     .rcvCurve   - Recovery curve from fitted model

% Hardcoded optimization options
opts = optimoptions('lsqcurvefit', 'Display', 'off', 'TolFun', 1e-6,...
    'MaxIter', 1000, 'FunctionTolerance', 1e-6);

% Perform the fit with error handling
try
    [pFit, resnorm, ~, exitflag] = lsqcurvefit(mdl.func, mdl.p0, t, fr(:), mdl.lb, mdl.ub, opts);
catch ME
    % If optimization fails due to numerical issues, return failed fit
    if contains(ME.message, 'Inf') || contains(ME.message, 'NaN') || contains(ME.message, 'UndefAtX0')
        pFit = mdl.p0;  % Use initial parameters
        resnorm = inf;
        exitflag = -1;  % Indicate failure
    else
        % Re-throw other errors
        rethrow(ME);
    end
end

% Calculate R-squared if fit was successful
if exitflag > 0
    % Calculate fitted values
    frFit = mdl.func(pFit, t);
    
    % Calculate R-squared
    ssTot = sum((fr - mean(fr)).^2);
    if ssTot > 0
        rsquare = 1 - (resnorm / ssTot);
    else
        rsquare = 1; % If total variance is 0, and resnorm is 0, R2 is 1
    end
    
    % Calculate AIC and BIC
    nData = length(fr);
    k = length(pFit); % Number of parameters
    if resnorm > 0
        AIC = 2*k + nData * log(resnorm / nData);
        BIC = k * log(nData) + nData * log(resnorm / nData);
        if strcmp(mdl.name, 'exp2')
            penality = 14 * k * log(nData);
            BIC = BIC + penality;
        end
    else % resnorm is 0, perfect fit
        AIC = -inf;
        BIC = -inf;
    end
    
    % Create recovery curve
    rcvCurve = frFit;
else
    rsquare = NaN;
    AIC = inf;
    BIC = inf;
    rcvCurve = nan(size(t));
end

% Create output structure
mdlFit = struct('pFit', pFit, 'resnorm', resnorm, 'exitflag', exitflag,...
    'rsquare', rsquare, 'AIC', AIC, 'BIC', BIC, 'rcvCurve', rcvCurve');

% Transfer all fields from input model struct
fields = fieldnames(mdl);
for iFld = 1:length(fields)
    mdlFit.(fields{iFld}) = mdl.(fields{iFld});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: PLOT_FRFIT_RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_frFit_results(pertOnset, frFit, frVec, mdls)
% PLOT_FRFIT_RESULTS Helper function to plot firing rate fit results
%
% INPUT:
%   pertOnset   - Index of perturbation onset
%   frFit       - Structure containing fit results from mea_frFit
%   frVec       - Vector of firing rate values [Hz]
%   tIdx        - Time index vector
%   mdls        - Model struct array used for fitting

% Create figure with proper sizing
[hFig, hAx] = plot_axSize('szOnly', false, 'axShape', 'wide', 'axHeight', 300);
hold on;

% Hard code time indices 
tIdx = (1:length(frVec)) * 30 / 60 / 60 * 3;  
tLimit = [0, 24];
xlim(tLimit)
xticks([0 : 6 : 24])

% Plot raw data
plot(tIdx, frVec, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5,...
    'HandleVisibility', 'off');

% Get successful models and sort by BIC
mdlsGood = find(frFit.mdlFits.exitflag > 0);
[~, sortIdx] = sort(frFit.mdlFits.BIC(mdlsGood));
mdlsSrtd = mdlsGood(sortIdx);

% Plot each successful model fit for recovery period (sorted by BIC)
mdlNames = cell(1, length(mdlsSrtd));
clrMdls = bone(length(sortIdx)) * 0.6;
tRcvIdx = tIdx(frFit.rcvOnset:end);

for iMdl = 1:length(mdlsSrtd)
    mdlIdx = mdlsSrtd(iMdl);
    rcvCurvePlot = frFit.mdlFits.rcvCurve(mdlIdx, :);
    
    % Capitalize first letter of model name and add BIC
    modelName = mdls(mdlIdx).name{1};
    bicValue = round(frFit.mdlFits.BIC(mdlIdx));
    mdlNames{iMdl} = sprintf('%s (BIC: %d)', ...
        [upper(modelName(1)), modelName(2:end)], bicValue);

    hPlt = plot(tRcvIdx, rcvCurvePlot, '--', 'Color', [clrMdls(iMdl, :), 0.8], ...
        'LineWidth', 2, 'DisplayName', mdlNames{iMdl}, 'HandleVisibility', 'off');
    if iMdl == 1
        hPlt.HandleVisibility = 'on';
        hPlt.Color = [0, 0, 0, 1];
        hPlt.LineWidth = 2.5;
        hPlt.LineStyle = '-';
    end
end

% Mark perturbation
clrPert = [0.5, 0, 0.5, 0.5];
yLimit = ylim;
yText = yLimit(2);
yLine = yLimit(2) - 0.05 * (yLimit(2) - yLimit(1)); % Line 5% below text
plot([tIdx(pertOnset), 24], [yLine, yLine], 'Color', clrPert,...
    'LineWidth', 4, 'HandleVisibility', 'off');
text(tIdx(pertOnset) + 0.5, yText, 'Baclofen', 'FontName', 'Arial', 'FontSize', 16, ...
    'Color', clrPert, 'HorizontalAlignment', 'left');

% Formatting
xlabel('Time (bins)');
ylabel('Firing Rate (Hz)');   
legend('show', 'Location', 'northeast');   
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300);

end

