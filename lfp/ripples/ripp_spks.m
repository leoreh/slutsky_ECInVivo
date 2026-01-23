function rippSpks = ripp_spks(spkTimes, rippTimes, ctrlTimes, varargin)
% RIPP_SPKS Analyzes spiking rate modulation during ripples.
%
%   rippSpks = RIPP_SPKS(spkTimes, rippTimes, ctrlTimes, varargin)
%
%   SUMMARY:
%       Calculates scalar modulation metrics comparing Ripple vs Control periods.
%       1. Instantaneous Firing Rates (FR) per event.
%       2. Mean FR differences (RippleMu vs ControlMu).
%       3. Z-Scored Gain.
%       4. Statistical Significance (Wilcoxon Sign-Rank / Rank-Sum).
%
%   INPUTS:
%       spkTimes    - (Cell) {N_units x 1} Spike times [s].
%       rippTimes   - (Mat)  [N x 2] Ripple start/end times [s].
%       ctrlTimes   - (Mat)  [N x 2] Control start/end times [s].
%       varargin    - Parameter/Value pairs:
%           'basepath' - (Char) Save location. (Default: pwd).
%           'flgSave'  - (Log)  Save output? (Default: true).
%
%   OUTPUTS:
%       rippSpks    - (Struct) Stats structure with [N_units x 1] fields:
%           .frRipp    - Mean Firing Rate during Ripples (Hz).
%           .frRand    - Mean Firing Rate during Control (Hz).
%           .frZ       - Z-scored modulation ((Ripp - Ctrl) / StdCtrl).
%           .frMod     - Modulation Index ((Ripp - Ctrl) / (Ripp + Ctrl)).
%           .pVal      - P-value from significance test.
%           .sigMod    - Boolean significance flag (p < 0.05).
%
%   DEPENDENCIES:
%       times2rate.
%
%   HISTORY:
%       Updated: 23 Jan 2026

% =========================================================================
%  ARGUMENTS
% =========================================================================
p = inputParser;
addRequired(p, 'spkTimes', @iscell);
addRequired(p, 'rippTimes', @isnumeric);
addRequired(p, 'ctrlTimes', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
parse(p, spkTimes, rippTimes, ctrlTimes, varargin{:});

basepath = p.Results.basepath;
flgSave = p.Results.flgSave;

% =========================================================================
%  CALCULATION
% =========================================================================
[~, basename] = fileparts(basepath);
savefile = fullfile(basepath, [basename, '.rippSpks.mat']);

nUnits = length(spkTimes);

% Initialize Output
rippSpks = struct();
rippSpks.pVal = nan(nUnits, 1);
rippSpks.sigMod = false(nUnits, 1);
rippSpks.frZ = nan(nUnits, 1);
rippSpks.frMod = nan(nUnits, 1);
rippSpks.frRipp = nan(nUnits, 1);
rippSpks.frRand = nan(nUnits, 1);

% Calculate Instantaneous Rates (Events x Units)
% 'binsize', Inf -> Returns one rate per event
rippSpks.rippRates = times2rate(spkTimes, 'winCalc', rippTimes, 'binsize', Inf);
rippSpks.ctrlRates = times2rate(spkTimes, 'winCalc', ctrlTimes, 'binsize', Inf);

% Per Unit Statistics
for iUnit = 1:nUnits

    rRate = rippSpks.rippRates(iUnit, :);
    cRate = rippSpks.ctrlRates(iUnit, :);

    if any(isnan(rRate)) || any(isnan(cRate))
        continue;
    end

    % Means & Std
    rMu = mean(rRate, 'omitnan');
    cMu = mean(cRate, 'omitnan');
    cSigma = std(cRate, 'omitnan');

    % Store Descriptive Stats
    rippSpks.frRipp(iUnit) = rMu;
    rippSpks.frRand(iUnit) = cMu;

    % Metrics
    if cSigma > 0
        rippSpks.frZ(iUnit) = (rMu - cMu) / cSigma;
    else
        rippSpks.frZ(iUnit) = NaN; % Avoid infinite gain
    end

    % Modulation Index (-1 to 1)
    if (rMu + cMu) > 0
        rippSpks.frMod(iUnit) = (rMu - cMu) / (rMu + cMu);
    end

    % Statistical Significance (Wilcoxon Sign Rank)
    % Paired test comparing distribution of ripple rates vs control rates
    % Note: Requires equal number of events. If N_ripp != N_ctrl,
    % we use ranksum (unpaired)
    if length(rRate) == length(cRate)
        rippSpks.pVal(iUnit) = signrank(rRate, cRate);
    else
        rippSpks.pVal(iUnit) = ranksum(rRate, cRate);
    end

    rippSpks.sigMod(iUnit) = rippSpks.pVal(iUnit) < 0.05;
end

% =========================================================================
%  SAVE
% =========================================================================
if flgSave
    save(savefile, 'rippSpks', '-v7.3');
    if flgSave
        save(savefile, 'rippSpks', '-v7.3');
    end
end

end