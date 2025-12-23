function rcv = mea_frRecovery(t, frMat, varargin)
% MEA_FRRECOVERY Calculates firing rate recovery metrics.
%
%   rcv = MEA_FRRECOVERY(T, FRMAT) calculates recovery metrics for firing
%   rate data aligned to a perturbation at t=0.
%
%   INPUTS:
%       t           - (vector) Time vector in seconds. Perturbation at t=0.
%       frMat       - (matrix) Firing rate matrix (units x time).
%
%   OPTIONAL (Key-Value Pairs):
%       uGood       - (logical) Vector of good units. {all true}
%       flgSave     - (logical) Whether to save the output struct. {false}
%       basepath    - (char) Base path for saving files. {pwd}
%       binSize     - (double) Bin size in seconds. {60}
%
%   OUTPUTS:
%       rcv         - (struct) Structure containing recovery metrics:
%                       .frBsl      - Baseline firing rate [Hz].
%                       .frTrough   - Trough firing rate [Hz].
%                       .frSs       - Steady-state firing rate [Hz].
%                       .pertDepth  - Log2 ratio of Baseline to Trough.
%                       .rcvErr     - Recovery error (abs log2 fold change).
%                       .rcvGain    - Recovery gain (log2 fold change).
%                       .rcvTime    - Time to detection of recovery [s].
%                       .rcvSlope   - Initial recovery slope [Hz/min].
%                       .normSlope  - Normalized slope (% range / min).
%                       .spkDfct    - Spike deficit metric.
%                       .uRcv       - Logical indicating successful recovery.
%                       .uPert      - Logical indicating significant perturbation.
%
%   See also: MEA_FRPREP, MEA_FRR

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 't', @isnumeric);
addRequired(p, 'frMat', @isnumeric);
addParameter(p, 'uGood', [], @(x) islogical(x) || isnumeric(x));
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'basepath', pwd, @ischar);

addParameter(p, 'binSize', 60, @isnumeric);

parse(p, t, frMat, varargin{:});
uGood = p.Results.uGood;
flgSave = p.Results.flgSave;
basepath = p.Results.basepath;
binSize = p.Results.binSize;

% Infer defaults
nUnits = size(frMat, 1);
if isempty(uGood)
    uGood = true(nUnits, 1);
end
uGood = logical(uGood);

% Ensure time is a row vector
if iscolumn(t)
    t = t';
end


%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Find perturbation index (t=0)
idxPert = find(t >= -1e-9, 1, 'first');

% Constants for numerical stability
% Hardcode c based on 1 spike in an hour
c(1) = 1 / 3600;
c(2) = c(1);

% Define Windows
% Window margins in bins
winMarg = round(5 * 60 / binSize);
winDur = round(60 * 60 / binSize);

% 1. Baseline: 1 hour window ending 5 min before perturbation
bslEnd = idxPert - winMarg;
bslStart = bslEnd - winDur;
bslWin = max(1, bslStart) : max(1, bslEnd);

% 2. Steady State: 1 hour window ending 5 min before end
ssEnd = length(t) - winMarg;
ssStart = ssEnd - winDur;
ssWin = max(idxPert, ssStart) : min(length(t), ssEnd);

% 3. Trough: Minimum in 20 min window after perturbation
troughWinLen = round(20 * 60 / binSize);
troughWin = idxPert : min(length(t), idxPert + troughWinLen);


%% ========================================================================
%  CALCULATE METRICS
%  ========================================================================

% Initialize output vectors
frBsl = nan(nUnits, 1);
frSs = nan(nUnits, 1);
frTrough = nan(nUnits, 1);
idxTrough = nan(nUnits, 1);
bslTime = nan(nUnits, 1);
rcvTime = nan(nUnits, 1);
rcvSlope = nan(nUnits, 1);
normSlope = nan(nUnits, 1);
spkDfct = nan(nUnits, 1);

% Numerical Rates
% -------------------------------------------------------------------------
% Baseline
if ~isempty(bslWin)
    frBsl = mean(frMat(:, bslWin), 2, 'omitnan');
end

% Steady State
if ~isempty(ssWin)
    frSs = mean(frMat(:, ssWin), 2, 'omitnan');
end

% Trough
if ~isempty(troughWin)
    [minVal, idxMin] = min(frMat(:, troughWin), [], 2);
    frTrough = minVal;

    % Refine trough by averaging around the minimum (model-free approach)
    % idxTrough is the index in the full time vector
    idxTrough = idxPert + idxMin - 1;

    for iUnit = 1:nUnits
        % Window around the minimum
        iTrough = idxTrough(iUnit);
        currWin = max(1, iTrough) : min(length(t), iTrough + winMarg);
        frTrough(iUnit) = mean(frMat(iUnit, currWin), 2, 'omitnan');
    end
end


% Recovery Metrics (Ratios)
% -------------------------------------------------------------------------

% Cap FR values used in ratios
frBslC = max([frBsl, repmat(c(2), nUnits, 1)], [], 2);
frTroughC = max([frTrough, repmat(c(1), nUnits, 1)], [], 2);
frSsC = max([frSs, repmat(c(2), nUnits, 1)], [], 2);

% Perturbation Depth
pertDepth = log2(frBslC ./ frTroughC);

% Recovery Fidelity
rcvErr = abs(log2(frSsC ./ frBslC));
rcvDiff = abs(frSsC - frBslC);

% Recovery Gain and Work
rcvGain = log2(frSsC ./ frTroughC);
rcvGain(rcvGain < 0) = 0; % Cap at 0
rcvWork = rcvGain ./ pertDepth;


% Success Determination
% -------------------------------------------------------------------------
thrPert = 1;
thrRcv = 0.67;
uPert = pertDepth >= thrPert;
uRcv = uPert & (frSs >= frTrough + thrRcv * (frBsl - frTrough));
uRcv(~uGood) = false; % Only good units


% Kinetics
% -------------------------------------------------------------------------
tPost = t(idxPert:end); % Post-perturbation time vector

for iUnit = 1:nUnits
    if ~uGood(iUnit)
        continue;
    end

    % SPIKE DEFICIT
    % AUC Expected vs Observed
    if ~isnan(frBsl(iUnit))
        aucExp = trapz(tPost, ones(size(tPost)) * frBsl(iUnit));
        aucObs = trapz(tPost, frMat(iUnit, idxPert:end));
        spkDfct(iUnit) = log2(aucExp / (aucObs + 1));
    end

    % RECOVERY TIME
    % Start looking from the trough
    if isnan(idxTrough(iUnit))
        continue;
    end

    currIdxTrough = idxTrough(iUnit);
    rcvCurve = frMat(iUnit, currIdxTrough:end);

    % Baseline target
    bslTrgt = frTrough(iUnit) + thrRcv * (frBsl(iUnit) - frTrough(iUnit));
    bslPnt = find(rcvCurve >= bslTrgt, 1, 'first');
    if ~isempty(bslPnt)
        bslTime(iUnit) = (currIdxTrough + bslPnt - 1 - idxPert) * binSize;
    else
        bslTime(iUnit) = t(end);
    end

    % Steady-state target
    rcvRange = frSs(iUnit) - frTrough(iUnit);
    rcvTrgt = frTrough(iUnit) + thrRcv * rcvRange;
    rcvPnt = find(rcvCurve >= rcvTrgt, 1, 'first');
    if ~isempty(rcvPnt)
        rcvTime(iUnit) = (currIdxTrough + rcvPnt - 1 - idxPert) * binSize;
    else
        rcvTime(iUnit) = t(end);
    end

    % SLOPE
    if ~isempty(rcvCurve) && length(rcvCurve) > 1
        rcvSlope(iUnit) = (max(diff(rcvCurve)) / binSize) * 60; % Hz/min
        if rcvRange ~= 0
            normSlope(iUnit) = (rcvSlope(iUnit) / rcvRange) * 100;
        end
    end
end


%% ========================================================================
%  OUTPUT STRUCT
%  ========================================================================

rcv.frBsl = frBsl;
rcv.frTrough = frTrough;
rcv.frSs = frSs;
rcv.pertDepth = pertDepth;
rcv.rcvErr = rcvErr;
rcv.rcvGain = rcvGain;
rcv.rcvWork = rcvWork;
rcv.rcvDiff = rcvDiff;
rcv.bslTime = bslTime;
rcv.rcvTime = rcvTime;
rcv.rcvSlope = rcvSlope;
rcv.normSlope = normSlope;
rcv.spkDfct = spkDfct;
rcv.uRcv = uRcv;
rcv.uPert = uPert;
rcv.idxTrough = idxTrough;
rcv.info.binSize = binSize;
rcv.info.tPert = 0;

if flgSave
    [~, basename] = fileparts(basepath);
    save(fullfile(basepath, [basename, '.frRcv.mat']), 'rcv', '-v7.3');
end

end
