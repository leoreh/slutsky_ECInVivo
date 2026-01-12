function [idxTrough, tAxis, debugData] = mcu_detectTrough(frMat, idxPert, varargin)
% MCU_DETECTTROUGH Detects global onset of recovery (end of silence).
%
%   [idxTrough, tAxis, debugData] = MCU_DETECTTROUGH(frMat, idxPert) calculates
%   the population median and finds the "wake-up" point after the perturbation
%   where activity consistently rises above a threshold.
%
%   INPUTS:
%       frMat       - (matrix) Firing rate matrix [Units x Time].
%       idxPert     - (scalar) Index of perturbation onset.
%
%   OPTIONAL (Key-Value Pairs):
%       binSize     - (double) Bin size in seconds. {60 s}
%       margMin   - (double) Margin iskip after pert. {5 bins}
%       flgPlot     - (logical) Whether to plot results. {false}
%
%   OUTPUTS:
%       idxTrough   - (scalar) Global index of recovery onset.
%       tAxis       - (vector) Time axis vector (hours).
%       debugData   - (struct) Intermediate data for debugging.
%
%   See also: MCU_DETECTPERT, MEA_FRFIT

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'frMat', @isnumeric);
addRequired(p, 'idxPert', @isnumeric);
addParameter(p, 'binSize', 60, @isnumeric);
addParameter(p, 'margMin', 5, @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);

parse(p, frMat, idxPert, varargin{:});
binSize = p.Results.binSize;
margMin = p.Results.margMin;
flgPlot = p.Results.flgPlot;

%% ========================================================================
%  DETECTION LOGIC
%  ========================================================================

% Population Trend (Median) - more robust to outliers
mfr = median(frMat, 1, 'omitnan');

% Denoise the trace
mfr = fr_denoise(mfr, 'flgPlot', false, 'frameLen', 10);

% Baseline (until idxPert - margin)
margBins = round(margMin * 60 / binSize);
frBsl = median(mfr(1 : idxPert - margBins), 'omitnan');

% Manual inspection shows FR finishes drop <10 bins (1hr abs) after pert
% onset. Note, this type of visual inspection should be done on the same
% trace as that which detected idxPert (eg, smoothed FR in 60s bins).
srchStart = idxPert + 10 + margBins;
srchEnd = srchStart + 60;
postWin = mfr(srchStart:srchEnd);

% Define Threshold
% Find minimal activity level (floor) in the post-perturbation window
% Use Xth percentile to avoid outlier zeros
floorVal = prctile(postWin, 1);

% Threshold: 10% of recovery range, bounded by absolute min
thr = max(floorVal + 0.1 * (frBsl - floorVal), floorVal + 0.1);

% Find "Wake-up" Point
% Look for first point where:
%   a) Value > Threshold
%   b) Stays > Threshold (or > 0.8*Threshold) for a window (e.g., 5 bins)
nBins = length(postWin);
winSize = 5;
isAbove = postWin > thr;

idxRcvRel = [];

for iBin = 1 : (nBins - winSize)
    if isAbove(iBin)
        % Check if sustained
        if mean(postWin(iBin : iBin+winSize)) > thr * 0.8
            idxRcvRel = iBin;
            break;
        end
    end
end

% Refine Trough (Backtrack)
if isempty(idxRcvRel)
    % No recovery found: use minimum in the window
    [~, idxMin] = min(postWin);
    idxRcvRel = idxMin;
else
    % Walk back from the "wake-up" point to find where the rise started.
    % We travel backwards until the signal drops close to the floor level.
    % Stop when fr <= floorVal + small_tolerance OR derivative is small/negative

    backtrackTol = floorVal; 

    % Default to start of search window if we can't find a better point
    newIdx = 1;

    for iBin = idxRcvRel : -1 : 2
        currVal = postWin(iBin);
        prevVal = postWin(iBin-1);

        % If we hit the floor level (or close to it)
        if currVal <= backtrackTol
            newIdx = iBin;
            break;
        end
    end
    idxRcvRel = newIdx;
end

% Convert back to global index
idxTrough = srchStart + idxRcvRel - 1;


%% ========================================================================
%  OUTPUT & PLOT
%  ========================================================================

% Time axis (hours)
nBinsTotal = length(mfr);
dt_hr = binSize / 3600;
tAxis = ((1:nBinsTotal) - idxPert) * dt_hr;

debugData.frPop = mfr;
debugData.frBsl = frBsl;
debugData.thr = thr;
debugData.floorVal = floorVal;

if flgPlot
    plot_axSize('szOnly', false, 'axShape', 'wide', 'axHeight', 300);
    hold on;

    plot(tAxis, mfr, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Median FR');
    yline(frBsl, 'b--', 'DisplayName', 'Baseline');
    yline(thr, 'g:', 'DisplayName', 'Threshold');

    tTrough = tAxis(idxTrough);
    xline(tTrough, 'r-', 'LineWidth', 2, 'DisplayName', 'Detected Trough');

    tPert = tAxis(idxPert);
    xline(tPert, '--', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Perturbation');

    % Mark Margin
    tMargin = tAxis(srchStart);
    xline(tMargin, ':', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Search Start');

    title('Population Recovery Detection');
    xlabel('Time (hr)');
    ylabel('Median FR (Hz)');
    legend('show');
end

end     % EOF
