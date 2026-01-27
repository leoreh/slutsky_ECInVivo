function [r, binEdges, binCents] = times2rate(times, varargin)
% TIMES2RATE Calculates event rates/counts in time bins.
%
%   [r, binEdges, binCents] = TIMES2RATE(TIMES, ...) calculates the rate
%   or count of events (e.g., spikes) falling into time bins. Bins can be
%   defined by a fixed 'binSize' within larger calculation windows ('winCalc').
%
%   INPUTS:
%       times       - (cell/vec) Timestamps per unit (cell array or single vector).
%
%   OPTIONAL (Key-Value Pairs):
%       'winCalc'   - (num) [M x 2] matrix defining M time windows [start, end].
%                     Default: [0, max_time_across_all_units].
%       'binsize'   - (num) Size of bins. Default: 60.
%                     If Inf, each 'winCalc' interval is a single bin.
%       'c2r'       - (log) True to convert counts to rate. Default: true.
%
%   OUTPUTS:
%       r           - (num) Matrix [nUnits x nBins] of rates or counts.
%       binEdges    - (num) Matrix [nBins x 2] of [start, end] times for each bin.
%       binCents    - (num) Vector [1 x nBins] of center times for each bin.
%
%   HISTORY:
%       02 Aug 20   Converted calcFR to times2rate. Implemented histcounts.
%       19 Nov 20   winCalc can be matrix (e.g. for states).
%       11 May 25   Implemented efficient 'histcounts' strategy.
%       23 Dec 23   Refactored for style and n2chunks integration.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'times');
addParameter(p, 'binsize', 60, @isnumeric);
addParameter(p, 'winCalc', [], @isnumeric);
addParameter(p, 'c2r', true, @islogical);

parse(p, times, varargin{:});
binSize = p.Results.binsize;
winCalc = p.Results.winCalc;
c2r     = p.Results.c2r;


%% ========================================================================
%  PREPARATIONS
%  ========================================================================

% Ensure 'times' is a cell array
if ~iscell(times)
    times = {times};
end
nUnits = length(times);

% Get max time across all units if winCalc is empty
if isempty(winCalc)
    winCalc = [0, Inf];
end
maxTime = max(cellfun(@(x) max([0; x(:)]), times));
if isinf(winCalc(end))
    winCalc(end) = maxTime;
end
if winCalc(1) < 0
    winCalc(1) = 0;
end
nWin = size(winCalc, 1);


%% ========================================================================
%  DEFINE BINS
%  ========================================================================
% Accumulate effective bin edges from all windows

if isinf(binSize)
    binEdges = winCalc;
else
    binEdges = zeros(0, 2);

    for iWin = 1:nWin
        winStart = winCalc(iWin, 1);
        winEnd   = winCalc(iWin, 2);
        winDur   = winEnd - winStart;

        if winDur <= binSize
            edges = [0, winDur];
        else
            % Use n2chunks to split duration
            % 'exclude' ensures we only keep full bins (dropping tail)
            chunks = n2chunks('n', winDur, 'chunksize', binSize, 'lastChunk', 'exclude');

            % Convert chunks (1-based integer intervals) to 0-based time offsets
            edges = [chunks(:,1)-1, chunks(:,2)];
        end

        % Adjust for window start and append
        binEdges = [binEdges; edges + winStart];
    end
end

nBins   = size(binEdges, 1);
binDurs = binEdges(:, 2) - binEdges(:, 1);
binCents = mean(binEdges, 2);


%% ========================================================================
%  CALCULATE RATES
%  ========================================================================

% Initialize Output (Default to NaN to match previous style, though typically should be 0)
r = nan(nUnits, nBins);

% Efficiency Strategy: "Odd Bin Trick"
% Create a single, interleaved edges vector.
histEdges = reshape(binEdges', [], 1);

% Check for Overlap / Non-Monotonicity
if issorted(histEdges)
    % --- Pathway A: Sorted, Non-Overlapping (Fastest) ---
    for iUnit = 1:nUnits
        unitTimes = times{iUnit};

        if isempty(unitTimes), continue; end

        cnts = histcounts(unitTimes, histEdges);
        cnts = cnts(1:2:end);

        if c2r
            r(iUnit, :) = cnts(:)' ./ binDurs(:)';
        else
            r(iUnit, :) = cnts;
        end
    end
else
    % --- Pathway B: Overlapping or Disordered (Robust) ---
    % Use cumulative counts on unique sorted edges
    [uniqEdges, ~, splitIdx] = unique(histEdges);

    % Map original start/end times to indices in uniqEdges
    binIdx = reshape(splitIdx, 2, [])';

    for iUnit = 1:nUnits
        unitTimes = times{iUnit};

        if isempty(unitTimes), continue; end

        % Count in unique atomic segments
        baseCnts = histcounts(unitTimes, uniqEdges);

        % Cumulative sum (prepend 0)
        cumCnts = [0, cumsum(baseCnts)];

        % Calculate bin counts by difference: Cum(End) - Cum(Start)
        cnts = cumCnts(binIdx(:,2)) - cumCnts(binIdx(:,1));

        if c2r
            r(iUnit, :) = cnts ./ binDurs(:)';
        else
            r(iUnit, :) = cnts;
        end
    end
end

end     % EOF
