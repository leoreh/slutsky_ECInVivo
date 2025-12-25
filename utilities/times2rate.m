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
%       'binSize'   - (num) Size of bins. Default: 60.
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
maxTime = max(cellfun(@(x) max([0; x(:)]), times));
if isempty(winCalc)
    winCalc = [0, maxTime];
end
if isinf(winCalc(end))
    winCalc(end) = maxTime;
end
nWin = size(winCalc, 1);


%% ========================================================================
%  DEFINE BINS
%  ========================================================================

% Accumulate effective bin edges from all windows
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

nBins   = size(binEdges, 1);
binDurs = binEdges(:, 2) - binEdges(:, 1);
binCents = mean(binEdges, 2);


%% ========================================================================
%  CALCULATE RATES
%  ========================================================================

% Initialize Output (Default to NaN to match previous style, though typically should be 0)
r = nan(nUnits, nBins);

% Efficiency Strategy: "Odd Bin Trick"
% Create a single, interleaved edges vector for histcounts.
histEdges = reshape(binEdges', [], 1);

for iUnit = 1:nUnits
    unitTimes = times{iUnit};

    % Skip if no spikes for this unit
    if isempty(unitTimes)
        continue;
    end

    % Get counts
    cnts = histcounts(unitTimes, histEdges);

    % Keep only our defined bins (odd-indexed intervals in the interleaved vector)
    cnts = cnts(1:2:end);

    if c2r
        r(iUnit, :) = cnts(:)' ./ binDurs(:)';
    else
        r(iUnit, :) = cnts;
    end
end

end     % EOF
