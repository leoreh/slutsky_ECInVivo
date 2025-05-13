function [r, binEdges, binCents] = times2rate(times, varargin)

% Calculates the rate or count of events (e.g., spikes) falling into time bins.
% Bins can be defined by a fixed 'binsize' within larger calculation windows
% ('winCalc'), or 'winCalc' intervals can themselves be treated as individual bins.
%
% INPUT:
%   times       Cell array of vectors (timestamps per unit), or single vector.
%   winCalc     (Optional) [M x 2] matrix defining M time windows [start, end].
%               Default: [0, max_time_across_all_units].
%   binsize     (Optional) Scalar, size of bins. Default: 60.
%               If Inf, each 'winCalc' interval is a single bin.
%   c2r         (Optional) Logical, true to convert counts to rate. Default: true.
%
% OUTPUT:
%   r           Matrix [nUnits x nBins] of rates or counts.
%   binEdges    [nBins x 2] matrix of [start, end] times for each bin.
%   binCents    Vector [1 x nBins] of center times for each bin.
%
% HISTORY:
% 02 Aug 20     Converted calcFR to times2rate. Implemented histcounts.
% 19 Nov 20     winCalc can be matrix (e.g. for states).
% 11 May 25     Implemented efficient 'histcounts' strategy.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Argument Parsing (Minimal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addRequired(p, 'times');
addParameter(p, 'binsize', 60);
addParameter(p, 'winCalc', []); % Default handled below
addParameter(p, 'c2r', true);
parse(p, times, varargin{:});

binSize     = p.Results.binsize;
winCalc     = p.Results.winCalc;
c2r         = p.Results.c2r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparations (Minimal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure 'times' is a cell array
if ~iscell(times)
    times = {times};
end
nUnits = length(times);

% Get max time across all units if winCalc is empty
maxTimes = max(cellfun(@(x) max([0; x(:)]), times));
if isempty(winCalc)  
    winCalc = [0, maxTimes];
end
if isinf(winCalc(end))
    winCalc(end) = maxTimes;
end
nWin = size(winCalc, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define All Effective Bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section determines all the actual [start, end] time bins that will be
% used for calculation, based on 'winCalc' and 'binsize'.

bins = cell(1, nWin);
for iWin = 1:nWin
    winDur = winCalc(iWin, 2) - winCalc(iWin, 1);
    
    if isinf(binSize) || winDur <= binSize
        bins{iWin} = winCalc(iWin, :);
    else
        bins{iWin} = winCalc(iWin, 1) : binSize : winCalc(iWin, 2);
        bins{iWin}(end) = winCalc(iWin, 2);       % correct last bin
    end
end

% Convert cell array of bin edges to matrix [nBins x 2]
binEdges = zeros(0, 2);
for iWin = 1:nWin
    winDur = winCalc(iWin, 2) - winCalc(iWin, 1);
    if isinf(binSize) || winDur <= binSize
        binEdges = [binEdges; winCalc(iWin, :)];
    else
        binStarts = bins{iWin}(1:end-1);
        binEnds = bins{iWin}(2:end);
        binEdges = [binEdges; [binStarts', binEnds']];
    end
end

nBins = size(binEdges, 1);
binDurs = binEdges(:,2) - binEdges(:,1); % Durations of each bin [nBins x 1]

% initialize
r = nan(nUnits, nBins);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Rates or Counts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses the "odd bin trick" with histcounts for efficiency.
% For each unit, counts spikes falling into all bins with a single histcounts call.

% Create a single, interleaved edges vector for histcounts:
% e.g., [bin1_start, bin1_end, bin2_start, bin2_end, ...]
histEdges = reshape(binEdges', [], 1); % [2*nBins x 1]

for iUnit = 1:nUnits
    unitTimes = times{iUnit};

    % Skip if no spikes for this unit 
    if isempty(unitTimes)
        continue; 
    end

    % Get counts; countsRaw(1) is between histEdges(1) and histEdges(2),
    % etc. Counts for our actual defined bins are in the odd-indexed
    % elements. [1 x nBins]
    cnts = histcounts(unitTimes, histEdges);   
    cnts = cnts(1:2:end); 

    if c2r % Convert to rate       
        r(iUnit, :) = cnts(:)' ./ binDurs(:)';
    else % Keep raw counts
        r(iUnit, :) = cnts;
    end
end

% Calculate bin centers from the actual bin edges used
binCents = mean(binEdges, 2);  % Row vector [1 x nBins] of bin centers

end 
