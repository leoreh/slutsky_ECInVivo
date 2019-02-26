function [fr, tstamps] = calcFR(spktimes, varargin)

% for each unit counts spikes in bins. can smooth result by a moving
% average (MA) or Gaussian kernel (GK) impleneted by multiple-pass MA.
% Default is to calculate firing rate in sliding 1-min windows of 20 s
% steps (Miyawaki et al., Sci. Rep., 2019). In practice this is done by
% setting binsize to 60 s and smoothing w/ moving average of 3 points. note
% that the output is given as spike count and thus to convert to Hz the
% output must be divided by the binsize. This is so that calcFR can replace
% times2binary (i.e. output a binary matrix if the binsize is small).
% 
% INPUT
% required:
%   spktimes    a cell array of vectors. each vector (unit) contains the
%               timestamps of spikes. for example {spikes.times{1:4}}
% optional:
%   winCalc     time window for calculation {[1 Inf]}. specified in s.
%   binsize     size in s of bins {60}.
%   smet        method for smoothing firing rate: moving average (MA) or
%               Gaussian kernel (GK) impleneted by multiple-pass MA.
% 
% OUTPUT
%   fr          matrix of units (rows) x firing rate in time bins (columns)
%   tstamps     vector of time stamps corresponding to beginning of bins
%
% 24 nov 18 LH. updates:
% 05 jan 18 LH  added normMethod and normWin
% 07 jan 18 LH  added disqualify units and debugging
% 11 jan 19 LH  split to various functions
% 14 jan 19 LH  added selection methods
% 24 feb 19 LH  debugging
%               replaced spikes w/ stimes as input
% 26 feb 19 LH  separated normalize and graphics to different functions
% 
% TO DO LIST
%               save which clusters pass the firing threshold and compare
%               with SU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'binsize', 60, @isscalar);
addOptional(p, 'winCalc', [1 Inf], validate_win);
addOptional(p, 'smet', 'MA', @ischar);

parse(p, varargin{:})
binsize = p.Results.binsize;
winCalc = p.Results.winCalc;
smet = p.Results.smet;

% validate window
if winCalc(1) < 1; winCalc(1) = 1; end
if winCalc(2) == Inf
    for i = 1 : length(spktimes)
        recDur(i) = max(spktimes{i}(:, 1));
    end
    winCalc(2) = max(recDur);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% count spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstamps = winCalc(1) * binsize : binsize : winCalc(2);
nbins = length(tstamps);
nunits = length(spktimes);

% count number of spikes in bins
fr = zeros(nunits, nbins);
for i = 1 : nunits
    for j = 1 : nbins - 1
        % correct for last bin
        if j > tstamps(end) / binsize
            binsize = mod(tstamps(end), binsize) * binsize;
        end
        fr(i, j) = sum(spktimes{i} > tstamps(j) & spktimes{i} < tstamps(j + 1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth FR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch smet
    case 'MA'
        fr = movmean(fr, 3);
    case 'GK'
        gk = gausswin(10);
        for i = 1 : nunits
            fr(i, :) = conv(fr(i, :), gk, 'same');
        end
end

end

% EOF