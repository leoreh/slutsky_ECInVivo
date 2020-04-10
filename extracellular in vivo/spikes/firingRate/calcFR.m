function [fr, binedges, bincents] = calcFR(spktimes, varargin)

% for each unit counts spikes in bins. can smooth result by a moving
% average (MA) or Gaussian kernel (GK) impleneted by multiple-pass MA.
% Default is to calculate firing rate in sliding 1-min windows of 20 s
% steps (Miyawaki et al., Sci. Rep., 2019). In practice this is done by
% setting binsize to 60 s and smoothing w/ moving average of 3 points.
% output may be given as spike count and thus to convert to Hz the output
% must be divided by the binsize. This is so that calcFR can replace
% times2binary (i.e. output a binary matrix if the binsize is small).
% spktimes and binsize must be the same units (e.g. seconds or samples).
% 
% INPUT
% required:
%   spktimes    a cell array of vectors where each vector (unit) contains the
%               timestamps of spikes. for example {spikes.times{1:4}}
% optional:
%   winCalc     time window for calculation {[1 Inf]}. 
%               Second elemet should be recording duration (e.g.
%               lfp.duration). if Inf will be the last spike.
%   binsize     size of bins {60}. can any units so long corresponse to
%               spktimes
%   c2r         convert counts to rate by dividing with binsize 
%   smet        method for smoothing firing rate: moving average (MA) or
%               Gaussian kernel (GK) impleneted by multiple-pass MA.
% 
% OUTPUT
%   fr          matrix of units (rows) x firing rate in time bins (columns)
%   binedegs    used for calculation
%   bincents    center of bins
%
% 24 nov 18 LH. updates:
% 05 jan 18 LH  added normMethod and normWin
% 07 jan 18 LH  added disqualify units and debugging
% 11 jan 19 LH  split to various functions
% 14 jan 19 LH  added selection methods
% 24 feb 19 LH  debugging
%               replaced spikes w/ stimes as input
% 26 feb 19 LH  separated normalize and graphics to different functions
% 02 jan 20 LH  adapted for burst suppression
%               better handling of bins
%               option to divide with binsize to obtain rate
%               handle vectors
% 
% TO DO LIST


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'binsize', 60, @isscalar);
addOptional(p, 'winCalc', [1 Inf], validate_win);
addParameter(p, 'c2r', true, @islogical);
addOptional(p, 'smet', 'MA', @ischar);

parse(p, varargin{:})
binsize = p.Results.binsize;
winCalc = p.Results.winCalc;
c2r = p.Results.c2r;
smet = p.Results.smet;

% arrange in cell
if ~iscell(spktimes)
    spktimes = {spktimes};
end

nunits = length(spktimes);

% validate window
recDur = zeros(1, nunits);
if winCalc(1) < 1; winCalc(1) = 1; end
if winCalc(2) == Inf
    for i = 1 : length(spktimes)
        recDur(i) = max(spktimes{i}(:, 1));
    end
    winCalc(2) = max(recDur);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binedges = winCalc(1) : binsize : winCalc(2);
% correct last bin
binmod = mod(winCalc(2), binsize);
binedges(end) = binedges(end) + binmod - 1;
nbins = length(binedges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% count spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% count number of spikes in bins
fr = zeros(nunits, nbins - 1);
bincents = zeros(1, nbins - 1);
for i = 1 : nunits
    for j = 1 : nbins - 1
        fr(i, j) = sum(spktimes{i} > binedges(j) &...
            spktimes{i} < binedges(j + 1));
        bincents(j) = binedges(j) + ceil(diff(binedges(j : j + 1)) / 2);
    end
end

% divide by binsize to produce rate
if c2r
    fr(:, 1 : end - 1) = fr(:, 1 : end - 1) / binsize;
    fr(:, end) = fr(:, end) / (binedges(end) - binedges(end - 1));
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

% validate orientation
if isvector(fr)
    fr = fr(:);
end

end

% EOF