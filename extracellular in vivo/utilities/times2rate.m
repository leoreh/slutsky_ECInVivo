function [r, binedges, bincents] = times2rate(itimes, varargin)

% for each unit counts spikes in bins. output may be given as spike count
% and thus to convert to Hz the output must be divided by the binsize. This
% is so that times2rate can replace times2binary (i.e. output a binary
% matrix if the binsize is small). spktimes and binsize must be the same
% units (e.g. seconds or samples).
% 
% INPUT
% required:
%   itimes      a cell array of vectors where each vector (unit) contains the
%               timestamps of spikes. for example {spikes.times{1:4}}
% optional:
%   winCalc     time window for calculation {[1 Inf]}. 
%               Second elemet should be recording duration (e.g.
%               lfp.duration). if Inf will be the last spike.
%   binsize     size of bins {60}. can be of any units so long corresponse
%               to spktimes. for example, spktimes from Wh (per tetrode) is
%               given in samples whereas spktimes from ce (per units) is
%               given in seconds
%   c2r         convert counts to rate by dividing with binsize 
% 
% OUTPUT
%   r           matrix of units (rows) x rate in time bins (columns)
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
% 02 aug 20 LH  converted calcFR to times2rate. implemented histcounts
% 
% TO DO LIST
%               update description


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'binsize', 60, @isscalar);
addOptional(p, 'winCalc', [0 Inf], validate_win);
addOptional(p, 'c2r', true, @islogical);

parse(p, varargin{:})
binsize     = p.Results.binsize;
winCalc     = p.Results.winCalc;
c2r         = p.Results.c2r;

% arrange in cell
if ~iscell(itimes)
    itimes = {itimes};
end

nunits = length(itimes);

% validate window
if winCalc(2) == Inf
    winCalc(2) = max(vertcat(itimes{:}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binedges = winCalc(1) : binsize : winCalc(2);
% correct last bin
binmod = mod(winCalc(2), binsize);
binedges(end) = binedges(end) + binmod - 1;
nbins = length(binedges);
bincents = zeros(1, nbins - 1);
for i = 1 : nbins - 1
    bincents(i) = binedges(i) + ceil(diff(binedges(i : i + 1)) / 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% count spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% count number of spikes in bins
r = zeros(nunits, nbins - 1);
for i = 1 : nunits
    if c2r
        r(i, :) = histcounts(itimes{i}, binedges,...
            'Normalization', 'countdensity');
    else
        r(i, :) = histcounts(itimes{i}, binedges,...
            'Normalization', 'count');
    end
end

% validate orientation
if isvector(r)
    r = r(:);
end

end

% EOF