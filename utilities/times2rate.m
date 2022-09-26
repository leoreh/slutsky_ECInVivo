function [r, binedges, bincents, binidx] = times2rate(times, varargin)

% gets a vector of times and converts to rate (counts the number of
% occurances) in a specified binsize. essentially a wrapper for histcounts.
% e.g., used to converts spikes.times to firing rate. can limit the bins to
% the boundries specified by winCalc.
% 
% INPUT
%   times       cell array of vectors where each vector (unit) contains the
%               timestamps of spikes. for example {spikes.times{1:4}}. can
%               also be a single vector
%   winCalc     time window for calculation {[0 Inf]}. 
%               last elemet should be recording duration (e.g.
%               lfp.duration). if Inf will be the last spike. can also be
%               an n x 2 matrix.
%   binsize     size of bins {60}. can be of any units so long corresponse
%               to times. for example, times from Wh (per tetrode) is
%               given in samples whereas times from ce (per unit) is
%               given in seconds
%   c2r         logical. convert counts to rate by dividing with binsize 
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
% 19 nov 20 LH  winClac can be matrix (e.g. for states)
% 
% TO DO LIST
%   # documentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'binsize', 60, @isscalar);
addOptional(p, 'winCalc', [0 Inf], @isnumeric);
addOptional(p, 'c2r', true, @islogical);

parse(p, varargin{:})
binsize     = p.Results.binsize;
winCalc     = p.Results.winCalc;
c2r         = p.Results.c2r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare spike times
if ~iscell(times)
    times = {times};
end
nunits = length(times);

% validate window
if isempty(winCalc)
    winCalc = [1 max(vertcat(times{:}))];
end
if winCalc(end) == Inf
    winCalc(end) = max(vertcat(times{:}));
end
nwin = size(winCalc, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% currently assumes bins are consecutive !!!
cnt = 1;
for iwin = 1 : nwin
    if winCalc(iwin, 2) - winCalc(iwin, 1) <= binsize
        binedges{iwin} = winCalc(iwin, :);
    else
        binedges{iwin} = winCalc(iwin, 1) : binsize : winCalc(iwin, 2);
        binedges{iwin}(end) = winCalc(iwin, 2);       % correct last bin
    end
 
    % count spikes
    nbins = length(binedges{iwin}) - 1;
    for iunit = 1 : nunits
        if c2r
            r(iunit, cnt : cnt + nbins - 1) = histcounts(times{iunit}, binedges{iwin},...
                'Normalization', 'countdensity');
        else
            r(iunit, cnt : cnt + nbins - 1) = histcounts(times{iunit}, binedges{iwin},...
                'Normalization', 'count');
        end
    end
    
    % find bin centers
    for ibin = 1 : nbins
        bincents(cnt) = binedges{iwin}(ibin) + ceil(diff(binedges{iwin}(ibin : ibin + 1)) / 2);
        binidx(cnt) = iwin;
        cnt = cnt + 1;
    end
end

% validate orientation
if isvector(r)
    r = r(:);
end

end

% EOF