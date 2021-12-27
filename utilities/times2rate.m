function [r, binedges, bincents, binidx] = times2rate(spktimes, varargin)

% for each unit counts spikes in bins. output may be given as spike count
% and thus to convert to Hz the output must be divided by the binsize. This
% is so that times2rate can replace times2binary (i.e. output a binary
% matrix if the binsize is small). spktimes, binsize and winCalc must be
% the same units (e.g. seconds or samples).
% 
% INPUT
% required:
%   spktimes    a cell array of vectors where each vector (unit) contains the
%               timestamps of spikes. for example {spikes.times{1:4}}. can
%               also be a single vector
% optional:
%   winCalc     time window for calculation {[0 Inf]}. 
%               last elemet should be recording duration (e.g.
%               lfp.duration). if Inf will be the last spike. can also be
%               an n x 2 matrix.
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

% arrange in cell
if ~iscell(spktimes)
    spktimes = {spktimes};
end

nunits = length(spktimes);

% validate window
if isempty(winCalc)
    winCalc = [1 max(vertcat(spktimes{:}))];
end
if winCalc(end) == Inf
    winCalc(end) = max(vertcat(spktimes{:}));
end
nwin = size(winCalc, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% currently assumes bins are consecutive !!!
k = 1;
for i = 1 : nwin
    if winCalc(i, 2) - winCalc(i, 1) <= binsize
        binedges{i} = winCalc(i, :);
    else
        binedges{i} = winCalc(i, 1) : binsize : winCalc(i, 2);
        binedges{i}(end) = winCalc(i, 2);       % correct last bin
    end
 
    % count spikes
    nbins = length(binedges{i}) - 1;
    for ii = 1 : nunits
        if c2r
            r(ii, k : k + nbins - 1) = histcounts(spktimes{ii}, binedges{i},...
                'Normalization', 'countdensity');
        else
            r(ii, k : k + nbins - 1) = histcounts(spktimes{ii}, binedges{i},...
                'Normalization', 'count');
        end
    end
    
    % find bin centers
    for ii = 1 : nbins
        bincents(k) = binedges{i}(ii) + ceil(diff(binedges{i}(ii : ii + 1)) / 2);
        binidx(k) = i;
        k = k + 1;
    end
end

% validate orientation
if isvector(r)
    r = r(:);
end

end

% EOF