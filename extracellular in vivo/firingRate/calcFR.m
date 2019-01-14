function fr = calcFR(spikes, varargin)

% for each unit calculates firing frequency in Hz, defined as spike count
% per binsize divided by binsize {default 1 min}. Also calculates average
% and std of normailzed firing frequency. Normalization is according to
% average or maximum FR within a specified window.
% 
% INPUT
% required:
%   spikes      struct (see getSpikes)
% optional:
%   graphics    plot figure {1}.
%   winCalc     time window for calculation {[1 Inf]}. specified in s.
%   binsize     size in s of bins {60}.
%   saveFig     save figure {1}.
%   basepath    recording session path {pwd}
%   metBL       calculate baseline as 'max' or {'avg'}.
%   winBL       window to calculate baseline FR {[1 Inf]}.
%               specified in s.
% 
% EXAMPLES      calcFR(spikes, 'metBL', 'avg', 'winBL', [90 Inf]);
%               will normalize FR according to the average FR between 90
%               min and the end of the recording.
%
% OUTPUT
% fr            struct with fields strd, norm, avg, std, bins, binsize,
%               normMethod, normWin
%
% 24 nov 18 LH. updates:
% 05 jan 18 LH  added normMethod and normWin
% 07 jan 18 LH  added disqualify units and debugging
% 11 jan 19 LH  split to various functions
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
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'winCalc', [1 Inf], validate_win);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'basepath', pwd);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'metBL', 'avg', @ischar);
addOptional(p, 'winBL', [1 Inf], validate_win);
addOptional(p, 'disqualify', false, @islogical);

parse(p, varargin{:})
binsize = p.Results.binsize;
winCalc = p.Results.winCalc / binsize;
winBL = p.Results.winBL / binsize;
metBL = p.Results.metBL;
disqualify = p.Results.disqualify;
basepath = p.Results.basepath;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;

% validate windows
if winCalc(1) == 0; winCalc(1) = 1; end
if winCalc(2) == Inf
    if iscell(spikes.spindices)
        for i = 1 : length(spikes.spindices)
            recDur(i) = max(spikes.spindices{i}(:, 1));
        end
    else
        recDur = max(spikes.spindices(:, 1));
    end
    winCalc(2) = max(recDur) / binsize;
end
if winBL(1) == 0; winBL(1) = 1; end
if winBL(2) == Inf; winBL(2) = recDur; end
winCalc = winCalc(1) : winCalc(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nunits = length(spikes.UID);
nbins = ceil(winCalc(end) - winCalc(1) + 1);

% count number of spikes in bins
fr.strd = zeros(nunits, nbins);
for i = 1 : nunits
    for j = 1 : nbins
        % correct for last minute
        if j > winCalc(end)
            binsize = mod(winCalc(2), binsize) * binsize;
        end
        fr.strd(i, j) = sum(ceil(spikes.times{i} / binsize) == winCalc(j)) / binsize;
    end
end

% select only units who fired above thr during baseline window
bl = avgFR(fr.strd, 'method', metBL, 'win', winBL);
stable = bl > 0.05;

bl = mean(fr.strd(:, 1:6), 2);
idx = max(diff(fr.strd(:, 1 : 6)'))' < bl;
stable = bl > 0.05;
x = stable | idx

% normalize
fr.norm = fr.strd(stable, :) ./ bl(stable);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    plotFRtime('fr', fr.norm, 'spktime', spikes.times, 'units', true, 'avg', true, 'saveFig', saveFig)  
    plotFRdistribution(bl, 'saveFig', saveFig) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar
    fr.winCalc = [winCalc(1) winCalc(end)];
    fr.winBL = winBL;
    fr.binsize = binsize;
    
    [~, filename] = fileparts(basepath);
    save([basepath, '\', filename, '.fr.mat'], 'fr')
end

end

% EOF