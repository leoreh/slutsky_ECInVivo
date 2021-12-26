function fr = firingRate(spktimes, varargin)

% wrapper for firing rate functions. calculates firing rate based on spike
% times and smoothes result by a moving average (MA) or Gaussian kernel
% (GK) impleneted by multiple-pass MA. Default is to calculate firing rate
% in sliding 1-min windows of 20 s steps (Miyawaki et al., Sci. Rep.,
% 2019). In practice this is done by setting binsize to 60 s and smoothing
% w/ moving average of 3 points.
%
% INPUT
%   spktimes    a cell array of vectors. each vector (unit / tetrode)
%               contains the timestamps of spikes. for example
%               {spikes.times{1:4}}
%   basepath    recording session path {pwd}
%   graphics    plot figure {1}.
%   saveVar     logical / char. save variable {true}. if char than variable
%               will be named saveVar.mat
%   winCalc     time window for calculation {[1 Inf]}. specified in s.
%   binsize     size bins {60}. specified in s.
%   winBL       window to calculate baseline FR {[1 Inf]}.
%               specified in s.
%   smet        method for smoothing firing rate: moving average (MA) or
%               Gaussian kernel (GK) impleneted by multiple-pass MA. {[]}.
%
% OUTPUT
% fr            struct 
%
% TO DO LIST
%               adjust winCalc to matrix
%
% 26 feb 19 LH  updates:
% 21 nov 19 LH  added active periods and mFR accordingly
% 09 may 20 LH  fixed c2r issues with calcFR
% 03 feb 21 LH  added states
% 26 dec 21 LH  gini coefficient and fano factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'binsize', 60, @isscalar);
addOptional(p, 'winCalc', [1 Inf], validate_win);
addOptional(p, 'winBL', [], validate_win);
addOptional(p, 'smet', 'none', @ischar);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveVar', true);

parse(p, varargin{:})
basepath = p.Results.basepath;
binsize = p.Results.binsize;
winCalc = p.Results.winCalc;
winBL = p.Results.winBL;
smet = p.Results.smet;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;

smfactor = 7;    % smooth factor
nunits = length(spktimes);
[~, basename] = fileparts(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calc fr across entire session
[fr.strd, ~, fr.tstamps] = times2rate(spktimes, 'binsize', binsize,...
    'winCalc', winCalc, 'c2r', true);
fr.mfr = mean(fr.strd, 2);

% calc fr according to states. note states, binsize, and spktimes must be
% the same units (typically sec)
if exist(fullfile(basepath, [basename '.AccuSleep_states.mat']))
    load(fullfile(basepath, [basename '.AccuSleep_states.mat']), 'ss')    
    fr.states.stateNames = ss.labelNames;
    nstates = length(ss.stateEpochs);

    % fit stateEpochs to winCalc and calc firing rate in states
    for istate = 1 : nstates - 1
        epochIdx = ss.stateEpochs{istate}(:, 2) < winCalc(2) &...
            ss.stateEpochs{istate}(:, 1) > winCalc(1);
        stateEpochs = ss.stateEpochs{istate}(epochIdx, :);
         if ~isempty(ss.stateEpochs{istate})
            [fr.states.fr{istate}, fr.states.binedges, fr.states.tstamps{istate}, fr.states.binidx] =...
                times2rate(spktimes, 'binsize', binsize, 'winCalc', stateEpochs, 'c2r', true);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modulate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% smooth
switch smet
    case 'MA'
        fr.strd = movmean(fr.strd, smfactor, 2);
    case 'GK'
        gk = gausswin(smfactor);
        gk = gk / sum(gk);
        for i = 1 : nunits
            fr.strd(i, :) = conv(fr.strd(i, :), gk, 'same');
        end
end

% normalize firing rate
if isempty(winBL)
    winBL = [1 size(fr.strd, 2)];
else
    winBL = winBL / binsize;
    if winBL(1) < 1
        winBL(1) = 1;
    end
    if winBL(2) == Inf
        winBL(2) = size(fr.strd, 2);
    end
end
bl_fr = fr.strd(:, winBL(1) : winBL(2));
bl_avg = mean(bl_fr, 2);
bl_std = std(bl_fr, [], 2);
fr.norm = fr.strd ./ bl_avg;

% apply criterions
bl_thr = 0.01;   % [Hz]
fr.bl_thr = bl_avg > bl_thr;
fr.stable = bl_std < bl_avg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% more params 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fano factor: variability in fr relative to mfr
fr.fanoFactor = var(bl_fr, [], 2) / mean(bl_fr, 2);

% Gini coefficient. calculated per unit (as in CE); the gini describes the
% inequality of fr bins throughout time. calculated across the popultion
% (mizuseki, cell rep., 2008), the gini describes the degree to which high
% mfr units accounted for most the spikes recorded
cum_fr = cumsum(sort(bl_fr, 2), 2);
cum_fr_norm = cum_fr ./ max(cum_fr, [], 2);
for iunit = 1 : nunits
    fr.gini_unit(iunit) = gini(ones(1, size(bl_fr, 2)), cum_fr_norm(iunit, :));
end
fr.gini_pop = gini(ones(1, size(bl_fr, 1)),...
    cumsum(sort(fr.mfr)) / max(cumsum(fr.mfr)));

% AR(1): auto-regressive
% plot(fr.strd(iunit, 1 : end-1), fr.strd(iunit, 2 : end), '*')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% struct
fr.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
fr.info.winBL = winBL;
fr.info.winCalc = winCalc;
fr.info.binsize = binsize;
fr.info.smoothMethod = smet;
fr.info.bl_thr = bl_thr;

% save
if saveVar
    if ischar(saveVar)
        save([basepath, filesep, basename, '.' saveVar '.mat'], 'fr')
    else
        save([basepath, filesep, basename, '.fr.mat'], 'fr')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    plot_FRtime_session('basepath', basepath)
end

return

% EOF

