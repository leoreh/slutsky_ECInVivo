function fr = FR(spktimes, varargin)

% wrapper for firing rate functions
% 
% INPUT
%   spktimes    a cell array of vectors. each vector (unit) contains the
%               timestamps of spikes. for example {spikes.times{1:4}}
%   basepath    recording session path {pwd}
%   graphics    plot figure {1}.
%   saveFig     save figure {1}.
%   saveVar     save variable {1}.
%   winCalc     time window for calculation {[1 Inf]}. specified in s.
%   binsize     size bins {60}. specified in s.
%   metBL       calculate baseline as 'max' or {'avg'}.
%   winBL       window to calculate baseline FR {[1 Inf]}.
%               specified in s.
%   select      cell array with strings expressing method to select units.
%               'thr' - units with fr > 0.05 during baseline
%               'stable' - units with std of fr < avg of fr during
%               baseline. default = none.
%   smet        method for smoothing firing rate: moving average (MA) or
%               Gaussian kernel (GK) impleneted by multiple-pass MA. {[]}.
% 
% OUTPUT
% fr            struct with fields strd, norm, bins, binsize,
%               normMethod, normWin
%
% 26 feb 19 LH. 
% 21 nov 19 LH  added active periods and mFR accordingly
% 09 may 20 LH  fixed c2r issues with calcFR

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
addOptional(p, 'metBL', 'avg', @ischar);
addOptional(p, 'select', {'thr'}, @iscell);
addOptional(p, 'smet', 'none', @ischar);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
binsize = p.Results.binsize;
winCalc = p.Results.winCalc;
winBL = p.Results.winBL;
metBL = p.Results.metBL;
select = p.Results.select;
smet = p.Results.smet;
graphics = p.Results.graphics;
saveFig = p.Results.saveFig;
saveVar = p.Results.saveVar;

% validate windows
if winCalc(1) < 1; winCalc(1) = 1; end
if winCalc(2) == Inf
    for i = 1 : length(spktimes)
        recDur(i) = max(spktimes{i}(:, 1));
    end
    winCalc(2) = max(recDur);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate spike counts
[fr.strd, ~, fr.tstamps] = calcFR(spktimes, 'binsize', binsize,...
    'winCalc', winCalc, 'c2r', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth FR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch smet
    case 'MA'
        r = movmean(r, 3);
    case 'GK'
        gk = gausswin(10);
        for i = 1 : nunits
            r(i, :) = conv(r(i, :), gk, 'same');
        end
end

% active periods and mFR
% threshold for active periods is determined according to the bimodel
% distribution of fr values. DOES NOT WORK FOR FR
% thr = sepBimodel('x', fr.strd(:), 'lognorm', true, 'graphics', true);
thr = 0.5;
act = fr.strd > thr;
act = logical(round(movmedian(act, 3)));
for i = 1 : size(fr.strd, 1)
    fr.mfr(i) = mean(fr.strd(i, act(i, :)));
end

% normalize firing rate
if isempty(winBL)
    winBL = [1 size(fr.strd, 1)];
else
    winBL = winBL / binsize;
    if winBL(1) < 1; winBL(1) = 1; end
    if winBL(2) == Inf; winBL(2) = size(fr.strd, 2); end
end
[fr.norm, fr.ithr, fr.istable] = normFR(fr.strd, 'metBL', 'avg', 'win', winBL, 'select', select);

% add params
fr.winBL = winBL;
fr.binsize = binsize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    plotFRtime('fr', fr, 'spktimes', spktimes, 'units', true,...
        'avg', true, 'raster', true, 'saveFig', saveFig)  
    
%     bl = blFR(fr.norm, 'method', metBL, 'win', winBL);
%     plotFRdistribution(bl, 'saveFig', saveFig) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar  
    [~, filename] = fileparts(basepath);
    save([basepath, filesep, filename, '.fr.mat'], 'fr')
end

end

% EOF