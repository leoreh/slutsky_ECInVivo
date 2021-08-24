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
%   saveFig     save figure {1}.
%   saveVar     logical / char. save variable {true}. if char than variable
%               will be named saveVar.mat
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
% TO DO LIST
%               adjust winCalc to matrix
%
% 26 feb 19 LH.
% 21 nov 19 LH  added active periods and mFR accordingly
% 09 may 20 LH  fixed c2r issues with calcFR
% 03 feb 21 LH  added states

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
addOptional(p, 'select', {''}, @iscell);
addOptional(p, 'smet', 'none', @ischar);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'saveVar', true);

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

% validate window
if winCalc(end) == Inf
    winCalc(end) = max(vertcat(spktimes{:}));
end

smfactor = 7;    % smooth factor
nunits = length(spktimes);
[~, basename] = fileparts(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calc fr across entire session
[fr.strd, ~, fr.tstamps] = times2rate(spktimes, 'binsize', binsize,...
    'winCalc', winCalc, 'c2r', true);

% calc fr according to states. note states, binsize, and spktimes must be
% the same units (typically sec)
if exist(fullfile(basepath, [basename '.AccuSleep_states.mat']))
    load(fullfile(basepath, [basename '.AccuSleep_states.mat']), 'ss')    
    fr.states.stateNames = ss.labelNames;
    nstates = length(ss.stateEpochs);
  
    % apply threshold for epoch leng to calc states
    thrBin = [60, 5, 5, 60, 8, 5];
    if length(thrBin) == 1
        thrBin = repmat(thrBin, 6, 1);
    elseif length(thrBin) ~= nstates - 1
        warning('thrBin length is different than the number of states')
    end

    % limit stateEpochs according to epoch length and fit to winCalc
    for istate = 1 : nstates - 1
        epochIdx = ss.stateEpochs{istate}(:, 2) < winCalc(2) &...
            ss.stateEpochs{istate}(:, 1) > winCalc(1);
        thrIdx =  ss.epLen{istate} > thrBin(istate);
        ss.stateEpochs{istate} = ss.stateEpochs{istate}(thrIdx & epochIdx, :);
         if ~isempty(ss.stateEpochs{istate})
            [fr.states.fr{istate}, fr.states.binedges, fr.states.tstamps{istate}, fr.states.binidx] =...
                times2rate(spktimes, 'binsize', binsize, 'winCalc', ss.stateEpochs{istate}, 'c2r', true);
        end
    end
        
    % buzsaki format
elseif exist(fullfile(basepath, [basename '.SleepState.states.mat']))
    load(fullfile(basepath, [basename '.SleepState.states.mat']))
    
    fr.states.statenames = {'WAKE', 'NREM', 'REM'};
    ss = struct2cell(SleepState.ints);
    nstates = length(ss);
    for i = 1 : nstates
        statetimes = ss{i};
        [fr.states.fr{i}, fr.states.binedges, fr.states.tstamps{i}, fr.states.binidx{i}] =...
            times2rate(spktimes, 'binsize', binsize, 'winCalc', statetimes, 'c2r', true);
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
fr.norm = fr.strd ./ mean(fr.strd(:, winBL(1) : winBL(2)), 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply criterions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bl = mean(fr.strd(:, winBL(1) : winBL(2)));

% select units who fired above thr
if any(strcmp(select, 'thr'))
    ithr = bl > 0.05;
else
    ithr = ones(nunits, 1);
end
% select units with low variability
if any(strcmp(select, 'stable'))
    bl_std = std(fr(:, win(1) : win(2)), [], 2);
    istable = bl_std < bl;
else
    istable = ones(nunits, 1);
end
idx = istable & ithr;

% add params
fr.winBL = winBL;
fr.binsize = binsize;
fr.mfr = mean(fr.strd, 2);

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
    if ischar(saveVar)
        save([basepath, filesep, basename, '.' saveVar '.mat'], 'fr')
    else
        save([basepath, filesep, basename, '.fr.mat'], 'fr')
    end
end

return

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% messaround with figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% session params
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', false, 'saveVar', false);
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

% recalculate firing rate
binsize = 60;
winBL = [1 * 60 120 * 60];
% winBL = [1 Inf];
fr = firingRate(spikes.times, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
    'binsize', binsize, 'saveVar', true, 'smet', 'MA', 'winBL', winBL);

% load vars
cd(basepath)
[~, basename] = fileparts(basepath);
load([basename '.cell_metrics.cellinfo.mat'])
load([basename '.spikes.cellinfo.mat'])
load([basename '.SleepState.states.mat'])
load([basename '.fr.mat'])
infoname = dir('*datInfo*');
load(infoname.name)

% states
states = {SleepState.ints.WAKEstate, SleepState.ints.NREMstate, SleepState.ints.REMstate};
for ii = 1 : length(states)
    tStates{ii} = InIntervals(fr.tstamps, states{ii});
    t{ii} = fr.tstamps(tStates{ii});
    frStates{ii} = mean(fr.strd(:, tStates{ii}), 2);
end

grp = [1 : 4];
pyr = strcmp(cell_metrics.putativeCellType, 'Pyramidal Cell');
int = strcmp(cell_metrics.putativeCellType, 'Narrow Interneuron');
% pyr = ones(1, length(spikes.ts)); % override
if isfield(spikes, 'su')
    su = spikes.su';
else
    su = ones(1, length(spikes.ts));
end
grpidx = zeros(1, length(spikes.shankID));
for ii = 1 : length(grp)
    grpidx = grpidx | spikes.shankID == grp(ii);
end

% fr vs. time -------------------------------------------------------------
% select units
units = int & su & grpidx;
data = fr.strd;
tstamps = fr.tstamps / 60;
nsamps = cumsum(datInfo.nsamps);
state = 2;      % 1 - awake; 2 - NREM

figure
plot(tstamps, (data(units, :))')
hold on
medata = median((data(units, :)), 1);
% plot(tstamps, medata, 'k', 'LineWidth', 5)
stdshade(data(units, :), 0.3, 'k', tstamps)
for i = 1 : length(nsamps) - 1
    plot([nsamps(i) nsamps(i)] / fs / 60, ylim, '--k')
end
axis tight
Y = ylim;
fill([states{state} fliplr(states{state})]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
    'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
% ylim([0 3])
ylabel('norm MFR')

% fr between two time periods----------------------------------------------

data = fr.strd;
state = 1;      % 1 - awake; 2 - NREM

win = floor([1, nsamps(2);...
    nsamps(4) nsamps(end)] / fs / 60);
if win(1, 1) == 0; win(1, 1) = 1; end

wint1 = InIntervals(tstamps, win(1, :));
wint2 = InIntervals(tstamps, win(2, :));
idx1 = wint1 & tStates{state};
idx2 = wint2 & tStates{state};

% average in periods
m1 = mean(data(units, idx1), 2)
m2 = mean(data(units, idx2), 2)

% histogram
figure
histogram(m1, 8)
hold on
histogram(m2, 8)

% scatter plot
units = su & grpidx;
m1 = mean(data(units, idx1), 2)
m2 = mean(data(units, idx2), 2)
c = repmat([0 0 1], length(spikes.ts), 1);
c(int, :) = repmat([1 0 0], sum(int), 1);

scatter(m1, m2, 50, c(units, :), 'filled')
hold on
X = xlim;
Y = ylim;
plot([0 X(2)], [0 Y(2)], '--k')


