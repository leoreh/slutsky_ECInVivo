function plot_fr_boutTimes(fr, varargin)

% plots state-dependent firing rates per bout across units and per unit
% across bouts

% INPUT
%   fr              struct. see calc_fr.m
%   basepath        recording session {pwd}
%   unitIdx         logical vec representing indices of units to include
%   sstates         numeric vec of selected state indices [1, 4, 5].
%                   assumes fr was calculated for these states only
%   saveFig         logical {true}
%
% 16 apr 24 LH 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'unitIdx', []);
addOptional(p, 'sstates', [1, 4, 5], @isnumeric);
addOptional(p, 'saveFig', false, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
unitIdx     = p.Results.unitIdx;
sstates     = p.Results.sstates;
saveFig     = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% state params
cfg = as_loadConfig;
clr = cfg.colors(sstates);
snames = cfg.names(sstates);

% graphics params
yLimit_fr = [0, 3];

% file
[~, basename] = fileparts(basepath);

% units
if isempty(unitIdx)
    unitIdx = 1 : size(fr.strd, 1);
end

% prep data
boutMat = cellfun(@(x) mean(x(unitIdx, :), 1)', fr.states.fr, 'uni', false);
boutMat = cell2padmat(boutMat, 2);
unitMat = fr.states.mfr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open figure
setMatlabGraphics(true)
fh = figure;
tlayout = [2, 3];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';

% comparison of bout MFR between states
axh = nexttile(th, 1, [1, 1]); cla; hold on
plot_boxMean('dataMat', boutMat, 'allpnts', true, 'axh', axh,...
    'plotType', 'bar', 'clr', vertcat(clr{:}))
ylim(yLimit_fr)
xticklabels(snames)
ylabel('Bout MFR [Hz]')
title(axh, 'All bouts')
title(th, basename, 'Interpreter', 'none')

% comparison of bout MFR between states, subsampled to the number of units
axh = nexttile(th, 2, [1, 1]); cla; hold on
clear tmp
tmp = cell(length(sstates), 1);
for istate = 1 : length(sstates)
    nbouts = sum(~isnan(boutMat(:, istate)));
    subIdx = randperm(nbouts, min([sum(unitIdx), nbouts]));
    tmp{istate} = boutMat(subIdx, istate);
end
tmp = cell2padmat(tmp, 2);
plot_boxMean('dataMat', tmp, 'allpnts', true, 'axh', axh,...
    'plotType', 'bar', 'clr', vertcat(clr{:}))
ylim(yLimit_fr)
xticklabels(snames)
ylabel('Bout MFR [Hz]')
title(axh, 'Random subset of bouts')

% comparison of unit MFR between states
axh = nexttile(th, 3, [1, 1]);
tmp = unitMat(unitIdx, :);
plot_boxMean('dataMat', tmp, 'allpnts', true, 'axh', axh,...
    'plotType', 'bar', 'clr', vertcat(clr{:}))
ylim(yLimit_fr)
xticklabels(snames)
ylabel('SU MFR [Hz]')
title(axh, 'Single Unit MFR')

% MFR per state bout across time
axh = nexttile(th, 4, [1, 2]); cla; hold on
for istate = 1 : length(sstates)
    tstamps = fr.states.tstamps{istate} / 60 / 60;
    yval = boutMat(:, istate);
    yval(isnan(yval)) = [];
    scatter(tstamps, yval, 30, 'filled', 'MarkerFaceColor', clr{istate})
    xlabel('Time [h]')
    ylabel('MFR [Hz]')
    legend(snames)
end
ylim(yLimit_fr)
title(axh, 'Bout MFR across time')

% correlation between bout dur and mfr
axh = nexttile(th, 6, [1, 1]); cla; hold on
for istate = 1 : length(sstates)
    boutDur = cellfun(@(x) diff(x), fr.states.binedges{istate}, 'uni', true);
    yval = boutMat(:, istate);
    yval(isnan(yval)) = [];
    scatter(boutDur, yval, 30, 'filled', 'MarkerFaceColor', clr{istate})
    xlabel('Bout Duration [s]')
    ylabel('MFR [Hz]')
    set(gca, 'xscale', 'log')
end
ylim(yLimit_fr)
xLimit = xlim;
xlim([1, max(xlim)])
title(axh, 'Correlation')

% save figure
if saveFig
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, [basename, '_', 'stateFRs']);
    savefig(fh, figname, 'compact')
end

% EOF