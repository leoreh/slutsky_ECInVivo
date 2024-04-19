function otl = get_specOutliers(varargin)

% gets spectrogram outliers, typically to exclude them as artifacts.
% spectrogram can be calculated or loaded from sSig.

% INPUT
%   basepath    char. fullpath to recording folder {pwd} ch      
%   ch          numeric vector depicting channels over which
%               the spectrogram should be calculated
%   stdThr      threshold in stds for defining outliers {5}
%   graphics    logical. plot figure {false}
%   saveVar     logical. organize and save struct {true}
%   flgCalc     logical. recalculate spectrogram or load {true} from sSig
%   flgForce    logical. force analysis or load if exists already {false}
%
% OUTPUT
%   otl         struct
%
% TO DO LIST
%
% 26 mar 24 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'flgCalc', false, @islogical)
addParameter(p, 'flgForce', false, @islogical)
addParameter(p, 'ch', [], @isnumeric)
addParameter(p, 'stdThr', 5, @isnumeric)
addParameter(p, 'graphics', false, @islogical)
addParameter(p, 'saveVar', true, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
flgCalc         = p.Results.flgCalc;
flgForce        = p.Results.flgForce;
ch              = p.Results.ch;
stdThr          = p.Results.stdThr;
saveVar         = p.Results.saveVar;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
[~, basename] = fileparts(basepath);
specfile = fullfile(basepath, [basename, '.spec.mat']);
ssigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
sessionfile = fullfile(basepath, [basename, '.session.mat']);
otlfile = fullfile(basepath, [basename, '.specOtl.mat']);

% load if exists
if ~flgForce & exist(otlfile)
    load(otlfile, 'otl')
    return
end

% channels for calculating spectogram
if isempty(ch)
    load(sessionfile, 'session')
    spkgrp = session.extracellular.spikeGroups.channels;
    ch = [spkgrp{:}];
end

% load tstamps
load(ssigfile, 'spec_tstamps')
tstamps = spec_tstamps;

% load / re-calc spec
if flgCalc
    spec = calc_spec('sig', [], 'fs', 1250, 'graphics', false, 'saveVar', false,...
    'padfft', 0, 'winstep', 1, 'logfreq', true, 'ftarget', [],...
    'ch', {ch}, 'force', true);
    s = spec.s;
    freq = spec.freq;
else
    load(ssigfile, 'spec')
    load(ssigfile, 'spec_freq')
    s = spec;
    freq = spec_freq;
    
    % organize spec struct for plotting
    clear spec
    spec.s = s; spec.freq = freq; spec.tstamps = tstamps;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detection of spec outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab only a specific frequency range to consider
[~, fidx1] = min(abs(freq - 15));
[~, fidx2] = min(abs(freq - 47));
fidx = fidx1 : fidx2;

[~, pc, ~, ~, expl] = pca(s(:, fidx), 'NumComponents', 3);
% dists = mahal(pc, pc);
dists = pc(:, 1);
specBoolean = dists > mean(dists) + std(dists) * stdThr;

% organize and save
otl.boolean = specBoolean;
otl.idx = find(otl.boolean);
otl.epochs = binary2epochs('vec', otl.boolean, 'minDur', [], 'maxDur', [],...
    'interDur', [], 'exclude', false, 'printFlag', false);
otl.tstamps = minutes(seconds(otl.idx));
otl.info.flgCalc = flgCalc;
otl.info.flgCalc = datetime(now, 'ConvertFrom', 'datenum');

if saveVar
    save(otlfile, 'otl')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(true)
fh = figure;
tlayout = [1, 3];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, basename, 'interpreter', 'none')
set(fh, 'DefaultAxesFontSize', 20);

% spectrogram marked with outliers
axh(1) = nexttile(th, 1, [1, 2]);
plot_spec(spec, 'ch', 1, 'logfreq', false, 'saveFig', false,...
    'axh', axh(1))
yLimit = ylim;
hold on
scatter(otl.idx / 3600, ones(length(otl.idx), 1) * yLimit(2), 20, 'filled', 'r')

% emg_rms
% axh(2) = nexttile(th, 4, [1, 2]);
% plot(tstamps / 3600, emg_rms, 'k', 'LineWidth', 0.5)
% 
% linkaxes(axh, 'x')
% axis tight

% spec projection in feature space
axh(3) = nexttile(th, 3, [1, 1]); cla; hold on
sh = scatter3(pc(:, 1), pc(:, 2), pc(:, 3),...
    30, 'k', 'filled');
scatter3(pc(otl.boolean, 1), pc(otl.boolean, 2), pc(otl.boolean, 3),...
    50, 'd', 'LineWidth', 1, 'MarkerEdgeColor', 'r');
view(-25, 25)
grid minor
sh.MarkerFaceAlpha = 0.6;
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
title(axh(3), 'Epoch PSD projection')

% EOF

