function artifacts = get_movArtifacts(varargin)

% gets movement artifacts from the spectrogram and validates with emg rms

% INPUT
%   basepath    char. fullpath to recording folder {pwd} ch      
%   ch          numeric vector depicting a groups of channels over which
%               the spectrogram should be calculated
%   graphics    logical. plot figure {false}
%   saveVar     logical. organize and save struct {true}
%   flgCalc     logical. recalculate spectrogram or load {true} from sSig
%   flgForce    logical. force analysis or load if exists already {false}
%
% OUTPUT
%   artIdx      numeric vectors if indices to seconds detected as artifacts
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
addParameter(p, 'graphics', false, @islogical)
addParameter(p, 'saveVar', true, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
flgCalc         = p.Results.flgCalc;
flgForce        = p.Results.flgForce;
ch              = p.Results.ch;
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
artfile = fullfile(basepath, [basename, '.artifacts.mat']);

% load if exists
if ~flgForce & exist(artfile)
    load(artfile, 'artifacts')
    return
end

% channels for calculating spectogram
if isempty(ch)
    load(sessionfile, 'session')
    spkgrp = session.extracellular.spikeGroups.channels;
    ch = [spkgrp{:}];
end

% thresholds for defining outliers. artifacts must follow both conditions
prctSpec = 99.9;
prctEmg = 0;

% exponent estimate of power low
powExp = 2;

% load data
load(ssigfile, 'emg_rms')
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
[~, fidx2] = min(abs(freq - 49));
fidx = fidx1 : fidx2;

% weights for average based on simple frequency power law
w = zeros(length(freq), 1);
w(fidx) = freq(fidx) .^ powExp;
w = w / sum(w);

% weighted average of spec
wpow = sum(s .* w', 2);

% get outliers
specOutliers = wpow > prctile(wpow, prctSpec);

% confirm outliers occur at times of high emg
emgOutliers = emg_rms' > prctile(emg_rms, prctEmg);

% organize and save
artifacts.boolean = specOutliers & emgOutliers;
artifacts.idx = find(artifacts.boolean);
artifacts.bouts = binary2bouts('vec', artifacts.boolean, 'minDur', [], 'maxDur', [],...
    'interDur', [], 'exclude', false, 'printFlag', false);
artifacts.tstamps = minutes(seconds(artifacts.idx));
artifacts.info.flgCalc = flgCalc;
artifacts.info.flgCalc = datetime(now, 'ConvertFrom', 'datenum');

if saveVar
    save(artfile, 'artifacts')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(true)
fh = figure;
tlayout = [2, 1];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, basename, 'interpreter', 'none')

axh(1) = nexttile(th, 1, [1, 1]);
plot_spec(spec, 'ch', 1, 'logfreq', false, 'saveFig', false,...
    'axh', axh(1))
yLimit = ylim;
hold on
scatter(artifacts.idx / 3600, ones(length(artifacts.idx), 1) * yLimit(2), 20, 'filled', 'r')

axh(2) = nexttile(th, 2, [1, 1]);
plot(tstamps / 3600, emg_rms, 'k', 'LineWidth', 0.5)

linkaxes(axh, 'x')
axis tight

% EOF

