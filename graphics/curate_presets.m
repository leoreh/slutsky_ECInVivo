function presets = curate_presets()
% CURATE_PRESETS Registry of curation presets for gui_curate.
%
%   presets = CURATE_PRESETS() returns a struct array, one entry per modality:
%       .name   (Char) display name shown in the Preset dropdown.
%       .var    (Char) file token used both to auto-detect the preset
%               (<basename>.<var>.mat) and as the curation target.
%       .load   (Fcn)  @(basepath, basename, pool) -> config struct with fields:
%                   .inputs   (Struct array) gui_curate inputs (signals + one
%                             eventTicks input). See gui_curate for the schema.
%                   .panels   (Struct array) default layout, each .source/.region.
%                   .saveFcn  (Fcn) @(accepted) persist the curation mask.
%                   .winPlot  (Num) default Bottom-window full width [s].
%
%   POOL: gui_curate keeps a persistent pool of already-loaded SIGNAL inputs and
%   passes it in. A loader reuses any signal already in the pool (poolGet) and
%   only loads what is missing, so revisiting a preset is instant and the heavy
%   signal prep (binary read, filtering) runs once. Signal names are unambiguous
%   across modalities (eeg vs rippleLfp) so the pool can hold both at once; the
%   shared 'raster' (spikes) is loaded once and reused.
%
%   To add a new default view (e.g. Sleep States), add one entry here with its
%   loader; gui_curate stays unchanged.
%
%   DEPENDENCIES:
%       ed_sigLoad, ripp_sigPrep, binary_load, basepaths2vars, as_loadConfig.
%
%   HISTORY:
%       Created: 23 Jun 2026 - presets extracted from ed_gui/ripp_curate.
%       Updated: 23 Jun 2026 - signal pool reuse; unambiguous signal names.

presets = struct('name', {}, 'var', {}, 'load', {});
presets(end + 1) = struct('name', 'EDs',     'var', 'ed',   'load', @presetEd);
presets(end + 1) = struct('name', 'Ripples', 'var', 'ripp', 'load', @presetRipp);
end

%% ========================================================================
%  PRESET: EDs (electrographic discharges)
%  ========================================================================
function cfg = presetEd(basepath, basename, pool)
% events from <basename>.ed.mat; signals from <basename>.sleep_sig.mat
edFile = fullfile(basepath, [basename, '.ed.mat']);
if ~isfile(edFile), error('curate_presets:noEd', 'missing %s', edFile); end
s = load(edFile, 'ed'); ed = s.ed;

I = {};

% LFP (eeg) / EMG / EMG RMS / spectrogram: reuse from pool, else load once
if hasInPool(pool, 'eeg')
    I{end+1} = poolGet(pool, 'eeg');
    I{end+1} = poolGet(pool, 'emg');
    if hasInPool(pool, 'emgRms'), I{end+1} = poolGet(pool, 'emgRms'); end
    if hasInPool(pool, 'spec'),   I{end+1} = poolGet(pool, 'spec');   end
else
    [~, ~, ~, ~, specAdapter, sSig] = ed_sigLoad(basepath, 'basename', basename);
    if isempty(specAdapter) && isfield(sSig, 'spec')
        specAdapter = struct('s', sSig.spec, 'freq', sSig.spec_freq(:), 'tstamps', sSig.spec_tstamps(:));
    end
    fs = ed.info.fs;
    I{end+1} = mkInput('eeg', 'trace', 'LFP', 1.2, sSig.eeg(:), fs, prc(sSig.eeg), 'k', 'narrow', 50);
    I{end+1} = mkInput('emg', 'trace', 'EMG', 0.8, sSig.emg(:), fs, prc(sSig.emg), 'k', 'narrow', 51);
    if isfield(sSig, 'emg_rms') && ~isempty(sSig.emg_rms)
        I{end+1} = mkInput('emgRms', 'trace', 'EMG RMS', 0.7, sSig.emg_rms(:), 1, [], 'k', 'wide', 30);
    end
    if ~isempty(specAdapter)
        I{end+1} = mkInput('spec', 'spec', 'Freq (Hz)', 1.4, specAdapter, NaN, [], [], 'wide', 20);
    end
end

% hypnogram (states) and raster (spikes): reuse from pool, else load once
if hasInPool(pool, 'hypnogram')
    I{end+1} = poolGet(pool, 'hypnogram');
else
    boutHr = loadStatesHr(basepath);
    if ~isempty(boutHr), I{end+1} = mkInput('hypnogram', 'hypnogram', 'State', 0.28, boutHr, NaN, [], [], 'wide', 10); end
end
I = appendRaster(I, pool, basepath);

% events (always fresh; cheap)
I{end+1} = mkInput('eventTicks', 'eventTicks', 'Events', 0.28, eventStruct(ed), NaN, [], [], 'wide', 40);
cfg.inputs = [I{:}];

% default layout (availability-filtered downstream)
cfg.panels = [pnl('hypnogram','wide'), pnl('spec','wide'), pnl('emgRms','wide'), pnl('eventTicks','wide'), ...
    pnl('eeg','narrow'), pnl('emg','narrow'), pnl('raster','narrow'), pnl('hypnogram','narrow')];
cfg.saveFcn = @(acc) saveAccepted(edFile, 'ed', acc);
cfg.winPlot = 1.0;
end

%% ========================================================================
%  PRESET: Ripples (SWR)
%  ========================================================================
function cfg = presetRipp(basepath, basename, pool)
% events from <basename>.ripp.mat; ripple-channel LFP recomputed from the binary
rippFile = fullfile(basepath, [basename, '.ripp.mat']);
if ~isfile(rippFile), error('curate_presets:noRipp', 'missing %s', rippFile); end
s = load(rippFile, 'ripp'); ripp = s.ripp;

I = {};
if hasInPool(pool, 'rippleLfp') && hasInPool(pool, 'rippleFilt')
    I{end+1} = poolGet(pool, 'rippleLfp');
    I{end+1} = poolGet(pool, 'rippleFilt');
else
    vses = basepaths2vars('basepaths', {basepath}, 'vars', {'session'}, 'flgPrnt', false);
    session = vses.session;
    nCh = session.extracellular.nChannels;
    fs  = session.extracellular.srLfp;
    if round(session.extracellular.sr) == 24414, bit2uv = 1; else, bit2uv = 0.195; end
    if isfield(session, 'channelTags') && isfield(session.channelTags, 'Ripple')
        rippCh = session.channelTags.Ripple;
    else
        rippCh = 1;
    end
    lfp = double(binary_load(fullfile(basepath, [basename, '.lfp']), 'duration', Inf, ...
        'fs', fs, 'nCh', nCh, 'start', 0, 'ch', rippCh, 'downsample', 1, 'bit2uv', bit2uv));
    if size(lfp, 2) > 1, lfp = mean(lfp, 2); end
    lfp = lfp(:);
    passband = [80 250];
    if isfield(ripp, 'info') && isfield(ripp.info, 'passband') && ~isempty(ripp.info.passband)
        passband = ripp.info.passband;
    end
    rs = ripp_sigPrep(lfp, fs, 'passband', passband, 'zMet', 'adaptive');
    I{end+1} = mkInput('rippleLfp',  'trace', 'LFP', 1.2, rs.lfp, fs, prc(rs.lfp), 'k', 'narrow', 50);
    I{end+1} = mkInput('rippleFilt', 'trace', 'Filtered LFP', 1.0, rs.filt, fs, prc(rs.filt), 'k', 'narrow', 51);
end
I = appendRaster(I, pool, basepath);

I{end+1} = mkInput('eventTicks', 'eventTicks', 'Ripples', 0.28, eventStruct(ripp), NaN, [], [], 'wide', 40);
cfg.inputs = [I{:}];

cfg.panels = [pnl('eventTicks','wide'), pnl('rippleLfp','narrow'), pnl('rippleFilt','narrow'), pnl('raster','narrow')];
cfg.saveFcn = @(acc) saveAccepted(rippFile, 'ripp', acc);
cfg.winPlot = 0.4;
end

%% ========================================================================
%  SHARED HELPERS
%  ========================================================================

function I = appendRaster(I, pool, basepath)
% raster (spikes) is shared across presets: reuse from pool, else load once
if hasInPool(pool, 'raster')
    I{end+1} = poolGet(pool, 'raster');
else
    spk = loadSpikes(basepath);
    if ~isempty(spk), I{end+1} = mkInput('raster', 'raster', 'Units', 1.2, spk, NaN, [], [], 'narrow', 60); end
end
end

function tf = hasInPool(pool, name)
tf = ~isempty(pool) && any(strcmp(name, {pool.name}));
end

function inp = poolGet(pool, name)
inp = pool(find(strcmp(name, {pool.name}), 1));
end

function inp = mkInput(name, type, label, height, data, fs, yl, clr, defRegion, defOrder)
% single-wrap data so a cell payload (raster/hypnogram) lands in one field
inp = struct('name', name, 'type', type, 'data', {data}, 'fs', fs, ...
    'ylim', yl, 'clr', clr, 'label', label, 'height', height, ...
    'defRegion', defRegion, 'defOrder', defOrder);
end

function p = pnl(source, region)
p = struct('source', source, 'region', region);
end

function y = prc(x)
y = prctile(x(:), [0.1, 99.9]);
end

function ev = eventStruct(s)
% compact generic events struct from a detection struct (ed or ripp): peak time
% plus optional start/stop, accepted seed, and vigilance state (for the status)
ev = struct('peakTime', s.peakTime(:));
if isfield(s, 'times') && ~isempty(s.times),  ev.times = s.times; end
if isfield(s, 'accepted') && ~isempty(s.accepted), ev.accepted = logical(s.accepted(:)); end
if isfield(s, 'state') && ~isempty(s.state),   ev.state = s.state(:); end
end

function spk = loadSpikes(basepath)
spk = {};
try
    v = basepaths2vars('basepaths', {basepath}, 'vars', {'spikes'}, 'flgPrnt', false);
    if isfield(v, 'spikes') && isfield(v.spikes, 'times'), spk = v.spikes.times(:); end
catch
end
end

function boutHr = loadStatesHr(basepath)
% sleep-state bouts as an hours cell of length cfg.nstates (or [] if absent)
boutHr = {};
try
    v = basepaths2vars('basepaths', {basepath}, 'vars', {'sleep_states'}, 'flgPrnt', false);
    if isfield(v, 'ss') && isfield(v.ss, 'bouts') && isfield(v.ss.bouts, 'times')
        bt  = v.ss.bouts.times;
        cfg = as_loadConfig();
        ns  = cfg.nstates;
        boutHr = cell(1, ns);
        for s = 1:ns
            if s <= numel(bt) && ~isempty(bt{s}), boutHr{s} = bt{s} / 3600; else, boutHr{s} = zeros(0, 2); end
        end
    end
catch
end
end

function saveAccepted(file, varName, accepted)
% set <var>.accepted and write all variables back, preserving file contents
S = load(file);
S.(varName).accepted = logical(accepted(:));
save(file, '-struct', 'S', '-v7.3');
end
