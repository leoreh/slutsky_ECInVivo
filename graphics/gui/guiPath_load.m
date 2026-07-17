function cfgData = guiPath_load(cfgData, basepath, ctx)
% GUIPATH_LOAD Populate a curation cfgData: load each panel's data into its field.
%
%   cfgData = GUIPATH_LOAD(cfgData, basepath) reads the address (src) of every
%   panel in the flat cfgData and writes the loaded data back INTO that panel
%   (.data, .fs, resolved .ylim). A panel that ALREADY has data is left untouched,
%   so a cfgData carried over from another preset (or from a running guiPath)
%   is only topped up, not reloaded. A panel whose address fails to load is
%   dropped. The returned cfgData is the same struct, now "full".
%
%   basepath defaults to pwd. ctx (optional) is a guiPath_ctx file cache; by
%   default a fresh one is made so each file is read at most once per call.
%
%   This is the "addresses -> data" half of the design (see guiPath_doc). Behaviour
%   (mode / window / save) lives in cfgGui, set by guiPath_presets, not here.
%
%   DEPENDENCIES:
%       guiPath_src, guiPath_ctx, guiPath_panel; as_loadConfig, basepaths2vars,
%       binary_load, ripp_sigPrep (for computed 'fn:' sources).
%
%   See also guiPath_presets, guiPath_panel, guiPath_src, guiPath, guiPath_doc.
%
%   HISTORY:
%       Created: 05 Jul 2026 - the loader (was guiPath_load); loads data
%                              only, skips already-loaded panels, drops failures.

if nargin < 2 || isempty(basepath), basepath = pwd; end
if nargin < 3 || isempty(ctx)
    [~, bn] = fileparts(basepath);
    ctx = guiPath_ctx(basepath, bn);
end

fns = fieldnames(cfgData);
keep = true(1, numel(fns));
for i = 1:numel(fns)
    p = normPanel(cfgData.(fns{i}), fns{i});
    if isfield(p, 'data') && ~isempty(p.data)      % already loaded -> reuse, do not reload
        cfgData.(fns{i}) = p; continue;
    end
    try
        p = loadPanel(p, ctx);
    catch
        keep(i) = false; continue;                 % missing file / bad address -> drop
    end
    cfgData.(fns{i}) = p;
end
if any(~keep), cfgData = rmfield(cfgData, fns(~keep)); end
end

% =========================================================================
%  PANEL NORMALIZATION + LOAD
% =========================================================================

function p = normPanel(entry, fieldName)
% guarantee the fields a panel needs (a guiPath_panel entry already has them; a
% hand-written plain struct may not)
p = entry;
if ~isfield(p, 'type')   || isempty(p.type),   error('guiPath_load:type', 'panel "%s" needs a type', fieldName); end
if ~isfield(p, 'region') || isempty(p.region), p.region = 'bottom'; end
if ~isfield(p, 'src'),    p.src   = []; end
if ~isfield(p, 'name')   || isempty(p.name),   p.name = fieldName; end
if ~isfield(p, 'fs')     || isempty(p.fs),     p.fs = []; end
if ~isfield(p, 'height') || isempty(p.height), p.height = 1; end
if ~isfield(p, 'clr')    || isempty(p.clr),    p.clr = 'k'; end
if ~isfield(p, 'label'),  p.label = ''; end
if ~isfield(p, 'order')  || isempty(p.order),  p.order = 99; end
if ~isfield(p, 'ylim'),   p.ylim = []; end
if strcmp(p.type, 'trace') && isempty(p.ylim), p.ylim = 'prc'; end
end

function p = loadPanel(p, ctx)
% resolve the panel's address into .data (+ .fs, resolved .ylim)
switch p.type
    case 'trace'
        % one signal, as a double column. Only a bin: address that named several
        % channels is a matrix to average (its long-standing behaviour); any
        % other 2-D data is a single signal that happens to be a row - flatten
        % it, do not average across it (averaging a [1 x N] row was the bug that
        % turned emg_rms into one scalar).
        [val, meta] = resolveData(p.src, ctx);
        val = double(val);
        if isfield(meta, 'kind') && strcmp(meta.kind, 'bin') && size(val, 2) > 1
            val = mean(val, 2);
        end
        p.data = val(:);
        if isempty(p.fs), p.fs = meta.fs; end
        p.ylim = resolveYlim(p.ylim, p.data);

    case 'traces'
        % a vertical stack: keep every channel, in its native class (int16 for
        % an .lfp), so N channels of a long session stay affordable. The draw
        % cuts and casts only the window slice. chInfo holds the fixed display
        % stats (spacing, per-channel baseline, labels), computed once from the
        % whole signal so the stack neither shifts nor rescales between windows.
        [val, meta] = resolveData(p.src, ctx);
        p.data = val;                              % [nSamples x nCh], native
        if isempty(p.fs), p.fs = meta.fs; end
        p.chInfo = traceStack(val, meta);

    case 'spec'
        p.data = resolveData(p.src, ctx);          % adapter struct .s/.freq/.tstamps
        p.fs = NaN;

    case 'hypnogram'
        bt = resolveData(p.src, ctx);              % cell of [start end] in seconds
        p.data = toHoursCell(bt, ctx);             % length nstates, in hours
        p.fs = NaN;

    case 'raster'
        spk = resolveData(p.src, ctx);             % cell of spike-time vectors [s]
        if ~iscell(spk), spk = {spk(:)}; end
        p.data = spk(:); p.fs = NaN;

    case 'eventTicks'
        raw = resolveData(p.src, ctx);             % ed/ripp struct or Nx1/2/3 matrix
        p.data = eventsFrom(raw); p.fs = NaN;

    case 'stateStrip'
        p = loadStrip(p, ctx);

    otherwise
        error('guiPath_load:type', 'unknown panel type "%s"', p.type);
end
end

function p = loadStrip(p, ctx)
% curation target for states: labels + epoch centres + names/colours/nstates
tsp = guiPath_src('sleep_sig:spec_tstamps', ctx);   % epoch centres [s]
epochT = tsp(:); nEp = numel(epochT);
if nEp == 0, error('guiPath_load:noSpec', 'states need a spectrogram for epoch times'); end
asCfg = getAsCfg(ctx); nstates = asCfg.nstates;
labels = [];
if ~isempty(p.src)                                 % explicit labels source (Load flow)
    try
        v = resolveData(p.src, ctx);
        if isnumeric(v) && numel(v) == nEp, labels = double(v(:)); end
    catch
    end
end
if isempty(labels)                                 % resume manual scoring, else classifier
    labels = loadInitialLabels(ctx.basepath, ctx.basename, nEp, nstates);
end
p.data = struct('labels', labels(:), 'epochT', epochT, ...
    'names', {asCfg.names}, 'colors', {asCfg.colors}, 'nstates', nstates);
p.fs = NaN;
end

% =========================================================================
%  DATA RESOLUTION (address strings, computed 'fn:' sources)
% =========================================================================

function [val, meta] = resolveData(src, ctx)
% an address (via guiPath_src) or a computed 'fn:NAME[.field]' source
if ischar(src) && numel(src) >= 3 && strcmp(src(1:3), 'fn:')
    [nm, path] = splitFirst(src(4:end));
    r = computed(nm, ctx);
    val = resolvePath(r.val, path);
    meta = struct('fs', r.fs, 'kind', 'fn');
    if isfield(r, 'ch'), meta.ch = r.ch; end    % channel labels for a traces panel
else
    [val, meta] = guiPath_src(src, ctx);
end
end

function r = computed(name, ctx)
% cached computed resources that are not a plain file read
key = ['fn:' name];
if isKey(ctx.cache, key), r = ctx.cache(key); return; end
switch name
    case 'spec'
        guiPath_src('sleep_sig', ctx);              % warm the sleep_sig cache
        ss = ctx.cache('sleepsig');
        if isempty(ss.spec), error('guiPath_load:noSpec', 'no spectrogram in sleep_sig'); end
        r = struct('val', ss.spec, 'fs', NaN);

    case 'ripple'
        vses = basepaths2vars('basepaths', {ctx.basepath}, 'vars', {'session'}, 'flgPrnt', false);
        session = vses.session;
        nCh = session.extracellular.nChannels;
        fs  = session.extracellular.srLfp;
        if round(session.extracellular.sr) == 24414, bit2uv = 1; else, bit2uv = 0.195; end
        rippCh = evt_rippCh(ctx.basepath, ctx.basename, session);
        lfp = double(binary_load(fullfile(ctx.basepath, [ctx.basename, '.lfp']), 'duration', Inf, ...
            'fs', fs, 'nCh', nCh, 'start', 0, 'ch', rippCh, 'downsample', 1, 'bit2uv', bit2uv));
        if size(lfp, 2) > 1, lfp = mean(lfp, 2); end
        lfp = lfp(:);
        passband = [80 250];
        try
            rp = guiPath_src('ripp', ctx);
            if isfield(rp, 'info') && isfield(rp.info, 'passband') && ~isempty(rp.info.passband)
                passband = rp.info.passband;
            end
        catch
        end
        rs = ripp_sigPrep(lfp, fs, 'passband', passband, 'zMet', 'adaptive');
        r = struct('val', rs, 'fs', fs);

    case 'rippStack'
        % the channels to stack for ripple curation = the channel(s) the ripple
        % pipeline detected on (ripp.info.rippCh, via evt_rippCh), shown
        % separately rather than averaged as detection does, so their morphology
        % can be compared. Follows the ripple output, not the session tag.
        ch = evt_rippCh(ctx.basepath, ctx.basename);
        ch = ch(:)';
        [val, meta] = guiPath_src(sprintf('bin:%s', mat2str(ch)), ctx);
        r = struct('val', val, 'fs', meta.fs, 'ch', meta.ch);

    case 'edLfp'
        % the LFP the EDs preset shows = the channel ED detection ran on, read
        % from ed.info so the display always follows detection (see ed_wrapper).
        % 'lfp' detection -> that raw channel; 'eeg' (or a pre-edCh ed) -> sSig.eeg
        src = 'eeg'; ch = [];
        try
            edStruct = guiPath_src('ed', ctx);
            if isfield(edStruct, 'info')
                if isfield(edStruct.info, 'sigSource') && ~isempty(edStruct.info.sigSource)
                    src = edStruct.info.sigSource;
                end
                if isfield(edStruct.info, 'edCh'), ch = edStruct.info.edCh; end
            end
        catch
        end
        if strcmpi(src, 'lfp') && ~isempty(ch)
            [val, meta] = guiPath_src(sprintf('bin:%s', mat2str(ch(:)')), ctx);
        else
            [val, meta] = guiPath_src('sleep_sig:eeg', ctx);
        end
        r = struct('val', val, 'fs', meta.fs);

    otherwise
        error('guiPath_load:fn', 'unknown computed source "fn:%s"', name);
end
ctx.cache(key) = r;
end

% =========================================================================
%  SMALL BUILDERS / HELPERS (pure)
% =========================================================================

function yl = resolveYlim(spc, data)
% Resolve a panel's ylim spec against its data:
%       [lo hi]     absolute limits, taken as given
%       <scalar> p  percentile clip to [p, 100-p]; 0 <= p < 50, 0 = full range
%       'prc'       percentile clip at prcDflt
%       'full', []  autoscale
%
% The percentile is what controls how much of the axis a trace fills, and a
% wide clip is what makes a trace look thin: at 0.1 a lone artifact sets the
% range and the signal collapses toward the midline. Raising p trades a few
% clipped extremes for amplitude on everything else, so a preset that cares
% about waveform shape should pass its own p rather than take the default.
%
% The percentile is estimated on a subsample (~1e5 points) rather than the whole
% signal: on a full-session trace the clip is visually identical but the sort is
% ~orders of magnitude cheaper (a full-session prctile is ~1.5 s; this is ~5 ms).
% A degenerate range (flat / NaN signal, or p >= 50) is dropped so the axis
% autoscales.
prcDflt = 0.1;

if isnumeric(spc) && numel(spc) == 2, yl = spc; return; end
yl = [];

prc = [];
if ischar(spc) && strcmp(spc, 'prc'),  prc = prcDflt; end
if isnumeric(spc) && isscalar(spc),    prc = spc;     end
if isempty(prc) || isempty(data), return; end

n = numel(data);
if n > 2e5, s = double(data(1:ceil(n / 1e5):end)); else, s = double(data(:)); end
yl = prctile(s, [prc, 100 - prc]);
if numel(yl) ~= 2 || ~all(isfinite(yl)) || yl(2) <= yl(1), yl = []; end
end

function info = traceStack(val, meta)
% fixed display stats for a channel stack, from a subsample of the whole signal:
%   .spacing  vertical gap between channels (data units)
%   .base     per-channel baseline (median), subtracted so each channel is
%             centred on its own row regardless of DC offset
%   .labels   channel numbers for the y-axis
% Computed once and held for the panel's life so the stack neither drifts nor
% rescales as the window moves. Spacing is robust (MAD-based) so one large
% channel does not blow the stack apart; 6 SD clears ordinary excursions while
% keeping the channels close enough to compare.
nCh = size(val, 2);
nr  = size(val, 1);
if nr > 2e5, s = double(val(1:ceil(nr / 1e5):end, :)); else, s = double(val); end
base = median(s, 1);
sd = median(abs(s - base), 1) / 0.6745;            % per-channel robust SD
sp = 6 * median(sd(isfinite(sd)));
if isempty(sp) || ~isfinite(sp) || sp <= 0, sp = 1; end
if isfield(meta, 'ch') && numel(meta.ch) == nCh, labels = meta.ch(:)';
else,                                             labels = 1:nCh;
end
info = struct('spacing', sp, 'base', base(:)', 'labels', labels);
end

function boutHr = toHoursCell(bt, ctx)
% sleep-state bouts (seconds cell) -> hours cell of length nstates
asCfg = getAsCfg(ctx); ns = asCfg.nstates;
if ~iscell(bt), bt = {bt}; end
boutHr = cell(1, ns);
for s = 1:ns
    if s <= numel(bt) && ~isempty(bt{s}), boutHr{s} = bt{s} / 3600; else, boutHr{s} = zeros(0, 2); end
end
end

function ev = eventsFrom(val)
% compact events struct from a detection struct (ed/ripp) or an Nx1/2/3 matrix
if isstruct(val)
    if ~isfield(val, 'peakTime') || isempty(val.peakTime)
        error('guiPath_load:events', 'an events struct needs a non-empty peakTime field');
    end
    ev = struct('peakTime', val.peakTime(:));
    if isfield(val, 'times') && ~isempty(val.times),  ev.times = val.times; end
    if isfield(val, 'accepted') && ~isempty(val.accepted), ev.accepted = logical(val.accepted(:)); end
    if isfield(val, 'state') && ~isempty(val.state),   ev.state = val.state(:); end
    return;
end
M = double(val);
if ~isnumeric(M) || isempty(M), error('guiPath_load:events', 'events must be a struct or numeric matrix'); end
if isvector(M)
    ev = struct('peakTime', M(:));
elseif size(M, 2) == 3
    ev = struct('peakTime', M(:, 2), 'times', [M(:, 1), M(:, 3)]);
elseif size(M, 2) == 2
    ev = struct('peakTime', mean(M, 2), 'times', M);
else
    error('guiPath_load:events', 'events matrix must be Nx1, Nx2, or Nx3');
end
end

function labels = loadInitialLabels(basepath, basename, nEp, nstates)
% initial labels for scoring: manual (sleep_labelsMan) if present, else the
% classifier (sleep_states ss.labels, or sleep_labels), else all undefined. Used
% only if its length matches the epoch count.
man    = fullfile(basepath, [basename, '.sleep_labelsMan.mat']);
states = fullfile(basepath, [basename, '.sleep_states.mat']);
plain  = fullfile(basepath, [basename, '.sleep_labels.mat']);
labels = tryLabels(man, 'labels', nEp);
if isempty(labels), labels = tryLabels(states, 'ss', nEp); end
if isempty(labels), labels = tryLabels(plain, 'labels', nEp); end
if isempty(labels), labels = ones(nEp, 1) * (nstates + 1); end
end

function labels = tryLabels(file, var, nEp)
% pull a length-nEp label vector from a file: 'labels' variable, or ss.labels
labels = [];
if ~isfile(file), return; end
try
    s = load(file, var);
    if ~isfield(s, var), return; end
    if strcmp(var, 'ss')
        if isfield(s.ss, 'labels') && numel(s.ss.labels) == nEp, labels = double(s.ss.labels(:)); end
    else
        if numel(s.(var)) == nEp, labels = double(s.(var)(:)); end
    end
catch
end
end

function cfg = getAsCfg(ctx)
if isKey(ctx.cache, 'ascfg'), cfg = ctx.cache('ascfg'); return; end
cfg = as_loadConfig();
ctx.cache('ascfg') = cfg;
end

function [head, tail] = splitFirst(path)
if isempty(path), head = ''; tail = ''; return; end
di = find(path == '.', 1);
if isempty(di), head = path; tail = ''; else, head = path(1:di - 1); tail = path(di + 1:end); end
end

function val = resolvePath(base, path)
val = base;
if isempty(path), return; end
parts = strsplit(path, '.');
for i = 1:numel(parts), val = val.(parts{i}); end
end
