function [varMap, guiMap] = guiPath_presets(name, basepath, basename, ctx)

% Return a curation preset for guiPath: a varMap (data) and a guiMap (view).
%
% A preset is two halves. The varMap is a struct of recipes (see var_recipe) -
% WHAT to load; var_load fills it. The guiMap is the arrangement - a struct of
% view panels (see guiPath_panel) plus behaviour (.name .mode .win .save). The
% two are joined by name: a panel's .var picks a varMap entry.
%
% Presets are session-aware: the ripple / ED channels and passband are resolved
% here (evt_rippCh, ripp.info, ed.info) and frozen into plain recipes, so the
% loader stays generic. The states preset composes its editable strip inline.
%
% EXAMPLES
% - [varMap, guiMap] = guiPath_presets('ripp', basepath)
%   the ripple preset's recipes + panels (slim; pass varMap to var_load).
% - list = guiPath_presets()
%   the preset list, a struct array {name, file}, for the dropdown / auto-detect.
%
% INPUTS
% - name            <char>(opt) preset name (case-insensitive; file tokens work):
%                   'EDs'|'ed', 'Ripples'|'ripp', 'States'|'sleep_states',
%                   'template'|'Custom'. No arg -> the {name, file} list.
% - basepath        <char>(opt) session folder. Default pwd.
% - basename        <char>(opt) file stem. Default: folder name of basepath.
% - ctx             <struct>(opt) var_ctx cache to share with the later var_load.
%
% OUTPUTS
% - varMap          <struct> name -> recipe.
% - guiMap          <struct> .panels (view panels) + .name .mode .win .save.
% - list            <struct> {name, file} (the no-arg form).
%
% SEE ALSO
% - var_recipe, var_load, guiPath_panel, guiPath, guiPath_doc.
%
% HISTORY
% - 260719          split into varMap (data) + guiMap (view); session-aware
%                   recipes; was [cfgData, cfgGui] over the address grammar.


%% ========================================================================
%  ENUMERATION (no arg -> the preset list)
%  ========================================================================
if nargin < 1
    varMap = struct('name', {'EDs', 'Ripples', 'States'}, ...
                    'file', {'ed',  'ripp',    'sleep_states'});
    return
end

if nargin < 2 || isempty(basepath), basepath = pwd; end
if nargin < 3 || isempty(basename), [~, basename] = fileparts(basepath); end
if nargin < 4 || isempty(ctx),      ctx = var_ctx(basepath, basename); end


%% ========================================================================
%  BUILD PRESET
%  ========================================================================
switch lower(char(name))

    case {'eds', 'ed'}
        [varMap, panels] = cfg_ed(ctx);
        guiMap = behaviour(panels, 'EDs', 'events', 1, 'ed');

    case {'ripples', 'ripp'}
        [varMap, panels] = cfg_ripp(ctx);
        guiMap = behaviour(panels, 'Ripples', 'events', 0.4, 'ripp');

    case {'states', 'sleep_states'}
        [varMap, panels] = cfg_states(ctx);
        guiMap = behaviour(panels, 'States', 'states', 10, 'labelsMan');

    case {'template', 'custom', ''}
        varMap = struct();
        guiMap = behaviour(struct(), 'Custom', 'events', 1, '');

    otherwise
        error('guiPath_presets:name', ...
            'unknown preset "%s" (EDs | Ripples | States | template)', name);
end

end


% =========================================================================
%  PRESET BUILDERS
% =========================================================================

function [vm, gm] = cfg_ed(ctx)
% EDs: events from <basename>.ed.mat, shown on the sleep context.

vm = sleepContext();                              % hyp / spec / emgRms / emg
vm.ed     = var_recipe('matvar', 'file', 'ed');
vm.raster = var_recipe('matvar', 'file', 'spikes', 'var', 'spikes', 'path', 'times');

% the LFP the EDs preset shows = the channel ED detection ran on (from ed.info):
% an 'lfp' source -> that raw channel; 'eeg' (default) -> sSig.eeg
[edSrc, edCh] = edSigSource(ctx);
if strcmpi(edSrc, 'lfp') && ~isempty(edCh)
    vm.lfp = var_recipe('bin', 'file', 'lfp', 'ch', edCh(:)', ...
        'average', true, 'outClass', 'native');
else
    vm.lfp = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'eeg');
end

% Top: state, spectrogram, EMG RMS, event ticks. Bottom: LFP, EMG, raster, state
gm = struct();
gm.hypT   = guiPath_panel('hypnogram', 'top', 'hyp', 'label', 'State');
gm.spec   = guiPath_panel('spec', 'top', 'spec');
gm.emgRms = guiPath_panel('trace', 'top', 'emgRms', 'label', 'EMG RMS', ...
    'height', 0.7, 'ylim', 'full');
gm.evt    = guiPath_panel('eventTicks', 'top', 'ed', 'label', 'Events');
gm.lfp    = guiPath_panel('trace', 'bottom', 'lfp', 'label', 'LFP', 'height', 1.2);
gm.emg    = guiPath_panel('trace', 'bottom', 'emg', 'label', 'EMG', 'height', 0.8);
gm.raster = guiPath_panel('raster', 'bottom', 'raster', 'label', 'Units');
gm.hypB   = guiPath_panel('hypnogram', 'bottom', 'hyp', 'label', 'State');
end


function [vm, gm] = cfg_ripp(ctx)
% Ripples: events from <basename>.ripp.mat; ripple-band LFP in the window.

vm = sleepContext();                              % hyp / spec / emgRms / emg
vm.ripp   = var_recipe('matvar', 'file', 'ripp');
vm.raster = var_recipe('matvar', 'file', 'spikes', 'var', 'spikes', 'path', 'times');

% the shank around the ripple channel (raw, stacked) + the filtered detection
% channel. Both follow the channel ripple detection ran on (evt_rippCh); the
% filtered trace is ripp_sigPrep on the bit2uv-scaled average.
session = getSession(ctx);
if round(session.extracellular.sr) == 24414, b2u = 1; else, b2u = 0.195; end
pb = rippPassband(ctx);
vm.rippStack = var_recipe('bin', 'file', 'lfp', ...
    'ch', evt_rippCh(ctx.basepath, ctx.basename), 'average', false, 'outClass', 'native');
vm.rippFilt  = var_recipe('bin', 'file', 'lfp', ...
    'ch', evt_rippCh(ctx.basepath, ctx.basename, session), 'average', true, ...
    'bit2uv', b2u, 'transform', {'rippPrep', {pb}}, 'path', 'filt');

gm = struct();
gm.hypT      = guiPath_panel('hypnogram', 'top', 'hyp', 'label', 'State');
gm.spec      = guiPath_panel('spec', 'top', 'spec');
gm.emgRms    = guiPath_panel('trace', 'top', 'emgRms', 'label', 'EMG RMS', ...
    'height', 0.7, 'ylim', 'full');
gm.evt       = guiPath_panel('eventTicks', 'top', 'ripp', 'label', 'Ripples');
gm.rippStack = guiPath_panel('traces', 'bottom', 'rippStack', 'label', 'LFP', ...
    'height', 2.0);
gm.rippFilt  = guiPath_panel('trace', 'bottom', 'rippFilt', 'label', 'Filtered LFP', ...
    'height', 1.0);
gm.emg       = guiPath_panel('trace', 'bottom', 'emg', 'label', 'EMG', ...
    'height', 0.8, 'ylim', 1);
gm.raster    = guiPath_panel('raster', 'bottom', 'raster', 'label', 'Units');
end


function [vm, gm] = cfg_states(ctx)
% States: the editable stateStrip is the target; top and bottom share it.

vm = struct();
vm.spec   = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', {'spec', 'spec_freq', 'spec_tstamps'});
vm.emgRms = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg_rms', 'fs', 1);
vm.eeg    = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'eeg');
vm.emg    = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg');
vm.states = var_recipe('value', 'data', stateStripData(ctx));

gm = struct();
gm.stripT = guiPath_panel('stateStrip', 'top', 'states', 'label', 'State');
gm.spec   = guiPath_panel('spec', 'top', 'spec');
gm.emgRms = guiPath_panel('trace', 'top', 'emgRms', 'label', 'EMG RMS', ...
    'height', 0.7, 'ylim', 'full');
gm.eeg    = guiPath_panel('trace', 'bottom', 'eeg', 'label', 'LFP', 'height', 1.2);
gm.emg    = guiPath_panel('trace', 'bottom', 'emg', 'label', 'EMG', 'height', 0.8);
gm.stripB = guiPath_panel('stateStrip', 'bottom', 'states', 'label', 'State');
end


% =========================================================================
%  SHARED RECIPES + BEHAVIOUR
% =========================================================================

function vm = sleepContext()
% the sleep-context signals every modality shows on top: hypnogram, spectrogram,
% EMG RMS, EMG. (The spec is read from the assembled sleep_sig, not recomputed.)
vm = struct();
vm.hyp    = var_recipe('matvar', 'file', 'sleep_states', 'var', 'ss', ...
    'path', 'bouts.times');
vm.spec   = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', {'spec', 'spec_freq', 'spec_tstamps'});
vm.emgRms = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg_rms', 'fs', 1);
vm.emg    = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg');
end


function guiMap = behaviour(panels, name, mode, win, save)
guiMap = struct('panels', panels, 'name', name, 'mode', mode, ...
    'win', win, 'save', save);
end


% =========================================================================
%  SESSION-AWARE RESOLVERS (read via the shared ctx)
% =========================================================================

function session = getSession(ctx)
session = var_fetch(var_recipe('matvar', 'file', 'session'), ctx);
end


function pb = rippPassband(ctx)
% ripp.info.passband if present, else the default detection band
pb = [80 250];
try
    pb = var_fetch(var_recipe('matvar', 'file', 'ripp', 'var', 'ripp', ...
        'path', 'info.passband'), ctx);
    if isempty(pb), pb = [80 250]; end
catch
end
end


function [src, ch] = edSigSource(ctx)
% ed.info.sigSource / edCh (which channel ED detection ran on)
src = 'eeg'; ch = [];
try
    ed = var_fetch(var_recipe('matvar', 'file', 'ed'), ctx);
    if isfield(ed, 'info')
        if isfield(ed.info, 'sigSource') && ~isempty(ed.info.sigSource)
            src = ed.info.sigSource;
        end
        if isfield(ed.info, 'edCh'), ch = ed.info.edCh; end
    end
catch
end
end


function strip = stateStripData(ctx)
% compose the editable state strip: epoch centres (spec_tstamps), the state
% config (as_loadConfig), and an initial label vector (manual scoring if present,
% else the classifier, else all-undefined)
tsp = var_fetch(var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', 'spec_tstamps'), ctx);
epochT = tsp(:);
nEp = numel(epochT);
if nEp == 0
    error('guiPath_presets:noSpec', 'states need a spectrogram for epoch times');
end
asCfg = as_loadConfig();
nstates = asCfg.nstates;
labels = loadInitialLabels(ctx.basepath, ctx.basename, nEp, nstates);
strip = struct('labels', labels(:), 'epochT', epochT, ...
    'names', {asCfg.names}, 'colors', {asCfg.colors}, 'nstates', nstates);
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
        if isfield(s.ss, 'labels') && numel(s.ss.labels) == nEp
            labels = double(s.ss.labels(:));
        end
    else
        if numel(s.(var)) == nEp, labels = double(s.(var)(:)); end
    end
catch
end
end

% EOF
