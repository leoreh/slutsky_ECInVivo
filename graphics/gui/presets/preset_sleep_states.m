function [varMap, guiMap] = preset_sleep_states(ctx)

% guiPath preset: manual sleep scoring on the AccuSleep signals.
%
% varMap = what to load (var_recipe); guiMap = how to show it (guiPath_panel).
% guiPath joins them by name: a panel's .var picks a varMap entry, so a top +
% bottom pair over one var loads once and draws twice. See guiPath_doc for the
% format and guiPath_preset for how presets are found, loaded and saved.
%
% Unlike the event presets, the curated target here is composed rather than
% read: the editable strip pairs the spectrogram's epoch centres with an
% initial label vector (manual scoring if present, else the classifier).
%
% INPUTS
% - ctx             <struct> var_ctx: basepath, basename, shared file cache.
%
% OUTPUTS
% - varMap          <struct> name -> recipe.
% - guiMap          <struct> .panels + .mode .win .save.
%
% SEE ALSO
% - guiPath_preset, var_recipe, guiPath_panel, guiPath_doc, as_loadConfig.
%
% HISTORY
% - 260719          split out of guiPath_presets (one file per preset).


%% ========================================================================
%  DATA (varMap)
%  ========================================================================
% the AccuSleep signals, then the editable strip (already materialized)
varMap = struct();
varMap.spec   = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', {'spec', 'spec_freq', 'spec_tstamps'});
varMap.emgRms = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', 'emg_rms', 'fs', 1);
varMap.eeg    = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'eeg');
varMap.emg    = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg');
varMap.states = var_recipe('value', 'data', stateStripData(ctx));


%% ========================================================================
%  VIEW (guiMap)
%  ========================================================================
% the editable strip is the target; top and bottom share it (one loaded input)
guiMap = struct('panels', struct(), 'mode', 'states', 'win', 10, ...
    'save', 'labelsMan');
guiMap.panels.stripT = guiPath_panel('stateStrip', 'top', 'states', ...
    'label', 'State');
guiMap.panels.spec   = guiPath_panel('spec', 'top', 'spec');
guiMap.panels.emgRms = guiPath_panel('trace', 'top', 'emgRms', ...
    'label', 'EMG RMS', 'height', 0.7, 'ylim', 'full');
guiMap.panels.eeg    = guiPath_panel('trace', 'bottom', 'eeg', ...
    'label', 'LFP', 'height', 1.2);
guiMap.panels.emg    = guiPath_panel('trace', 'bottom', 'emg', ...
    'label', 'EMG', 'height', 0.8);
guiMap.panels.stripB = guiPath_panel('stateStrip', 'bottom', 'states', ...
    'label', 'State');

end


% =========================================================================
%  SESSION RESOLVERS (read via the shared ctx cache)
% =========================================================================

function strip = stateStripData(ctx)
% compose the editable state strip: epoch centres (spec_tstamps), the state
% config (as_loadConfig), and an initial label vector
tsp = var_fetch(var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', 'spec_tstamps'), ctx);
epochT = tsp(:);
nEp = numel(epochT);
if nEp == 0
    error('preset_sleep_states:noSpec', ...
        'states need a spectrogram for epoch times');
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
