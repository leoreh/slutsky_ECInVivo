function [varMap, guiMap] = preset_sleep_states(ctx)

% guiPath preset: manual sleep scoring on the AccuSleep signals.
%
% varMap = what to load (var_recipe); guiMap = how to show it (guiPath_panel).
% guiPath joins them by name: a panel's .var picks a varMap entry, so a top +
% bottom pair over one var loads once and draws twice. See guiPath_doc for the
% format and guiPath_preset for how presets are found, loaded and saved.
%
% The state set (stateSet) is the same one every preset shows, so scoring can be
% done from any of them. What this preset adds is the signals you score AGAINST -
% the raw EEG and EMG in the window - and a mode that opens on the strip.
%
% INPUTS
% - ctx             <struct> var_ctx: basepath, basename, shared file cache.
%
% OUTPUTS
% - varMap          <struct> name -> recipe.
% - guiMap          <struct> .panels + .mode .win.
%
% SEE ALSO
% - guiPath_preset, var_recipe, guiPath_panel, guiPath_doc, stateSet.
%
% HISTORY
% - 260719          split out of guiPath_presets (one file per preset).
% - 260720          the strip comes from the shared stateSet (was composed
%                   here, from sleep_labelsMan alone).


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
varMap.states = var_recipe('value', 'data', stateSet(ctx));


%% ========================================================================
%  VIEW (guiMap)
%  ========================================================================
% the editable strip is the target; top and bottom share it (one loaded input)
guiMap = struct('panels', struct(), 'mode', 'states', 'win', 10);
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

% EOF
