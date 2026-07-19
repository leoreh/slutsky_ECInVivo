function [varMap, guiMap] = preset_ed(ctx)

% guiPath preset: epileptiform discharges from <basename>.ed.mat.
%
% varMap = what to load (var_recipe); guiMap = how to show it (guiPath_panel).
% guiPath joins them by name: a panel's .var picks a varMap entry, so a top +
% bottom pair over one var loads once and draws twice. See guiPath_doc for the
% format and guiPath_preset for how presets are found, loaded and saved.
%
% The LFP shown is the channel ED detection ran on, resolved here and frozen
% into a plain recipe, so the loader stays generic and a saved view built on
% this preset re-resolves it on the next session.
%
% INPUTS
% - ctx             <struct> var_ctx: basepath, basename, shared file cache.
%
% OUTPUTS
% - varMap          <struct> name -> recipe.
% - guiMap          <struct> .panels + .mode .win .save.
%
% SEE ALSO
% - guiPath_preset, var_recipe, guiPath_panel, guiPath_doc.
%
% HISTORY
% - 260719          split out of guiPath_presets (one file per preset).


%% ========================================================================
%  DATA (varMap)
%  ========================================================================
% the sleep context every modality shows on top, then the EDs + units
varMap = struct();
varMap.hyp    = var_recipe('matvar', 'file', 'sleep_states', 'var', 'ss', ...
    'path', 'bouts.times');
varMap.spec   = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', {'spec', 'spec_freq', 'spec_tstamps'});
varMap.emgRms = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', 'emg_rms', 'fs', 1);
varMap.emg    = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg');
varMap.ed     = var_recipe('matvar', 'file', 'ed');
varMap.raster = var_recipe('matvar', 'file', 'spikes', 'var', 'spikes', ...
    'path', 'times');

% the LFP this preset shows = the channel ED detection ran on (from ed.info):
% an 'lfp' source -> that raw channel; 'eeg' (the default) -> sSig.eeg
[edSrc, edCh] = edSigSource(ctx);
if strcmpi(edSrc, 'lfp') && ~isempty(edCh)
    varMap.lfp = var_recipe('bin', 'file', 'lfp', 'ch', edCh(:)', ...
        'average', true, 'outClass', 'native');
else
    varMap.lfp = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'eeg');
end


%% ========================================================================
%  VIEW (guiMap)
%  ========================================================================
% Top: state, spectrogram, EMG RMS, event ticks. Bottom: LFP, EMG, units, state
guiMap = struct('panels', struct(), 'mode', 'events', 'win', 1, 'save', 'ed');
guiMap.panels.hypT   = guiPath_panel('hypnogram', 'top', 'hyp', ...
    'label', 'State');
guiMap.panels.spec   = guiPath_panel('spec', 'top', 'spec');
guiMap.panels.emgRms = guiPath_panel('trace', 'top', 'emgRms', ...
    'label', 'EMG RMS', 'height', 0.7, 'ylim', 'full');
guiMap.panels.evt    = guiPath_panel('eventTicks', 'top', 'ed', ...
    'label', 'Events');
guiMap.panels.lfp    = guiPath_panel('trace', 'bottom', 'lfp', ...
    'label', 'LFP', 'height', 1.2);
guiMap.panels.emg    = guiPath_panel('trace', 'bottom', 'emg', ...
    'label', 'EMG', 'height', 0.8);
guiMap.panels.raster = guiPath_panel('raster', 'bottom', 'raster', ...
    'label', 'Units');
guiMap.panels.hypB   = guiPath_panel('hypnogram', 'bottom', 'hyp', ...
    'label', 'State');

end


% =========================================================================
%  SESSION RESOLVERS (read via the shared ctx cache)
% =========================================================================

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

% EOF
