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
% - guiMap          <struct> .panels + .mode .win.
%
% SEE ALSO
% - guiPath_preset, var_recipe, guiPath_panel, guiPath_doc, stateSet.
%
% HISTORY
% - 260719          split out of guiPath_presets (one file per preset).
% - 260720          the state context is the shared, curatable stateSet (was a
%                   read-only hypnogram over ss.bouts.times).
% - 260720b         follows the rebuilt ED pipeline: the Bottom now carries the
%                   detection band-pass alongside the raw trace, because that
%                   filtered trace is what ed_detect thresholds - a candidate
%                   that looks unconvincing raw is judged on it.
% - 260721          ED detection reads one auto-picked .lfp channel, so the
%                   chMode branch and the sleep_sig eeg fallback are gone.


%% ========================================================================
%  DATA (varMap)
%  ========================================================================
% the sleep context every modality shows on top, then the EDs + units
varMap = struct();
% the state strip is omitted on a session with no sleep scoring; guiPath then
% drops the panels naming it, so the preset still opens on its signals
sSet = stateSet(ctx);
if ~isempty(sSet)
    varMap.states = var_recipe('value', 'data', sSet);
end
varMap.spec   = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', {'spec', 'spec_freq', 'spec_tstamps'});
varMap.emgRms = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', 'emg_rms', 'fs', 1);
varMap.emg    = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg');
varMap.ed     = var_recipe('matvar', 'file', 'ed');
varMap.raster = var_recipe('matvar', 'file', 'spikes', 'var', 'spikes', ...
    'path', 'times');

% the LFP this preset shows = the channel ED detection ran on (ed.info.edCh).
% The filtered twin is the same channel through the detection band-pass, which
% is the trace ed_detect actually thresholds.
[edCh, edBand] = edSigSource(ctx);
varMap.lfp = var_recipe('bin', 'file', 'lfp', 'ch', edCh(:)', ...
    'average', true, 'outClass', 'native');
varMap.edFilt = var_recipe('bin', 'file', 'lfp', 'ch', edCh(:)', ...
    'average', true, 'transform', {'bandpass', {edBand}});


%% ========================================================================
%  VIEW (guiMap)
%  ========================================================================
guiMap = struct('panels', struct(), 'base', 'ed', 'mode', 'events', 'win', 0.3);
guiMap.panels.states = guiPath_panel('stateStrip', 'top', 'states');
guiMap.panels.spec = guiPath_panel('spec', 'top', 'spec');
guiMap.panels.emgRms = guiPath_panel('trace', 'top', 'emgRms', ...
    'height', 0.7, 'label', 'EMG RMS', 'ylim', 'full');
guiMap.panels.ed = guiPath_panel('eventTicks', 'top', 'ed');
guiMap.panels.lfp = guiPath_panel('trace', 'bottom', 'lfp', ...
    'height', 1.2, 'label', 'LFP', 'yAdjust', 0.512);
guiMap.panels.edFilt = guiPath_panel('trace', 'bottom', 'edFilt', ...
    'label', 'Filtered LFP');
guiMap.panels.emg = guiPath_panel('trace', 'bottom', 'emg', ...
    'height', 0.8, 'label', 'EMG');
guiMap.panels.raster = guiPath_panel('raster', 'bottom', 'raster');
guiMap.panels.states_2 = guiPath_panel('stateStrip', 'bottom', 'states');
guiMap.panels.ed_2 = guiPath_panel('eventTicks', 'bottom', 'ed');

end


% =========================================================================
%  SESSION RESOLVERS (read via the shared ctx cache)
% =========================================================================

function [ch, band] = edSigSource(ctx)
% ed.info edCh / passband - the channel ED detection ran on and the band it
% thresholded. Defaults match ed_methods, so the preset still opens on a
% session whose ed.mat predates these fields.
ch = 1; band = [60 150];
try
    ed = var_fetch(var_recipe('matvar', 'file', 'ed'), ctx);
    if isfield(ed, 'info')
        if isfield(ed.info, 'edCh') && ~isempty(ed.info.edCh)
            ch = ed.info.edCh;
        end
        if isfield(ed.info, 'passband') && ~isempty(ed.info.passband)
            band = ed.info.passband;
        end
    end
catch
end
end

% EOF

