function [varMap, guiMap] = preset_ripp(ctx)

% guiPath preset: ripples from <basename>.ripp.mat, on the sleep context.
%
% varMap = what to load (var_recipe); guiMap = how to show it (guiPath_panel).
% guiPath joins them by name: a panel's .var picks a varMap entry, so a top +
% bottom pair over one var loads once and draws twice. See guiPath_doc for the
% format and guiPath_preset for how presets are found, loaded and saved.
%
% Session-aware values (the detection channel, its scaling, the passband) are
% resolved here and frozen into plain recipes, so the loader stays generic and a
% saved view built on this preset re-resolves them on the next session.
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
% the sleep context every modality shows on top, then the ripples + units
varMap = struct();
varMap.hyp    = var_recipe('matvar', 'file', 'sleep_states', 'var', 'ss', ...
    'path', 'bouts.times');
varMap.spec   = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', {'spec', 'spec_freq', 'spec_tstamps'});
varMap.emgRms = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', 'emg_rms', 'fs', 1);
varMap.emg    = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg');
varMap.ripp   = var_recipe('matvar', 'file', 'ripp');
varMap.raster = var_recipe('matvar', 'file', 'spikes', 'var', 'spikes', ...
    'path', 'times');

% the shank around the ripple channel (raw, stacked) + the filtered detection
% channel. Both follow the channel ripple detection ran on (ripp_pickCh); the
% filtered trace is ripp_sigPrep on the bit2uv-scaled average.
session = getSession(ctx);
if round(session.extracellular.sr) == 24414, b2u = 1; else, b2u = 0.195; end
pb = rippPassband(ctx);
varMap.rippStack = var_recipe('bin', 'file', 'lfp', ...
    'ch', ripp_pickCh(ctx.basepath, 'basename', ctx.basename), ...
    'average', false, 'outClass', 'native');
varMap.rippFilt  = var_recipe('bin', 'file', 'lfp', ...
    'ch', ripp_pickCh(ctx.basepath, 'basename', ctx.basename, 'session', session), ...
    'average', true, ...
    'bit2uv', b2u, 'transform', {'rippPrep', {pb}}, 'path', 'filt');


%% ========================================================================
%  VIEW (guiMap)
%  ========================================================================
% Top: state, spectrogram, EMG RMS, ripple ticks. Bottom: the shank, the
% filtered channel, EMG, units.
guiMap = struct('panels', struct(), 'mode', 'events', 'win', 0.4, ...
    'save', 'ripp');
guiMap.panels.hypT      = guiPath_panel('hypnogram', 'top', 'hyp', ...
    'label', 'State');
guiMap.panels.spec      = guiPath_panel('spec', 'top', 'spec');
guiMap.panels.emgRms    = guiPath_panel('trace', 'top', 'emgRms', ...
    'label', 'EMG RMS', 'height', 0.7, 'ylim', 'full');
guiMap.panels.evt       = guiPath_panel('eventTicks', 'top', 'ripp', ...
    'label', 'Ripples');
guiMap.panels.rippStack = guiPath_panel('traces', 'bottom', 'rippStack', ...
    'label', 'LFP', 'height', 2.0);
guiMap.panels.rippFilt  = guiPath_panel('trace', 'bottom', 'rippFilt', ...
    'label', 'Filtered LFP', 'height', 1.0);
guiMap.panels.emg       = guiPath_panel('trace', 'bottom', 'emg', ...
    'label', 'EMG', 'height', 0.8, 'ylim', 1);
guiMap.panels.raster    = guiPath_panel('raster', 'bottom', 'raster', ...
    'label', 'Units');

end


% =========================================================================
%  SESSION RESOLVERS (read via the shared ctx cache)
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

% EOF
