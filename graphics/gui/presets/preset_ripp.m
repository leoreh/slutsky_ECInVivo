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
% - 260720          rippFilt band-passes directly instead of running the full
%                   ripp_sigPrep and keeping only .filt.
% - 260720b         the overview EMG panel shows the gate's own score
%                   (evt_emgScore on a regular grid) instead of AccuSleep's
%                   emg_rms, on fixed y-limits, so a curation threshold is
%                   readable off the trace. Other presets keep emg_rms.


%% ========================================================================
%  DATA (varMap)
%  ========================================================================
% the sleep context every modality shows on top, then the ripples + units
varMap = struct();
varMap.hyp    = var_recipe('matvar', 'file', 'sleep_states', 'var', 'ss', ...
    'path', 'bouts.times');
varMap.spec   = var_recipe('matfield', 'file', 'sleep_sig', ...
    'field', {'spec', 'spec_freq', 'spec_tstamps'});
% The overview EMG is the ripple gate's own metric, not AccuSleep's log-RMS: the
% trace's y-value IS ripp.emg, so a threshold set in ripp_curate can be read off
% this panel. AccuSleep's emg_rms is a 1 s statistic - one bin spans ~20 ripples,
% so it can never say whether a given event is contaminated. The other presets
% (sleep, ed) keep emg_rms, which is the right signal for scoring states.
varMap.emgScore = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg', ...
    'transform', {'emgScore', {nremBouts(ctx), rippDur(ctx)}});
varMap.emg    = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg');
varMap.ripp   = var_recipe('matvar', 'file', 'ripp');
varMap.raster = var_recipe('matvar', 'file', 'spikes', 'var', 'spikes', ...
    'path', 'times');

% the shank around the ripple channel (raw, stacked) + the filtered detection
% channel. Both follow the channel ripple detection ran on (ripp_pickCh); the
% filtered trace is the detection band-pass of the bit2uv-scaled average - the
% same filter ripp_sigPrep applies for its .filt, without the hilbert /
% detection-signal / z-score work whose output no panel draws.
session = getSession(ctx);
if round(session.extracellular.sr) == 24414, b2u = 1; else, b2u = 0.195; end
pb = rippPassband(ctx);
varMap.rippStack = var_recipe('bin', 'file', 'lfp', ...
    'ch', ripp_pickCh(ctx.basepath, 'basename', ctx.basename), ...
    'average', false, 'outClass', 'native');
varMap.rippFilt  = var_recipe('bin', 'file', 'lfp', ...
    'ch', ripp_pickCh(ctx.basepath, 'basename', ctx.basename, 'session', session), ...
    'average', true, ...
    'bit2uv', b2u, 'transform', {'bandpass', {pb}});


%% ========================================================================
%  VIEW (guiMap)
%  ========================================================================
guiMap = struct('panels', struct(), ...
    'base', 'ripp', 'mode', 'events', 'win', 1, 'save', 'ripp');
guiMap.panels.hyp = guiPath_panel('hypnogram', 'top', 'hyp');
guiMap.panels.spec = guiPath_panel('spec', 'top', 'spec');
% absolute y-limits, not autoscale: a fixed scale is the whole point here, so
% the same height means the same score in every mouse. 0 is the resting NREM
% level and the WAKE/NREM boundary sits near 1. Shift+scroll (or shift +/-/0)
% over the panel widens or tightens this about its centre for a closer look.
guiMap.panels.emgScore = guiPath_panel('trace', 'top', 'emgScore', ...
    'height', 0.7, 'label', 'EMG score', 'ylim', [-2, 6]);
guiMap.panels.ripp = guiPath_panel('eventTicks', 'top', 'ripp', ...
    'label', 'Ripples');
guiMap.panels.rippStack = guiPath_panel('traces', 'bottom', 'rippStack', ...
    'height', 2, 'yAdjust', 1.953125);
guiMap.panels.rippFilt = guiPath_panel('trace', 'bottom', 'rippFilt', ...
    'label', 'Filtered LFP');
guiMap.panels.emg = guiPath_panel('trace', 'bottom', 'emg', ...
    'height', 0.8, 'label', 'EMG', 'ylim', 1);
guiMap.panels.raster = guiPath_panel('raster', 'bottom', 'raster');
guiMap.panels.ripp_2 = guiPath_panel('eventTicks', 'bottom', 'ripp', ...
    'label', 'Ripples');

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


function bt = nremBouts(ctx)
% the NREM baseline the detector used, taken through evt_boutTimes so the state
% convention is not restated here. Empty degrades gracefully - evt_emgScore then
% standardizes against the whole recording.
bt = [];
try
    ss = var_fetch(var_recipe('matvar', 'file', 'sleep_states', 'var', 'ss'), ctx);
    [~, ~, bt] = evt_boutTimes(struct('ss', ss), [0, Inf], Inf);
catch
end
end


function w = rippDur(ctx)
% median ripple duration: the window the score is evaluated over, so the trace
% sits on the ripples' own timescale. 50 ms before anything is detected.
w = 0.05;
try
    t = var_fetch(var_recipe('matvar', 'file', 'ripp', 'var', 'ripp', ...
        'path', 'times'), ctx);
    if ~isempty(t), w = median(diff(t, 1, 2)); end
catch
end
if ~isfinite(w) || w <= 0, w = 0.05; end
end

% EOF

