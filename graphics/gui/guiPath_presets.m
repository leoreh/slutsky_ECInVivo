function [cfgData, cfgGui] = guiPath_presets(varargin)

% Returns a curation preset (panels and behaviour) for guiPath_curate.
%
% EXAMPLES
% - [cfgRipp, cfgGui] = guiPath_presets('ripp')
%   returns the typical panels used when visualizing ripples. The panels
%   come back slim (addresses only, no data); pass cfgData to guiPath_load
%   to fill them with data.
%
% - [cfgData, cfgGui] = guiPath_presets('ed', cfgRipp)
%   reuses loaded panels. It will take the previous cfgData and add to it
%   panels relevant to epileptiform discharges. Thus, switching presets
%   loads only what is new.
%
% - list = guiPath_presets()
%   returns the preset list as a struct array with fields {name, file}.
%   guiPath_curate uses it for the dropdown and file auto-detection.
%
% INPUTS
% - name            <char>(opt) which preset to build. Case-insensitive.
%                   File tokens also work.
%                       'EDs' or 'ed'.
%                       'Ripples' or 'ripp'.
%                       'States' or 'sleep_states'.
%                       'template' (empty cfgData).
% - cfgPrev         <struct>(opt) A cfgData from a previous call whose
%                   already-loaded panels are reused.
%
% OUTPUTS
% - cfgData         <struct> one field per panel. See guiPath_panel for a
%                   panel's fields.
% - cfgGui          <struct> the behaviour guiPath_curate needs, with fields:
%                       .name  display name.
%                       .file  session file token.
%                       .mode  'events' | 'states'.
%                       .win   default window width [s].
%                       .save  target token: ('ed' | 'ripp' | 'labelsMan'),
%                              '', or a save(x) handle.
% - list            <struct> preset list with fields {name, file}.
%
% SEE ALSO
% - guiPath_curate
% - guiPath_panel
% - guiPath_load
% - guiPath_src
% - guiPath_doc
%
% HISTORY
% - 260623          extracted from ed_gui / ripp_curate.
% - 260705          one preset by name -> [cfgData, cfgGui];
%                   reuse of loaded panels
%                   no bundled registry of data.
% - 260706          dropped the handle registry for an explicit switch.
% - 260706          renamed to guiPath_presets; GUI package flattened to gui_*.


%% ========================================================================
%  ARGUMENTS
%  ========================================================================

% Enumeration: with no name, return the preset list {name, file}. guiPath_curate
% uses it for the dropdown and to auto-detect which preset's file exists.
% Keep this list in sync with the switch below.
if nargin < 1
    cfgData = struct('name', {'EDs', 'Ripples', 'States'}, ...
                     'file', {'ed',  'ripp',    'sleep_states'});
    return
end

p = inputParser;
addOptional(p, 'name',    '', @(x) ischar(x) || isstring(x));
addOptional(p, 'cfgPrev', [], @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});

name    = char(p.Results.name);
cfgPrev = p.Results.cfgPrev;


%% ========================================================================
%  BUILD PRESET
%  ========================================================================
% Each case builds its panels (cfgData) and its gui behaviour (cfgGui).
switch lower(name)

    case {'eds', 'ed'}
        cfgData = cfg_ed();
        cfgGui  = struct('name', 'EDs', 'file', 'ed', ...
            'mode', 'events', 'win', 1, 'save', 'ed');

    case {'ripples', 'ripp'}
        cfgData = cfg_ripp();
        cfgGui  = struct('name', 'Ripples', 'file', 'ripp', ...
            'mode', 'events', 'win', 0.4, 'save', 'ripp');

    case {'states', 'sleep_states'}
        cfgData = cfg_states();
        cfgGui  = struct('name', 'States', 'file', 'sleep_states', ...
            'mode', 'states', 'win', 10, 'save', 'labelsMan');

    case {'template', 'custom', ''}
        cfgData = struct();
        cfgGui  = struct('name', 'Custom', 'file', '', ...
            'mode', 'events', 'win', 1, 'save', '');

    otherwise
        error('guiPath_presets:name', ...
            'unknown preset "%s" (EDs | Ripples | States | template)', name);

end


%% ========================================================================
%  REUSE LOADED PANELS
%  ========================================================================
% If a previous cfgData is given, carry its data into panels that match, so
% switching presets only loads what is new.
if ~isempty(cfgPrev)
    cfgData = reuse_loaded(cfgData, cfgPrev);
end

end


%% ========================================================================
%  REUSE HELPER
%  ========================================================================
function cfgData = reuse_loaded(cfgData, cfgPrev)
% Carry loaded data from cfgPrev into panels matching by field name and
% address, so a panel already loaded is not loaded again.

fldNames = fieldnames(cfgData);

for iFld = 1:numel(fldNames)

    fldName = fldNames{iFld};
    if ~isfield(cfgPrev, fldName), continue; end

    pnlNew = cfgData.(fldName);
    pnlOld = cfgPrev.(fldName);
    if ~isfield(pnlOld, 'data') || isempty(pnlOld.data), continue; end

    % Same panel type and source: the loaded data still applies.
    if strcmp(pnlNew.type, pnlOld.type) && isequal(pnlNew.src, pnlOld.src)
        pnlNew.data = pnlOld.data;
        if isfield(pnlOld, 'fs'),   pnlNew.fs   = pnlOld.fs;   end
        if isfield(pnlOld, 'ylim'), pnlNew.ylim = pnlOld.ylim; end
        cfgData.(fldName) = pnlNew;
    end

end
end


%% ========================================================================
%  PANEL BUILDERS
%  ========================================================================
% One field per panel. 'name' ties panels that share an input or a target.
function c = cfg_ed()
% EDs: events from <basename>.ed.mat, shown on the sleep context.

c = struct();

% top region: state, spectrogram, EMG RMS, event ticks
c.hypT = guiPath_panel('hypnogram', 'top', 'sleep_states:ss.bouts.times', ...
    'name', 'hypnogram', 'label', 'State');

c.spec = guiPath_panel('spec', 'top', 'fn:spec');

c.emgRms = guiPath_panel('trace', 'top', 'sleep_sig:emg_rms', 'fs', 1, ...
    'label', 'EMG RMS', 'height', 0.7, 'ylim', 'full');

c.evt = guiPath_panel('eventTicks', 'top', 'ed', 'name', 'eventTicks', ...
    'label', 'Events');

% bottom region: LFP (the channel ED was detected on, via ed.info), EMG,
% unit raster, state
c.lfp = guiPath_panel('trace', 'bottom', 'fn:edLfp', 'label', 'LFP', ...
    'height', 1.2);

c.emg = guiPath_panel('trace', 'bottom', 'sleep_sig:emg', 'label', 'EMG', ...
    'height', 0.8);

c.raster = guiPath_panel('raster', 'bottom', 'spikes:spikes.times', ...
    'label', 'Units');

c.hypB = guiPath_panel('hypnogram', 'bottom', 'sleep_states:ss.bouts.times', ...
    'name', 'hypnogram', 'label', 'State');

end


function c = cfg_ripp()
% Ripples: events from <basename>.ripp.mat; ripple-band LFP in the window.

c = struct();

% top region: state, spectrogram, EMG RMS, ripple ticks
c.hypT = guiPath_panel('hypnogram', 'top', 'sleep_states:ss.bouts.times', ...
    'name', 'hypnogram', 'label', 'State');

c.spec = guiPath_panel('spec', 'top', 'fn:spec');

c.emgRms = guiPath_panel('trace', 'top', 'sleep_sig:emg_rms', 'fs', 1, ...
    'label', 'EMG RMS', 'height', 0.7, 'ylim', 'full');

c.evt = guiPath_panel('eventTicks', 'top', 'ripp', 'name', 'eventTicks', ...
    'label', 'Ripples');

% bottom region: ripple LFP, filtered LFP, EMG, unit raster
c.rippleLfp = guiPath_panel('trace', 'bottom', 'fn:ripple.lfp', ...
    'label', 'LFP', 'height', 1.2);

c.rippleFilt = guiPath_panel('trace', 'bottom', 'fn:ripple.filt', ...
    'label', 'Filtered LFP', 'height', 1.0);

c.emg = guiPath_panel('trace', 'bottom', 'sleep_sig:emg', 'label', 'EMG', ...
    'height', 0.8);

c.raster = guiPath_panel('raster', 'bottom', 'spikes:spikes.times', ...
    'label', 'Units');

end


function c = cfg_states()
% States: the editable stateStrip is the target; top and bottom share input.

c = struct();

% top region: editable state strip, spectrogram, EMG RMS
c.stripT = guiPath_panel('stateStrip', 'top', '', 'name', 'states', ...
    'label', 'State');

c.spec = guiPath_panel('spec', 'top', 'fn:spec');

c.emgRms = guiPath_panel('trace', 'top', 'sleep_sig:emg_rms', 'fs', 1, ...
    'label', 'EMG RMS', 'height', 0.7, 'ylim', 'full');

% bottom region: LFP, EMG, editable state strip
c.eeg = guiPath_panel('trace', 'bottom', 'sleep_sig:eeg', 'label', 'LFP', ...
    'height', 1.2);

c.emg = guiPath_panel('trace', 'bottom', 'sleep_sig:emg', 'label', 'EMG', ...
    'height', 0.8);

c.stripB = guiPath_panel('stateStrip', 'bottom', '', 'name', 'states', ...
    'label', 'State');

end
