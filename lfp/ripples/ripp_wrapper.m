function ripp = ripp_wrapper(varargin)
% RIPP_WRAPPER Detect, characterise, and curate sharp-wave ripples (SWR).
%
%   ripp = RIPP_WRAPPER(varargin)
%
%   SUMMARY:
%       Orchestrates the ripple pipeline for one session:
%       1. Setup: resolve the method (met) and file paths.
%       2. Detect-or-load: if <basename>.ripp.mat exists and flgForce is false,
%          load it (the cheap re-curate path); otherwise call ripp_detect (the
%          shared core: signal, detection, params, LFP maps, QA metrics, state,
%          and the .accepted mask), then run the spike analysis on the ACCEPTED
%          events, convert to absolute time, save, and plot.
%       3. Optionally export a NeuroScope event file (evt2ns) and launch the
%          curation GUI (guiPath, 'Ripples'). The export sits outside the
%          detect-or-load branch, so a re-run refreshes it against the stored
%          .accepted mask.
%
%       QA marks, it does not remove: <basename>.ripp.mat holds every detected
%       event with an .accepted flag (NREM/valid state, low EMG, above the MUA
%       gain gate). Downstream analysis filters .accepted; wake events remain for
%       inspection; and a gate threshold can be re-screened from the saved
%       .emg/.spkGain without re-detecting. The heavier products (maps, spikes,
%       phase, per-bout rates) are built on the accepted events.
%
%   INPUTS (Parameter/Value):
%       'basepath'   - (Char)   Session directory. {pwd}
%       'basename'   - (Char)   File stem. {folder name}
%       'met'        - (Struct) Detection + QA config. {ripp_methods('default')}
%       'rippCh'     - (Num)    Explicit 1-indexed channel; else resolved. {[]}
%       'win'        - (Vec)    Analysis window [start end] (s). {[0 Inf]}
%       'mapDur'     - (Vec)    PETH / map window [pre post] (s). {[-0.1 0.1]}
%       'flgPlot'    - (Log)    Generate the summary figure? {true}
%       'flgSave'    - (Log)    Save output .mat files? {false}
%       'flgNS'      - (Log)    Write the NeuroScope event file? {false}
%       'flgCurate'  - (Log)    Launch the curation GUI? {false}
%       'flgForce'   - (Log)    Re-detect even if .ripp.mat exists? {false}
%       'verbose'    - (Log)    Print progress? {true}
%
%   OUTPUT:
%       ripp         - (Struct) Events + per-event metrics. ALL detected events;
%           .times .peakTime .state .accepted .ctrlTimes .info, ripple params
%           (.amp .freq .freqEvent .freqPeak .peakProm .energy .dur .skew), QA
%           metrics (.emg .spkGain), and per-event population spike metrics
%           (.spks, aligned to the ACCEPTED events).
%
%   FILES SAVED (when flgSave = true):
%       basename.ripp.mat        - all events + per-event metrics + .accepted
%       basename.rippStates.mat  - per-bout rate/density over accepted events
%       basename.rippMaps.mat    - per-event LFP maps (accepted events)
%       basename.rippSpks.mat    - per-unit spike stats + PETH (light)
%       basename.rippSpkMaps.mat - 3D spike raster [unit x event x bin] (heavy)
%       basename.rippSpkLfp.mat  - spike-LFP phase coupling
%
%   FILES SAVED (when flgNS = true; independent of flgSave):
%       basename.rip.evt         - NeuroScope events. Tagged 'rip', not 'ripp':
%       NeuroScope loads an event file only when the id is exactly three
%       characters. See evt2ns.
%
%   DEPENDENCIES:
%       ripp_methods, ripp_detect, spklfp_phase, guiPath; and the shared event
%       layer (lfp/events): evt_files, evt_spks, evt_saveSpks, evt_plotSpks,
%       evt2ns.
%
%   HISTORY:
%       Updated: 260706 (shared spine via evt_* helpers; drop the steps modes).
%       Updated: 260715 (restore the NeuroScope export as flgNS).
%       Updated: 260719 (met-driven; detect core factored to ripp_detect; QA is
%                now an .accepted mask, not a removal).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'met', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'rippCh', [], @isnumeric);
addParameter(p, 'win', [0, Inf], @isnumeric);
addParameter(p, 'mapDur', [-0.1 0.1], @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgNS', false, @islogical);
addParameter(p, 'flgCurate', false, @islogical);
addParameter(p, 'flgForce', false, @islogical);
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});

basepath  = p.Results.basepath;
met       = p.Results.met;
rippCh    = p.Results.rippCh;
win       = p.Results.win;
mapDur    = p.Results.mapDur;
flgPlot   = p.Results.flgPlot;
flgSave   = p.Results.flgSave;
flgNS     = p.Results.flgNS;
flgCurate = p.Results.flgCurate;
flgForce  = p.Results.flgForce;
verbose   = p.Results.verbose;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end
if isempty(met), met = ripp_methods('default'); end

files = evt_files(basepath, basename, 'ripp');
files.phase  = fullfile(basepath, [basename, '.rippSpkLfp.mat']);
files.states = fullfile(basepath, [basename, '.rippStates.mat']);

if verbose, fprintf('[RIPP]: Session %s\n', basename); end

%% ========================================================================
%  DETECT OR LOAD
%  ========================================================================
if isfile(files.evt) && ~flgForce
    if verbose, fprintf('[RIPP]: Loading existing %s.ripp.mat\n', basename); end
    S = load(files.evt, 'ripp');
    ripp = S.ripp;

else
    % ---- Core: signal, detect, params, maps, QA metrics, state, accepted ----
    [ripp, aux] = ripp_detect(basepath, 'met', met, 'win', win, ...
        'rippCh', rippCh, 'mapDur', mapDur, 'verbose', verbose);
    rippMaps = aux.rippMaps;
    acc = ripp.accepted;

    % ---- Spiking on the accepted events (SU/MU modulation + PETH; phase) ----
    hasSpks = ~isempty(aux.spkTimes) && any(acc);
    if hasSpks
        if verbose, fprintf('[RIPP]: Analysing spikes...\n'); end
        rippSpks = evt_spks(aux.spkTimes, aux.muTimes, ripp.times(acc, :), ...
            ripp.ctrlTimes(acc, :), ripp.peakTime(acc), 'unitType', aux.uType, ...
            'mapDur', mapDur);
        ripp.spks = rippSpks.events;      % per-event metrics (accepted-aligned)
        rippSpks = rmfield(rippSpks, 'events');

        spkLfp = spklfp_phase(aux.sig.rippSig.filt, aux.spkTimes, aux.fs, ...
            'lfpTimes', ripp.times(acc, :), 'nPerms', 0);
    end

    % maps are analysis-facing: keep the accepted events only
    mapFlds = setdiff(fieldnames(rippMaps), {'tstamps'});
    for iFld = 1:numel(mapFlds)
        rippMaps.(mapFlds{iFld}) = rippMaps.(mapFlds{iFld})(acc, :);
    end

    % ---- Finalise: absolute time + provenance ---------------------------
    ripp.times     = ripp.times + win(1);
    ripp.peakTime  = ripp.peakTime + win(1);
    ripp.ctrlTimes = ripp.ctrlTimes + win(1);
    ripp.info.basename = basename;
    ripp.info.win      = win;
    ripp.info.runtime  = datetime('now');

    % ---- Save -----------------------------------------------------------
    if flgSave
        if verbose, fprintf('[RIPP]: Saving output files...\n'); end
        save(files.evt, 'ripp', '-v7.3');
        save(files.maps, 'rippMaps', '-v7.3');
        rippStates = aux.rippStates;
        save(files.states, 'rippStates', '-v7.3');
        if hasSpks
            evt_saveSpks(rippSpks, basepath, basename, 'ripp');
            save(files.phase, 'spkLfp', '-v7.3');
        end
    end

    % ---- Plot -----------------------------------------------------------
    if flgPlot && hasSpks
        if verbose, fprintf('[RIPP]: Generating summary plot...\n'); end
        evt_plotSpks(rippSpks, 'basepath', basepath, 'flgSaveFig', true, ...
            'name', 'ripp', 'lbl', 'Ripple');
    end
end

%% ========================================================================
%  NEUROSCOPE
%  ========================================================================
% Placed outside the branch so a re-run refreshes the export against the stored
% .accepted mask; rejected events stay in the file under their own palette entry.
if flgNS
    if verbose, fprintf('[RIPP]: Writing NeuroScope events...\n'); end
    accepted = true(size(ripp.times, 1), 1);
    if isfield(ripp, 'accepted'), accepted = ripp.accepted; end
    evt2ns(ripp.times, ripp.peakTime, 'basepath', basepath, ...
        'basename', basename, 'fileTag', 'rip', 'lbl', 'Ripple', ...
        'accepted', accepted);
end

%% ========================================================================
%  CURATE
%  ========================================================================
if flgCurate
    if isfile(files.evt)
        if verbose, fprintf('[RIPP]: Launching curation GUI...\n'); end
        guiPath(basepath, 'preset', 'Ripples', 'basename', basename);
    elseif verbose
        fprintf('[RIPP]: No %s.ripp.mat; set flgSave=true to curate.\n', ...
            basename);
    end
end

if verbose, fprintf('[RIPP]: Done (%s).\n', basename); end

end     % EOF
