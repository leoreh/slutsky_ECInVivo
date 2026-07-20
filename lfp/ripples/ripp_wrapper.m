function ripp = ripp_wrapper(varargin)
% RIPP_WRAPPER Run the ripple pipeline end to end (detect -> curate -> analyze).
%
%   ripp = RIPP_WRAPPER(varargin)
%
%   SUMMARY:
%       A thin orchestrator over the three composable ripple stages, for the
%       batch (no-manual-curation) case:
%         1. ripp_detect  - signal -> events + per-event features (light).
%         2. ripp_curate  - apply the met.qa filter -> .accepted (headless here).
%         3. ripp_analyze - spikes / phase / maps on the ACCEPTED events (heavy).
%       Each stage is also a standalone function, so a session that needs manual
%       curation is run à la carte: ripp_wrapper(...,'flgCurate',true) detects,
%       saves, and opens the curation GUI, then STOPS - you curate, save in the
%       GUI, and call ripp_analyze(basepath) yourself once the mask is right.
%
%       Re-run model: if <basename>.ripp.mat exists and flgForce is false, the
%       stored (already detected + curated + analysed) result is loaded and
%       returned - the cheap path. flgForce re-detects and re-runs the chain.
%
%   INPUTS (Parameter/Value):
%       'basepath'   - (Char)   Session directory. {pwd}
%       'basename'   - (Char)   File stem. {folder name}
%       'met'        - (Struct) Detection + QA config. {ripp_methods('default')}
%       'rippCh'     - (Num)    Explicit 1-indexed channel(s); else resolved. {[]}
%       'win'        - (Vec)    Analysis window [start end] (s). {[0 Inf]}
%       'mapDur'     - (Vec)    PETH / map window [pre post] (s). {[-0.1 0.1]}
%       'flgPlot'    - (Log)    Spike-modulation summary figure (analyze)? {true}
%       'flgSave'    - (Log)    Save the output .mat files? {false}
%       'flgNS'      - (Log)    Write the NeuroScope event file? {false}
%       'flgCurate'  - (Log)    Detect + save, then open the curation GUI and
%                               stop (run ripp_analyze afterwards)? {false}
%       'flgDetectOnly' - (Log) Detect + save only, then stop - for the staged
%                               workflow (detect all mice, curate, then analyze
%                               all mice). {false}
%       'flgForce'   - (Log)    Re-detect even if .ripp.mat exists? {false}
%       'verbose'    - (Log)    Print progress? {true}
%
%   OUTPUT:
%       ripp - (Struct) The events struct. After the full chain it carries all
%           detected events with per-event metrics, the curated .accepted mask,
%           and the accepted-aligned .spks. Analysis filters .accepted.
%
%   FILES SAVED (when flgSave = true):
%       basename.ripp.mat        - all events + metrics + .accepted (+ .spks)
%       basename.rippMaps.mat    - per-event LFP maps, ALL detected events, row-
%                                  aligned to ripp; written at detect (single)
%       basename.rippStates.mat  - per-bout rate/density over accepted (curate)
%       basename.rippSpks.mat    - per-unit spike stats + PETH (analyze)
%       basename.rippSpkMaps.mat - 3D spike raster (analyze)
%       basename.rippSpkLfp.mat  - spike-LFP phase coupling (analyze)
%
%   DEPENDENCIES:
%       ripp_methods, ripp_detect, ripp_curate, ripp_analyze, evt_gate,
%       evt_files, evt_maps, evt2ns, backup_file, ripp_invalidate.
%
%   HISTORY:
%       Updated: 260719 (met-driven; detect core factored to ripp_detect).
%       Updated: 260719b (split into detect -> curate -> analyze stages; QA gate
%                moved to ripp_curate/evt_gate; heavy analysis runs post-curation
%                on the accepted set only).
%       Updated: 260720 (rippMaps written here at detect, over ALL events and
%                row-aligned to ripp, so the curation GUI reads it instead of
%                rebuilding the signal; readers subset it by .accepted).

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
addParameter(p, 'flgDetectOnly', false, @islogical);
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
flgDetectOnly = p.Results.flgDetectOnly;
flgForce  = p.Results.flgForce;
verbose   = p.Results.verbose;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end
if isempty(met), met = ripp_methods('default'); end

files = evt_files(basepath, basename, 'ripp');
if verbose, fprintf('[RIPP]: Session %s\n', basename); end

%% ========================================================================
%  STAGE 1: DETECT (fresh) OR LOAD (cheap re-run)
%  ========================================================================
aux   = [];
fresh = flgForce || ~isfile(files.evt);

if ~fresh
    if verbose, fprintf('[RIPP]: Loading existing %s.ripp.mat\n', basename); end
    S = load(files.evt, 'ripp');
    ripp = S.ripp;
else
    [ripp, aux] = ripp_detect(basepath, 'met', met, 'win', win, ...
        'rippCh', rippCh, 'verbose', verbose);
    ripp = finalizeRipp(ripp, win, basename);
    if flgSave || flgCurate || flgDetectOnly
        backup_file(files.evt);                 % preserve any prior curated mask
        save(files.evt, 'ripp', '-v7.3');       % detect output (accepted all-true)
        ripp_invalidate(basepath, basename, 'flgStates', true);   % prior products stale
        saveMaps(files.maps, aux, ripp, win, mapDur);   % after: maps are not stale
    end
end

% detect-only stops here (curate + analyze are run as separate steps/loops)
if flgDetectOnly
    if verbose, fprintf('[RIPP]: detect only (%s).\n', basename); end
    return;
end

%% ========================================================================
%  STAGE 2: CURATE (manual GUI, or the headless gate)
%  ========================================================================
if flgCurate
    if isfile(files.evt)
        if verbose, fprintf('[RIPP]: Opening curation GUI...\n'); end
        ripp_curate(basepath, 'basename', basename, 'qa', met.qa);
        if verbose
            fprintf(['[RIPP]: curate + save in the GUI, then run ' ...
                'ripp_analyze(''%s'').\n'], basepath);
        end
    elseif verbose
        fprintf('[RIPP]: No %s.ripp.mat; set flgSave/flgForce to curate.\n', ...
            basename);
    end
    return;                                     % analyze is a separate manual step
end

if fresh
    if flgSave
        % headless gate; analyze runs next and overwrites the products, so no
        % need to invalidate them here
        ripp_curate(basepath, 'basename', basename, 'qa', met.qa, ...
            'flgGui', false, 'flgInvalidate', false, 'verbose', verbose);
        S = load(files.evt, 'ripp');            % reload the curated mask
        ripp = S.ripp;
    else
        ripp.accepted = evt_gate(ripp, met.qa);    % in-memory preview
    end
end

%% ========================================================================
%  STAGE 3: ANALYZE (accepted events; needs the saved struct)
%  ========================================================================
if fresh && flgSave
    ripp = ripp_analyze(basepath, 'basename', basename, 'aux', aux, ...
        'mapDur', mapDur, 'flgSave', true, 'flgPlot', flgPlot, 'verbose', verbose);
end

%% ========================================================================
%  NEUROSCOPE
%  ========================================================================
if flgNS
    if verbose, fprintf('[RIPP]: Writing NeuroScope events...\n'); end
    accepted = true(size(ripp.times, 1), 1);
    if isfield(ripp, 'accepted'), accepted = ripp.accepted; end
    evt2ns(ripp.times, ripp.peakTime, 'basepath', basepath, ...
        'basename', basename, 'fileTag', 'rip', 'lbl', 'Ripple', ...
        'accepted', accepted);
end

if verbose, fprintf('[RIPP]: Done (%s).\n', basename); end

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function ripp = finalizeRipp(ripp, win, basename)
% Shift the window-relative detect output to absolute time and stamp provenance.
ripp.times     = ripp.times + win(1);
ripp.peakTime  = ripp.peakTime + win(1);
ripp.ctrlTimes = ripp.ctrlTimes + win(1);
ripp.info.basename = basename;
ripp.info.win      = win;
ripp.info.runtime  = datetime('now');
end     % finalizeRipp


function saveMaps(file, aux, ripp, win, mapDur)
% Per-event LFP maps for EVERY detected event, row-aligned to ripp.peakTime.
%
% Written at detection because this is the one stage that already holds the
% prepared signal - the maps themselves are a fraction of a second, while
% rebuilding them later costs a full signal load + prep. The curation GUI needs
% one row per event before any mask exists, so an accepted-only product could
% not serve it. Being all-events makes the file mask-INDEPENDENT: a curate save
% cannot stale it (readers subset by .accepted, as they do for ripp), which is
% why ripp_invalidate leaves it alone and only a re-detection rewrites it.
%
% Stored single: display and averaging data, where ~7 significant digits sit far
% below the measurement, and it halves a file that holds every event.
if isempty(aux) || ~isfield(aux, 'sig') || ~isfield(aux.sig, 'rippSig')
    return;
end
w0 = win(1);
if ~isfinite(w0), w0 = 0; end
rippMaps = evt_maps(aux.sig.rippSig, ripp.peakTime - w0, aux.fs, ...
    'mapDur', mapDur);

fn = fieldnames(rippMaps);
for iFld = 1 : numel(fn)
    if ~strcmp(fn{iFld}, 'tstamps') && isnumeric(rippMaps.(fn{iFld}))
        rippMaps.(fn{iFld}) = single(rippMaps.(fn{iFld}));
    end
end
save(file, 'rippMaps', '-v7.3');
end     % saveMaps
