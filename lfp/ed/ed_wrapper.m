function ed = ed_wrapper(varargin)
% ED_WRAPPER Run the ED pipeline end to end (detect -> curate).
%
%   ed = ED_WRAPPER(varargin)
%
%   SUMMARY:
%       A thin orchestrator over the two ED stages for the batch case:
%       ed_detect (candidates + features) then ed_curate (the .accepted gate,
%       headless here). Each stage is also a standalone function.
%
%       CURATING ONE SESSION is a two-step pass, and both steps matter:
%           ed_wrapper(..., 'flgCurate', true)   % bulk: set the thresholds
%           guiPath(basepath, 'preset', 'ed')    % per-event: accept / reject
%       The first cuts the thousands of candidates a permissive detector
%       proposes down to the set worth looking at; the second walks that set one
%       event at a time. Skipping the first means stepping through everything,
%       which is not a task anyone finishes.
%
%       Re-run model: if <basename>.ed.mat exists and flgForce is false, the
%       stored result is loaded and returned - the cheap path, which touches no
%       signals. flgForce re-detects, backing the previous file up first so a
%       curated mask is never lost silently.
%
%       There is no analyze stage. The ED question is how many discharges a
%       mouse has and how they distribute over vigilance states, which edStates
%       and ed_tbl answer; spike, PETH and phase products are deliberately not
%       built. Adding one later means adding a stage, not unpicking this.
%
%   INPUTS (Parameter/Value):
%       'basepath'  - (Char)   Session directory. {pwd}
%       'basename'  - (Char)   File stem. {folder name}
%       'met'       - (Struct) Detection + QA config. {ed_methods('default')}
%       'edCh'      - (Num)    Explicit channel for met.chMode 'ripp'. {[]}
%       'win'       - (Vec)    Analysis window [start end] (s). {[0 Inf]}
%       'flgSave'   - (Log)    Save the output .mat files? {false}
%       'flgCurate' - (Log)    Detect + save, then open the bulk curation GUI
%                              and stop? {false}
%       'flgForce'  - (Log)    Re-detect even if .ed.mat exists? {false}
%       'verbose'   - (Log)    Print progress? {true}
%
%   OUTPUT:
%       ed - (Struct) All candidates with their per-event metrics and the
%            curated .accepted mask. Analysis filters .accepted (see ed_tbl).
%
%   FILES SAVED (when flgSave = true):
%       basename.ed.mat       - all candidates + metrics + .accepted
%       basename.edMaps.mat   - per-event LFP maps, ALL candidates, row-aligned
%                               to ed; written at detect (single precision)
%       basename.edStates.mat - per-bout rate / density over accepted (curate)
%
%   DEPENDENCIES:
%       ed_methods, ed_detect, ed_curate, evt_gate, evt_files, evt_maps,
%       backup_file.
%
%   HISTORY:
%       Created: 260622
%       Updated: 260706 (shared spine via the evt_* helpers).
%       Updated: 260720 (rebuilt as a thin detect -> curate chain mirroring
%                ripp_wrapper. QA marks instead of removing; the detection
%                arguments moved into ed_methods; the spike analyses were
%                dropped as unused by the ED question.)

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'met', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'edCh', [], @isnumeric);
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgCurate', false, @islogical);
addParameter(p, 'flgForce', false, @islogical);
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});

basepath  = p.Results.basepath;
met       = p.Results.met;
edCh      = p.Results.edCh;
win       = p.Results.win;
flgSave   = p.Results.flgSave;
flgCurate = p.Results.flgCurate;
flgForce  = p.Results.flgForce;
verbose   = p.Results.verbose;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end
if isempty(met), met = ed_methods('default'); end

files = evt_files(basepath, basename, 'ed');
if verbose, fprintf('[ED]: Session %s\n', basename); end

%% ========================================================================
%  STAGE 1: DETECT (fresh) OR LOAD (cheap re-run)
%  ========================================================================
fresh = flgForce || ~isfile(files.evt);

if ~fresh
    if verbose, fprintf('[ED]: Loading existing %s.ed.mat\n', basename); end
    S  = load(files.evt, 'ed');
    ed = S.ed;
else
    [ed, aux] = ed_detect(basepath, 'basename', basename, 'met', met, ...
        'win', win, 'edCh', edCh, 'verbose', verbose);
    relPeak = ed.peakTime;              % maps are built against the same signal

    ed.times    = ed.times + win(1);
    ed.peakTime = ed.peakTime + win(1);
    ed = rmfield(ed, {'pos', 'bouts'});     % sample indices; now meaningless
    ed.info.basename = basename;
    ed.info.win      = win;
    ed.info.runtime  = datetime('now');

    if flgSave || flgCurate
        backup_file(files.evt);         % preserve any prior curated mask
        save(files.evt, 'ed', '-v7.3');
        saveMaps(files.maps, aux, relPeak);
    end
end

%% ========================================================================
%  STAGE 2: CURATE (bulk GUI, or the headless gate)
%  ========================================================================
if flgCurate
    if verbose, fprintf('[ED]: Opening bulk curation GUI...\n'); end
    ed_curate(basepath, 'basename', basename, 'qa', met.qa);
    if verbose
        fprintf(['[ED]: set the thresholds and Save, then step the events ' ...
            'with guiPath(''%s'', ''preset'', ''ed'').\n'], basepath);
    end
    return;                             % the per-event pass is manual
end

if fresh
    if flgSave
        ed_curate(basepath, 'basename', basename, 'qa', met.qa, ...
            'flgGui', false, 'verbose', verbose);
        S  = load(files.evt, 'ed');     % reload the curated mask
        ed = S.ed;
    else
        ed.accepted = evt_gate(ed, met.qa);     % in-memory preview
    end
end

if verbose
    fprintf('[ED]: Done (%s): %d detected, %d accepted.\n', ...
        basename, numel(ed.accepted), nnz(ed.accepted));
end

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function saveMaps(file, aux, relPeak)
% Per-event LFP maps for EVERY candidate, row-aligned to ed.
%
% Written at detection because this is the one stage that already holds the
% prepared signal. Being all-events makes the file mask-INDEPENDENT: a curate
% save cannot stale it (readers subset by .accepted, as they do for ed), and the
% curation GUI needs one row per candidate before any mask exists. Stored
% single: display and averaging data, where ~7 significant digits sit far below
% the measurement, and it halves a file that holds every event.
MAPDUR = [-0.1 0.1];

edMaps = evt_maps(struct('lfp', aux.edSig.lfp, 'filt', aux.edSig.filt), ...
    relPeak, aux.fs, 'mapDur', MAPDUR);

fn = fieldnames(edMaps);
for iFld = 1 : numel(fn)
    if ~strcmp(fn{iFld}, 'tstamps')
        edMaps.(fn{iFld}) = single(edMaps.(fn{iFld}));
    end
end
save(file, 'edMaps', '-v7.3');

end     % saveMaps
