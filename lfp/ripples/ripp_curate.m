function [ripp, hFig] = ripp_curate(basepath, varargin)
% RIPP_CURATE Post-detection QA gate for ripples: headless or interactive (stage 2).
%
%   [ripp, hFig] = RIPP_CURATE(basepath, varargin)
%
%   SUMMARY:
%       The curation stage of the ripple pipeline (detect -> curate -> analyze).
%       It loads the saved <basename>.ripp.mat and sets the per-event .accepted
%       mask from a QA filter - which vigilance states to keep and per-metric
%       [lo hi] ranges (EMG, MUA gain, ...). The gate itself is evt_gate; this
%       function is the two ways to drive it:
%
%       - Headless (flgGui = false): apply the given qa spec, save .accepted (and
%         the spec in ripp.info.qa) back to ripp.mat, and rebuild the per-bout
%         rate/density table. This is the automatic gate - it replaces the old
%         evt_qa call, so a batch run needs no human.
%
%       - Interactive (flgGui = true, default): a GUI seeded from the same qa
%         spec. State checkboxes and metric thresholds recompute the kept/removed
%         split live - counts per state and the kept-vs-removed mean waveform
%         update as you move them - so a mouse can be curated by its own
%         judgement (states differ in scoring quality between mice). Save writes
%         the same way. Nothing is destroyed: rejected events keep their row with
%         .accepted = false, and Reset restores the default spec.
%
%       Saving only touches .accepted, ripp.info.qa, and rippStates (all cheap
%       and derived from the mask), and backs the file up first. The heavy
%       products (spikes, phase, averaged maps) are the analyze stage, run AFTER
%       curation on the accepted set - so they are never stale against the mask.
%
%   INPUTS:
%       basepath - <char> session directory (must hold <basename>.ripp.mat).
%       varargin - Parameter/Value:
%           'basename' - <char>   file stem. {folder name}
%           'qa'       - <struct> filter spec (see ripp_methods '.qa'); the
%                                 headless gate and the GUI's initial state.
%                                 {ripp_methods('default').qa}
%           'flgGui'   - <log>    open the GUI (true) or apply headless (false).
%                                 {true}
%           'flgInvalidate' - <log> when the mask changes, delete the stale
%                                 accepted-aligned analyze products (rippMaps /
%                                 rippSpks / rippSpkMaps / rippSpkLfp) so they
%                                 cannot be read stale before ripp_analyze reruns.
%                                 The batch turns this off (analyze overwrites
%                                 them next). {true}
%           'Visible'  - <char>   'on' | 'off' for headless GUI tests. {'on'}
%           'verbose'  - <log>    print progress? {true}
%
%   OUTPUTS:
%       ripp - <struct> the loaded events (headless: with .accepted updated).
%       hFig - <handle> the GUI figure ([] when headless).
%
%   DEPENDENCIES:
%       evt_files, evt_gate, evt_detrend, ripp_methods, backup_file,
%       ripp_invalidate, evt_states, evt_boutTimes, basepaths2vars;
%       GUI: gui_layout, gui_labeledControl, gui_filterPanel,
%       gui_selectedCats, guiTbl_xy, ripp_sigLoad, ripp_sigPrep, evt_maps,
%       as_loadConfig.
%
%   HISTORY:
%       260719b the curation stage; absorbs evt_qa's ripple role (headless gate)
%               and evolves ripp_gateGui into a disk-based bulk curator.
%       260720  waveform view tiled by state (one panel per state, kept vs
%               removed overlaid), so a threshold's effect is read per state.
%               The per-state count lines are gone - each tile's legend carries
%               them - and the state filter is sized to its list, not scrolled.
%       260720b the waveform maps are READ from the detect-stage rippMaps
%               (all events) instead of rebuilt from the signal; a session
%               without a matching file still falls back to rebuilding.

%% ========================================================================
%  ARGUMENTS + LOAD
%  ========================================================================

p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'qa', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'flgGui', true, @islogical);
addParameter(p, 'flgInvalidate', true, @islogical);
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(char(x), {'on', 'off'})));
addParameter(p, 'verbose', true, @islogical);
parse(p, basepath, varargin{:});
qa       = p.Results.qa;
flgGui   = p.Results.flgGui;
flgInval = p.Results.flgInvalidate;
vis      = char(p.Results.Visible);
verbose  = p.Results.verbose;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end
if isempty(qa), met = ripp_methods('default'); qa = met.qa; end

files = evt_files(basepath, basename, 'ripp');
if ~isfile(files.evt)
    error('ripp_curate:noFile', ...
        '%s not found; run detection first (ripp_wrapper).', files.evt);
end
S = load(files.evt, 'ripp');
ripp = S.ripp;
hFig = [];

%% ========================================================================
%  HEADLESS: apply the gate, save, rebuild states
%  ========================================================================

if ~flgGui
    ripp.accepted = evt_gate(ripp, qa);
    changed = saveCurated(files.evt, ripp.accepted, qa);
    buildStates(basepath, ripp, ripp.accepted);
    if changed && flgInval
        nDel = ripp_invalidate(basepath, basename);
        if verbose && nDel > 0
            fprintf(['[RIPP_CURATE] %s : removed %d stale analyze product(s); ' ...
                'rerun ripp_analyze\n'], basename, nDel);
        end
    end
    if verbose
        fprintf('[RIPP_CURATE] %s : %d / %d accepted (headless)\n', ...
            basename, nnz(ripp.accepted), numel(ripp.accepted));
    end
    return;
end

%% ========================================================================
%  GUI: state filter + metric thresholds -> live kept/removed
%  ========================================================================

% all-events LFP maps for the waveform view (best-effort; disabled if no signal)
st = struct();
st.maps = [];
st.xt = [];
try
    [st.maps, st.xt] = loadMaps(basepath, basename, ripp);
catch ME
    warning('ripp_curate:maps', 'waveform view disabled (%s)', ME.message);
end

hFig = uifigure('Name', ['Ripple curation: ' basename], ...
    'Position', [80 80 1500 850], 'Visible', vis);
[~, gPlot, gCtrl, gActions] = gui_layout(hFig, 'CtrlWidth', 250);

% present states, plus an explicit "(unscored)" entry for events whose peak
% falls in no scored bout (state <undefined>) - so they are visible and can be
% kept or dropped, rather than silently excluded by the categorical filter
realCats  = categories(removecats(ripp.state));
stateCats = realCats;
if any(isundefined(ripp.state))
    stateCats = [realCats(:); {'(unscored)'}];
end
defNames  = specStateNames(qa);           % default-checked (real states only)
initState = ismember(stateCats, defNames);
dflt      = specDefaults(qa);

% sized to the state list so it never scrolls (gui_filterPanel lays out 22 px
% per checkbox + 2 px spacing, inside 4 px of padding)
pnlState = gui_labeledControl(gCtrl, 'panel', 'Keep states:', ...
    'RowHeight', 24 * numel(stateCats) + 4);
st.chkState = gui_filterPanel(pnlState, stateCats, @(~,~) refresh(hFig), ...
    'InitVal', initState);

st.edGain = gui_labeledControl(gCtrl, 'editnum', 'MUA gain >=', ...
    'Value', dflt.gain, 'ValueChangedFcn', @(~,~) refresh(hFig));
st.edEmg = gui_labeledControl(gCtrl, 'editnum', 'EMG <=', ...
    'Value', dflt.emg, 'ValueChangedFcn', @(~,~) refresh(hFig));
st.edProm = gui_labeledControl(gCtrl, 'editnum', 'prominence >=', ...
    'Value', dflt.prom, 'ValueChangedFcn', @(~,~) refresh(hFig));

st.lblCount = gui_labeledControl(gCtrl, 'label', '');

gui_labeledControl(gActions, 'button', '', 'Text', 'Reset to default', ...
    'ButtonPushedFcn', @(~,~) onReset(hFig));
gui_labeledControl(gActions, 'button', '', 'Text', 'Save', ...
    'ButtonPushedFcn', @(~,~) doSave(hFig));

st.hPanel = uipanel(gPlot, 'BorderType', 'none');

st.ripp      = ripp;
st.qa0       = qa;
st.files     = files;
st.basepath  = basepath;
st.basename  = basename;
st.realStates = realCats;
st.flgInval  = flgInval;
st.accepted  = ripp.accepted;
st.stateCol  = plotState(ripp.state);   % tiling variable; fixed for the session
hFig.UserData = st;

refresh(hFig);

end     % EOF


% =========================================================================
%  GUI CALLBACKS
% =========================================================================
function refresh(hFig)
% Recompute accepted from the current controls; redraw counts + waveform.
st = hFig.UserData;
qa = buildSpec(st);
st.accepted = evt_gate(st.ripp, qa);
hFig.UserData = st;

% overall count; the per-state split is read off each tile's legend
st.lblCount.Text = sprintf('kept %d / %d', nnz(st.accepted), numel(st.accepted));

% kept-vs-removed waveform, one tile per state. Only the kept/removed flag
% changes as a threshold moves, so the widget is updated in place - rebuilding
% it would reset the Y / tile / group selections on every keystroke.
if ~isempty(st.maps)
    tbl = curateTbl(st);
    ud = st.hPanel.UserData;
    if isstruct(ud) && isfield(ud, 'setDataFcn')
        ud.setDataFcn(tbl);
    else
        guiTbl_xy(st.xt, tbl, 'Parent', st.hPanel, 'yVar', 'lfp', ...
            'tileVar', 'state', 'grpVar', 'status', 'xLbl', 'time (ms)');
    end
end

end     % refresh


function tbl = curateTbl(st)
% One row per detected event: its LFP waveform, whether the current filter keeps
% it (the group), and the state it falls in (the tile).
status = repmat("removed", numel(st.accepted), 1);
status(st.accepted) = "kept";
status = categorical(status, {'removed', 'kept'});
tbl = table(st.maps.lfp, status, st.stateCol, ...
    'VariableNames', {'lfp', 'status', 'state'});

end     % curateTbl


function s = plotState(state)
% ripp.state as a tiling variable: <undefined> is promoted to its own 'unscored'
% level. A categorical comparison never matches <undefined>, so without this the
% unscored events would get no tile and vanish from the view without a word.
s = removecats(state(:));
if any(isundefined(s))
    s = addcats(s, {'unscored'});
    s(isundefined(s)) = 'unscored';
end

end     % plotState


function onReset(hFig)
% Restore the controls (states + thresholds) to the default qa spec.
st = hFig.UserData;
dflt = specDefaults(st.qa0);
st.edGain.Value = dflt.gain;
st.edEmg.Value  = dflt.emg;
st.edProm.Value = dflt.prom;
defNames = specStateNames(st.qa0);
for iChk = 1:numel(st.chkState)
    st.chkState(iChk).Value = ismember(st.chkState(iChk).Text, defNames);
end
refresh(hFig);

end     % onReset


function doSave(hFig)
% Persist accepted + the spec, and rebuild the per-bout rate/density table.
st = hFig.UserData;
qa = buildSpec(st);
changed = saveCurated(st.files.evt, st.accepted, qa);
buildStates(st.basepath, st.ripp, st.accepted);
msg = sprintf('Saved: %d / %d accepted (+ rippStates)', ...
    nnz(st.accepted), numel(st.accepted));
if changed && st.flgInval
    nDel = ripp_invalidate(st.basepath, st.basename);
    if nDel > 0
        msg = sprintf('%s. Removed %d stale product(s) - rerun ripp_analyze.', ...
            msg, nDel);
    end
end
gui_notify(hFig, msg, 'success');

end     % doSave


% =========================================================================
%  SPEC <-> CONTROLS
% =========================================================================
function qa = buildSpec(st)
% Read the current controls into a qa filter spec. The "(unscored)" pseudo-state
% maps to qa.unscored (keep <undefined> events); the rest are real state labels.
% Guard the all-unchecked case: with real states present but none checked, keep
% NO real state - a sentinel that matches no event - rather than evt_gate's
% []="any state" escape hatch that would silently keep every state.
sel = gui_selectedCats(st.chkState);
qa.unscored = ismember('(unscored)', sel);
realSel = setdiff(sel, {'(unscored)'}, 'stable');
if isempty(realSel) && ~isempty(st.realStates)
    realSel = {'<none>'};
end
qa.states = realSel;
qa.ranges = struct('spkGain', [st.edGain.Value, Inf], ...
    'emg', [-Inf, st.edEmg.Value], ...
    'peakProm', [st.edProm.Value, Inf]);

end     % buildSpec


function d = specDefaults(qa)
% Pull the control default values (gain/emg/prom) out of a qa spec.
d.gain = -Inf;
d.emg  = Inf;
d.prom = -Inf;
if isfield(qa, 'ranges')
    if isfield(qa.ranges, 'spkGain'),  d.gain = qa.ranges.spkGain(1); end
    if isfield(qa.ranges, 'emg'),      d.emg  = qa.ranges.emg(2);     end
    if isfield(qa.ranges, 'peakProm'), d.prom = qa.ranges.peakProm(1); end
end

end     % specDefaults


function names = specStateNames(qa)
% Default-checked state labels from a qa spec (indices resolved via config).
names = {};
if ~isfield(qa, 'states') || isempty(qa.states)
    return;
end
if isnumeric(qa.states)
    cfg = as_loadConfig([]);
    idx = qa.states(qa.states >= 1 & qa.states <= numel(cfg.names));
    names = cfg.names(idx);
else
    names = cellstr(qa.states);
end

end     % specStateNames


% =========================================================================
%  PERSISTENCE
% =========================================================================
function changed = saveCurated(file, accepted, qa)
% Back up, then overwrite .accepted + .info.qa in the saved struct, and strip a
% now-stale accepted-aligned .spks (ripp_analyze rebuilds it). Returns whether
% the mask actually changed vs what is on disk.
backup_file(file);
S = load(file);
accepted = logical(accepted(:));
changed = ~isfield(S.ripp, 'accepted') || ...
    ~isequal(logical(S.ripp.accepted(:)), accepted);
S.ripp.accepted = accepted;
if isfield(S.ripp, 'spks')
    S.ripp = rmfield(S.ripp, 'spks');   % accepted-aligned; rebuilt by ripp_analyze
end
if ~isfield(S.ripp, 'info') || ~isstruct(S.ripp.info)
    S.ripp.info = struct();
end
S.ripp.info.qa = qa;
save(file, '-struct', 'S', '-v7.3');

end     % saveCurated


function buildStates(basepath, ripp, accepted)
% Rebuild + save the per-bout rate/density table for the current mask (cheap;
% no signal or spikes). Skips silently when sleep states are unavailable.
win = [0 Inf];
if isfield(ripp, 'info') && isfield(ripp.info, 'win'), win = ripp.info.win; end

v = basepaths2vars('basepaths', {basepath}, 'vars', {'session', 'sleep_states'});
fs = v.session.extracellular.srLfp;
if isinf(win(2)), win(2) = v.session.extracellular.nSamples / fs; end
sigDur = win(2) - win(1);
w0 = win(1);
if ~isfinite(w0), w0 = 0; end

boutTimes = evt_boutTimes(v, win, sigDur);
if isempty(boutTimes)
    return;                     % no scoring -> no rate/density table
end
evt_states(ripp.times - w0, ripp.peakTime - w0, boutTimes, ...
    'accepted', logical(accepted(:)), 'basepath', basepath, ...
    'flgSave', true, 'flgPlot', false, 'name', 'ripp', 'lbl', 'Ripple');

end     % buildStates


% =========================================================================
%  WAVEFORM MAPS (all events; loaded once for the GUI)
% =========================================================================
function [maps, xt] = loadMaps(basepath, basename, ripp)
% The per-event LFP maps behind the waveform view, cropped to DISPDUR.
%
% Detection writes rippMaps over ALL detected events (ripp_wrapper), so this is
% normally a file read - a second or so instead of the full signal load + prep.
% A session detected before that convention, or one whose file no longer matches
% the event list, falls back to rebuilding from the signal.
%
% The LFP is DETRENDED per event before the crop, and the order matters.
% evt_detrend fits its baseline on the flanks of whatever window it is given;
% on the SAVED map those flanks are far enough out to be background, while on
% the ±60 ms display crop they would still be inside the sharp wave and the
% detrend would eat a slice of it. Every event rides on its own drift, so
% without this the kept-vs-removed averages differ partly by whatever the
% drifts happened to do.
dispDur = [-0.06 0.06];
nEv = numel(ripp.peakTime);

files = evt_files(basepath, basename, 'ripp');
if isfile(files.maps)
    S = load(files.maps, 'rippMaps');
    if isfield(S, 'rippMaps') && isfield(S.rippMaps, 'lfp') && ...
            size(S.rippMaps.lfp, 1) == nEv
        S.rippMaps.lfp = evt_detrend(double(S.rippMaps.lfp), ...
            S.rippMaps.tstamps);
        [maps, xt] = cropMaps(S.rippMaps, dispDur);
        return;
    end
end

win = [0 Inf];
if isfield(ripp, 'info') && isfield(ripp.info, 'win'), win = ripp.info.win; end

v = basepaths2vars('basepaths', {basepath}, ...
    'vars', {'session', 'sleep_states'});
fs = v.session.extracellular.srLfp;
if isinf(win(2)), win(2) = v.session.extracellular.nSamples / fs; end
sigDur = win(2) - win(1);
w0 = win(1);
if ~isfinite(w0), w0 = 0; end

[~, ~, nremTimes] = evt_boutTimes(v, win, sigDur);
lfp = ripp_sigLoad(basepath, 'win', win, 'session', v.session, ...
    'basename', basename, 'rippCh', ripp.info.rippCh, 'bit2uv', []);
% rebuild the detection signal exactly as detection did, artifact mask included;
% a pre-260720 ripp.mat carries no otlThr, hence the default
otlThr = 8;
if isfield(ripp.info, 'otlThr'), otlThr = ripp.info.otlThr; end
rippSig = ripp_sigPrep(lfp, fs, 'detectMet', ripp.info.detectMet, ...
    'passband', ripp.info.passband, 'zMet', ripp.info.zMet, ...
    'nremTimes', nremTimes, 'otlThr', otlThr);

% rebuilt at the display width, so the detrend has only these flanks to work
% with - see the note above; a session with a saved rippMaps gets the better one
maps = evt_maps(rippSig, ripp.peakTime - w0, fs, 'mapDur', dispDur);
maps.lfp = evt_detrend(double(maps.lfp), maps.tstamps);
xt = maps.tstamps * 1000;               % ms

end     % loadMaps


function [maps, xt] = cropMaps(maps, dur)
% Keep the columns within DUR. The saved maps span the analyze window, which is
% wider than the view needs; cropping is a column index, not a recomputation.
keep = maps.tstamps >= dur(1) & maps.tstamps <= dur(2);
fn = fieldnames(maps);
for iFld = 1 : numel(fn)
    if strcmp(fn{iFld}, 'tstamps'), continue; end
    if isnumeric(maps.(fn{iFld})) && size(maps.(fn{iFld}), 2) == numel(keep)
        maps.(fn{iFld}) = maps.(fn{iFld})(:, keep);
    end
end
maps.tstamps = maps.tstamps(keep);
xt = maps.tstamps * 1000;               % ms

end     % cropMaps
