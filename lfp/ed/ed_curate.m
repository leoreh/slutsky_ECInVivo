function [ed, hFig] = ed_curate(basepath, varargin)
% ED_CURATE Post-detection QA gate for EDs: headless or interactive (stage 2).
%
%   [ed, hFig] = ED_CURATE(basepath, varargin)
%
%   SUMMARY:
%       Stage 2 of the ED pipeline. Loads <basename>.ed.mat and sets the
%       per-event .accepted mask from a QA filter - per-metric [lo hi] ranges.
%       The gate itself is evt_gate, shared with the ripple pipeline; this is
%       the two ways to drive it:
%
%       - Headless (flgGui = false): apply the spec, save .accepted and the
%         spec in ed.info.qa, rebuild the per-bout rate table. The automatic
%         gate, so a batch run needs no human.
%       - Interactive (default): the same spec as thresholds. The counts and
%         the kept-vs-removed mean waveform, tiled by state, update live as you
%         move them. Three knobs carry the decision - is it sharp (fastZ), does
%         it go up (posZ), does it stand alone (isoZ) - and EMG is there for a
%         session where movement artifact is the problem.
%
%       Read the waveform view knowing what it is: a MEAN. It only shows the
%       discharge shape once the kept set is mostly discharges. If the counts
%       are in the hundreds the average is whatever the bulk happens to be, and
%       the tile says nothing - raise fastZ until the count is plausible for a
%       day of recording (tens), then judge the shape.
%
%       This is the BULK pass, and it is what makes the per-event pass
%       possible: detection is permissive by design and a 24 h recording yields
%       thousands of candidates, far too many to step through one at a time.
%       Set the thresholds here, then walk the survivors in guiPath.
%
%       Nothing is destroyed - a rejected event keeps its row with .accepted =
%       false, and Reset restores the default spec. Saving touches only
%       .accepted, ed.info.qa and edStates, and backs the file up first. There
%       is no invalidation step as ripples have: edMaps holds one row per
%       DETECTED event and is row-aligned to ed, so a mask cannot stale it.
%
%   INPUTS:
%       basepath - <char> session directory (must hold <basename>.ed.mat).
%       varargin - Parameter/Value:
%           'basename' - <char>   file stem. {folder name}
%           'qa'       - <struct> filter spec (see ed_methods '.qa').
%                                 {ed_methods('default').qa}
%           'flgGui'   - <log>    open the GUI (true) or apply headless.{true}
%           'Visible'  - <char>   'on' | 'off', for headless GUI tests. {'on'}
%           'verbose'  - <log>    print progress? {true}
%
%   OUTPUTS:
%       ed   - <struct> the loaded events (headless: with .accepted updated).
%       hFig - <handle> the GUI figure ([] when headless).
%
%   DEPENDENCIES:
%       evt_files, evt_gate, ed_methods, evt_states, evt_boutTimes,
%       basepaths2vars, backup_file; GUI: gui_layout, gui_labeledControl,
%       gui_notify, guiTbl_xy.
%
%   HISTORY:
%       260720 created as the ED twin of ripp_curate, replacing the destructive
%              evt_qa + evt_subset filter the old ed_wrapper ran at detection.
%              No state filter: ED asks how discharges distribute over states,
%              so restricting them is a question for ed_tbl, not for curation.

%% ========================================================================
%  ARGUMENTS + LOAD
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'qa', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'flgGui', true, @islogical);
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(char(x), {'on', 'off'})));
addParameter(p, 'verbose', true, @islogical);
parse(p, basepath, varargin{:});
qa      = p.Results.qa;
flgGui  = p.Results.flgGui;
vis     = char(p.Results.Visible);
verbose = p.Results.verbose;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end
if isempty(qa), qa = ed_methods('default').qa; end

files = evt_files(basepath, basename, 'ed');
if ~isfile(files.evt)
    error('ed_curate:noFile', ...
        '%s not found; run detection first (ed_wrapper).', files.evt);
end
S  = load(files.evt, 'ed');
ed = S.ed;
hFig = [];

%% ========================================================================
%  HEADLESS
%  ========================================================================
if ~flgGui
    ed.accepted = evt_gate(ed, qa);
    saveCurated(files.evt, ed.accepted, qa);
    buildStates(basepath, basename, ed, ed.accepted);
    if verbose
        fprintf('[ED_CURATE] %s : %d / %d accepted (headless)\n', ...
            basename, nnz(ed.accepted), numel(ed.accepted));
    end
    return;
end

%% ========================================================================
%  GUI
%  ========================================================================
hFig = uifigure('Name', ['ED curation: ' basename], ...
    'Position', [80 80 1500 850], 'Visible', vis);
[~, gPlot, gCtrl, gActions] = gui_layout(hFig, 'CtrlWidth', 250);

dflt = specDefaults(qa);
st = struct();
st.edFast = gui_labeledControl(gCtrl, 'editnum', 'sharp  fastZ >=', ...
    'Value', dflt.fastZ, 'ValueChangedFcn', @(~,~) refresh(hFig));
st.edPos = gui_labeledControl(gCtrl, 'editnum', 'upward posZ >=', ...
    'Value', dflt.posZ, 'ValueChangedFcn', @(~,~) refresh(hFig));
st.edIso = gui_labeledControl(gCtrl, 'editnum', 'alone  isoZ >=', ...
    'Value', dflt.isoZ, 'ValueChangedFcn', @(~,~) refresh(hFig));
st.edEmg = gui_labeledControl(gCtrl, 'editnum', 'EMG <=', ...
    'Value', dflt.emg, 'ValueChangedFcn', @(~,~) refresh(hFig));
st.lblCount = gui_labeledControl(gCtrl, 'label', '');

gui_labeledControl(gActions, 'button', '', 'Text', 'Reset to default', ...
    'ButtonPushedFcn', @(~,~) onReset(hFig));
gui_labeledControl(gActions, 'button', '', 'Text', 'Save', ...
    'ButtonPushedFcn', @(~,~) doSave(hFig));

st.hPanel   = uipanel(gPlot, 'BorderType', 'none');
st.ed       = ed;
st.qa0      = qa;
st.files    = files;
st.basepath = basepath;
st.basename = basename;
st.accepted = ed.accepted;
st.stateCol = plotState(ed.state);      % tiling variable; fixed for the session

% the waveform view is best-effort: without the detect-stage maps the
% thresholds still work, they just have no picture behind them
st.maps = [];
st.xt   = [];
try
    [st.maps, st.xt] = loadMaps(files.maps, ed);
catch ME
    warning('ed_curate:maps', 'waveform view disabled (%s)', ME.message);
end

hFig.UserData = st;
refresh(hFig);

end     % EOF


% =========================================================================
%  GUI CALLBACKS
% =========================================================================
function refresh(hFig)
% Recompute accepted from the current thresholds; redraw counts + waveform.
st = hFig.UserData;
st.accepted = evt_gate(st.ed, buildSpec(st));
hFig.UserData = st;

st.lblCount.Text = sprintf('kept %d / %d', nnz(st.accepted), ...
    numel(st.accepted));

% Only the kept / removed flag changes as a threshold moves, so the widget is
% updated in place - rebuilding it would reset the Y / tile / group selections
% on every keystroke.
if ~isempty(st.maps)
    status = repmat("removed", numel(st.accepted), 1);
    status(st.accepted) = "kept";
    tbl = table(st.maps.lfp, categorical(status, {'removed', 'kept'}), ...
        st.stateCol, 'VariableNames', {'lfp', 'status', 'state'});

    ud = st.hPanel.UserData;
    if isstruct(ud) && isfield(ud, 'setDataFcn')
        ud.setDataFcn(tbl);
    else
        guiTbl_xy(st.xt, tbl, 'Parent', st.hPanel, 'yVar', 'lfp', ...
            'tileVar', 'state', 'grpVar', 'status', 'xLbl', 'time (ms)');
    end
end

end     % refresh


function onReset(hFig)
% Restore the thresholds to the default spec.
st = hFig.UserData;
dflt = specDefaults(st.qa0);
st.edFast.Value = dflt.fastZ;
st.edPos.Value  = dflt.posZ;
st.edIso.Value  = dflt.isoZ;
st.edEmg.Value  = dflt.emg;
refresh(hFig);

end     % onReset


function doSave(hFig)
% Persist accepted + the spec, and rebuild the per-bout rate table.
st = hFig.UserData;
saveCurated(st.files.evt, st.accepted, buildSpec(st));
buildStates(st.basepath, st.basename, st.ed, st.accepted);
gui_notify(hFig, sprintf('Saved: %d / %d accepted (+ edStates)', ...
    nnz(st.accepted), numel(st.accepted)), 'success');

end     % doSave


% =========================================================================
%  SPEC <-> CONTROLS
% =========================================================================
function qa = buildSpec(st)
qa.ranges = struct('fastZ', [st.edFast.Value, Inf], ...
    'posZ', [st.edPos.Value, Inf], ...
    'isoZ', [st.edIso.Value, Inf], ...
    'emg', [-Inf, st.edEmg.Value]);

end     % buildSpec


function d = specDefaults(qa)
% The control values a spec implies; an absent bound is an open one.
d = struct('fastZ', -Inf, 'posZ', -Inf, 'isoZ', -Inf, 'emg', Inf);
if ~isfield(qa, 'ranges'), return; end
lo = {'fastZ', 'posZ', 'isoZ'};
for iFld = 1 : numel(lo)
    if isfield(qa.ranges, lo{iFld}), d.(lo{iFld}) = qa.ranges.(lo{iFld})(1); end
end
if isfield(qa.ranges, 'emg'), d.emg = qa.ranges.emg(2); end

end     % specDefaults


% =========================================================================
%  PERSISTENCE
% =========================================================================
function saveCurated(file, accepted, qa)
% Back up, then overwrite .accepted + .info.qa in the saved struct.
backup_file(file);
S = load(file);
S.ed.accepted = logical(accepted(:));
S.ed.info.qa  = qa;
save(file, '-struct', 'S', '-v7.3');

end     % saveCurated


function buildStates(basepath, basename, ed, accepted)
% Rebuild + save the per-bout rate / density table for the current mask (cheap;
% no signal). Skips silently when sleep states are unavailable.
win = ed.info.win;
w0  = win(1);
if ~isfinite(w0), w0 = 0; end
if isinf(win(2)), winDur = Inf; else, winDur = win(2) - win(1); end

v = basepaths2vars('basepaths', {basepath}, 'vars', {'sleep_states'});
boutTimes = evt_boutTimes(v, win, winDur);
if isempty(boutTimes)
    return;
end
evt_states(ed.times - w0, ed.peakTime - w0, boutTimes, ...
    'accepted', logical(accepted(:)), 'basepath', basepath, ...
    'basename', basename, 'flgSave', true, 'flgPlot', false, ...
    'name', 'ed', 'lbl', 'ED');

end     % buildStates


% =========================================================================
%  WAVEFORMS
% =========================================================================
function s = plotState(state)
% ed.state as a tiling variable: <undefined> is promoted to its own 'unscored'
% level. A categorical comparison never matches <undefined>, so without this the
% unscored events would get no tile and vanish from the view without a word.
s = removecats(state(:));
if any(isundefined(s))
    s = addcats(s, {'unscored'});
    s(isundefined(s)) = 'unscored';
end

end     % plotState


function [maps, xt] = loadMaps(file, ed)
% The per-event LFP maps behind the waveform view, cropped to DISPDUR. A file
% that no longer matches the event list is refused rather than silently drawing
% the wrong waveforms.
DISPDUR = [-0.1 0.1];
if ~isfile(file)
    error('no edMaps file; re-run detection with flgSave');
end
S = load(file, 'edMaps');
if ~isfield(S, 'edMaps') || size(S.edMaps.lfp, 1) ~= numel(ed.peakTime)
    error('edMaps does not match the event list; re-run detection');
end

maps = S.edMaps;
keep = maps.tstamps >= DISPDUR(1) & maps.tstamps <= DISPDUR(2);
fn = fieldnames(maps);
for iFld = 1 : numel(fn)
    if ~strcmp(fn{iFld}, 'tstamps') && size(maps.(fn{iFld}), 2) == numel(keep)
        maps.(fn{iFld}) = maps.(fn{iFld})(:, keep);
    end
end
maps.tstamps = maps.tstamps(keep);
xt = maps.tstamps * 1000;               % ms

end     % loadMaps
