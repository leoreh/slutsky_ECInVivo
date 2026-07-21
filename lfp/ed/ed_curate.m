function [ed, hFig] = ed_curate(basepath, varargin)
% ED_CURATE Curate discharges by waveform TYPE, over an adjustable filter.
%
%   [ed, hFig] = ED_CURATE(basepath, varargin)
%
%   SUMMARY:
%       Stage 2 of the ED pipeline. Loads <basename>.ed.mat, applies the two
%       noise thresholds as a POOL, groups the pool into waveform clusters
%       (ed_clust) and lets you accept whole clusters. Saving writes .accepted,
%       .clustId and the choice into ed.info.
%
%       A 24 h recording proposes thousands of candidates and holds a few dozen
%       discharges, so the unit of curation here is a TYPE, not an event.
%
%       WHAT IS ACCEPTED AND WHAT IS SHOWN ARE SEPARATE CONTROLS. Mixing them
%       means you cannot inspect the events you rejected without rejecting or
%       accepting something by accident.
%
%       ACCEPT (checkboxes). An event is accepted when its CLUSTER and its
%       STATE are both ticked - that is the mask Save writes. Both lists start
%       TICKED: curation here is REJECTION, so everything the filter passed is
%       accepted until you rule something out.
%           clusters  which waveform types are discharges. Untick one whose
%                     median waveform is a sharp wave, a step or noise.
%           states    which vigilance states count. Untick one to drop a
%                     stretch of the recording wholesale (movement artifact in
%                     WAKE, say) without touching the shape decision.
%
%       SHOW (dropdown): 'both' | 'accepted' | 'removed'. Chooses which rows
%       reach the plot and nothing else - it cannot change the mask.
%
%       THE VIEW is guiTbl_xy over those rows, carrying four variables to pivot
%       on: lfp (the waveform, Y), cluster, state and status. So "Plot By
%       (Tiles)" switches between a per-CLUSTER and a per-STATE view, "Group By
%       (Colors)" overlays the other, and Dispersion + Median give a robust
%       central trace rather than a mean.
%
%       The two thresholds are live knobs: the label updates as you type, so a
%       pool can be sized before paying for a fit. Reject then Re-cluster is a
%       refinement loop (see onCluster), and reopening resumes the saved
%       curation rather than clustering afresh (see restoreSaved). The cluster
%       count follows the pool unless the box overrides it (0 = auto), so a
%       refinement round over a narrowed pool gets fewer, tighter groups.
%
%       Headless (flgGui = false) applies only the thresholds, for a batch run
%       that has no human. That mask is NOT an answer - it is the pool.
%
%   INPUTS:
%       basepath - <char> session directory (must hold <basename>.ed.mat).
%       varargin - Parameter/Value:
%           'basename' - <char>   file stem. {folder name}
%           'met'      - <struct> config; reads .qa and .clust.
%                                 {ed_methods('default')}
%           'flgGui'   - <log>    open the GUI (true) or filter headless.{true}
%           'Visible'  - <char>   'on' | 'off', for headless GUI tests. {'on'}
%           'verbose'  - <log>    print progress? {true}
%
%   OUTPUTS:
%       ed   - <struct> the loaded events (headless: with .accepted updated).
%       hFig - <handle> the GUI figure ([] when headless).
%
%   DEPENDENCIES:
%       evt_files, evt_gate, evt_states, evt_boutTimes, ed_methods, ed_clust,
%       basepaths2vars, backup_file; GUI: gui_layout, gui_labeledControl,
%       gui_notify, guiTbl_xy.
%
%   HISTORY:
%       260720 created as the ED twin of ripp_curate (threshold knobs over a
%              kept-vs-removed mean waveform).
%       260721 rebuilt on waveform clustering: accept TYPES, not events. The
%              view is guiTbl_xy, which already does tiles, grouping and a
%              median-with-spread trace, so one pivotable view replaces two
%              hand-drawn ones. Accept and view became separate controls, every
%              cluster starts accepted, Re-cluster refines by fitting only the
%              accepted events, and reopening resumes the saved curation. See
%              dev/ed_pipeline_rebuild.md.

%% ========================================================================
%  ARGUMENTS + LOAD
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'met', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'flgGui', true, @islogical);
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(char(x), {'on', 'off'})));
addParameter(p, 'verbose', true, @islogical);
parse(p, basepath, varargin{:});
met     = p.Results.met;
flgGui  = p.Results.flgGui;
vis     = char(p.Results.Visible);
verbose = p.Results.verbose;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end
if isempty(met), met = ed_methods('default'); end

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
    pool = evt_gate(ed, met.qa);
    ed.accepted = pool;
    saveCurated(files.evt, pool, nan(numel(pool), 1), [], met.qa, {});
    buildStates(basepath, basename, ed, pool);
    if verbose
        fprintf('[ED_CURATE] %s : %d / %d pass the filter\n', ...
            basename, nnz(pool), numel(pool));
    end
    return;
end

%% ========================================================================
%  GUI
%  ========================================================================
st = struct();
st.ed       = ed;
st.met      = met;
st.files    = files;
st.basepath = basepath;
st.basename = basename;
st.state    = plotState(ed);
[st.wv, st.tst] = loadMaps(files.maps, ed);

hFig = uifigure('Name', ['ED curation: ' basename], ...
    'Position', [60 60 1600 850], 'Visible', vis);
[~, gPlot, gCtrl, gActions] = gui_layout(hFig, 'CtrlWidth', 240);

dflt = specDefaults(met.qa);
st.edFast = gui_labeledControl(gCtrl, 'editnum', 'sharp  fastZ >=', ...
    'Value', dflt.fastZ, 'ValueChangedFcn', @(~,~) onKnob(hFig));
st.edIso = gui_labeledControl(gCtrl, 'editnum', 'alone  isoZ >=', ...
    'Value', dflt.isoZ, 'ValueChangedFcn', @(~,~) onKnob(hFig));
kInit = met.clust.nClust;
if isempty(kInit), kInit = 0; end       % 0 = scale with the pool
st.edK = gui_labeledControl(gCtrl, 'editnum', 'clusters (0 = auto)', ...
    'Value', kInit);
st.lblPool = gui_labeledControl(gCtrl, 'label', '');
gui_labeledControl(gCtrl, 'button', '', 'Text', 'Re-cluster', ...
    'ButtonPushedFcn', @(~,~) onCluster(hFig, false));
gui_labeledControl(gCtrl, 'button', '', 'Text', 'Reset to filter', ...
    'ButtonPushedFcn', @(~,~) onCluster(hFig, true));

st.gClust = gui_labeledControl(gCtrl, 'panel', 'accept clusters', ...
    'RowHeight', '1x');
st.gState = gui_labeledControl(gCtrl, 'panel', 'accept states', ...
    'RowHeight', 'fit');
st.lblKeep = gui_labeledControl(gCtrl, 'label', '');

% independent of the above: which rows reach the plot, nothing else
st.ddShow = gui_labeledControl(gCtrl, 'dropdown', 'show', ...
    'Items', {'both', 'accepted', 'removed'}, 'Value', 'accepted', ...
    'ValueChangedFcn', @(~,~) refresh(hFig));

gui_labeledControl(gActions, 'button', '', 'Text', 'Save', ...
    'ButtonPushedFcn', @(~,~) doSave(hFig));

st.hPanel = uipanel(gPlot, 'BorderType', 'none');
st.chk    = gobjects(0);
st.cid    = nan(numel(ed.peakTime), 1);
st.nClust = 0;
st.acc    = true(numel(ed.peakTime), 1);    % nothing rejected yet
hFig.UserData = st;

buildStateChecks(hFig, prevStates(ed));
if ~restoreSaved(hFig, ed)
    onCluster(hFig, true);
end

end     % EOF


% =========================================================================
%  GUI CALLBACKS
% =========================================================================
function onKnob(hFig)
% A threshold moved: show the pool it implies. The clustering is NOT redone -
% a fit over a few hundred events is not something to run on every keystroke,
% and the point of the knob is to choose a pool size before paying for it.
st = hFig.UserData;
gate = evt_gate(st.ed, buildSpec(st));
st.lblPool.Text = sprintf('gate %d | %d accepted  (press Re-cluster)', ...
    nnz(gate), nnz(gate & st.acc));

end     % onKnob


function onCluster(hFig, flgReset)
% Fit waveform clusters over the pool and adopt the result.
%
% Re-cluster fits only what is CURRENTLY ACCEPTED, intersected with the
% thresholds so the knobs still bite. Having thrown out WAKE, or a cluster of
% step artifacts, you do not want the groups spent describing events already
% rejected - you want them over what is left, which splits the survivors finer
% each round. The accepted SET does not move: the input IS what was accepted,
% and every new cluster starts ticked, so a rejection lives in the input rather
% than in the tick pattern. Carrying the old ticks over would be wrong anyway -
% the partition is new, so index 3 no longer means what it did.
%
% Because it only ever narrows, FLGRESET goes back to the whole pool the
% thresholds imply; without it a mis-click would be unrecoverable short of
% reopening the session.
st = hFig.UserData;
c  = st.met.clust;

gate = evt_gate(st.ed, buildSpec(st));
st.pool = gate;
if ~flgReset, st.pool = gate & st.acc; end
iPool = find(st.pool);

k = st.edK.Value;
if k < 2, k = c.nClust; end     % 0 in the box, or an empty default: auto

st.cid = nan(numel(st.pool), 1);
st.nClust = 0;
if ~isempty(iPool)
    % the per-event measures join the waveform components: they are shape
    % descriptors already computed, so withholding them from the clustering
    % would only throw information away
    scalar = [st.ed.fastZ(iPool), st.ed.isoZ(iPool), st.ed.posZ(iPool), ...
        st.ed.amp(iPool), st.ed.dur(iPool)];
    [cid, cInfo] = ed_clust(st.wv(iPool, :), st.tst, 'win', c.win, ...
        'nPC', c.nPC, 'nClust', k, 'scalar', scalar);
    st.cid(iPool) = cid;
    st.nClust = cInfo.nClust;
end
hFig.UserData = st;

adoptClust(hFig, 1 : st.nClust, ...
    sprintf('gate %d | clustered %d -> %d groups', nnz(gate), ...
    numel(iPool), st.nClust));

end     % onCluster


function ok = restoreSaved(hFig, ed)
% Put back the partition and the choices a previous session saved, instead of
% clustering afresh. The saved labels ARE the partition that was judged: a
% re-fit would cost a fit and hand back a different one, leaving the ticks
% pointing at groups nobody looked at. Anything that no longer lines up with
% the event list is refused, so a re-detection falls through to a fresh fit
% rather than drawing labels that belong to other events.
ok = isfield(ed, 'clustId') && numel(ed.clustId) == numel(ed.peakTime) ...
    && ~all(isnan(ed.clustId));
if ~ok, return; end

st = hFig.UserData;
st.cid    = ed.clustId(:);
st.pool   = ~isnan(st.cid);
st.nClust = max(st.cid);

sel = 1 : st.nClust;
if isfield(ed.info, 'clustSel') && ~isempty(ed.info.clustSel)
    sel = ed.info.clustSel;
end

dflt = specDefaults(ed.info.qa);
st.edFast.Value = dflt.fastZ;
st.edIso.Value  = dflt.isoZ;
hFig.UserData = st;

adoptClust(hFig, sel, sprintf('restored: %d clustered -> %d groups', ...
    nnz(st.pool), st.nClust));

end     % restoreSaved


function adoptClust(hFig, sel, msg)
% Take on the labelling now in UserData: report it, rebuild the cluster
% checkboxes with SEL ticked, and redraw. Shared by a fresh fit and a restore,
% which differ only in where the labels came from.
%
% The count box is NOT touched. It holds the count the user ASKED for, and the
% label reports the count that came back - writing the resolved count into the
% box would turn "0 = auto" into a fixed number after the first fit, so a
% refinement round over a narrowed pool would keep splitting it into as many
% groups as the pool it came from.
st = hFig.UserData;
st.lblPool.Text = msg;

% one checkbox per cluster, labelled with its size. A pool too small to hold
% types (ed_clust returns no labels) simply gets none.
delete(st.gClust.Children);
st.chk = gobjects(0);
if st.nClust > 0
    g = uigridlayout(st.gClust, [st.nClust + 1, 1], 'Padding', 2, ...
        'RowHeight', [repmat({'fit'}, 1, st.nClust), {'1x'}], ...
        'RowSpacing', 1, 'Scrollable', 'on');
    st.chk = gobjects(st.nClust, 1);
    for iK = 1 : st.nClust
        st.chk(iK) = uicheckbox(g, 'Value', ismember(iK, sel), ...
            'Text', sprintf('%d   (n = %d)', iK, nnz(st.cid == iK)), ...
            'ValueChangedFcn', @(~,~) refresh(hFig));
    end
end
hFig.UserData = st;

refresh(hFig);

end     % adoptClust


function refresh(hFig)
% Recompute the accepted mask and push the chosen rows to the view.
%
% The widget is built ONCE and fed rows thereafter. guiTbl_xy freezes its Y and
% Plot By / Group By item lists at construction, and the variable NAMES never
% change here - only the cluster categories do, which its setDataFcn reconciles
% by name. So even a re-cluster keeps whatever view the user set, which
% rebuilding would throw away on every press.
st = hFig.UserData;
accepted = acceptMask(st);
st.acc = accepted;              % what a Re-cluster will be fitted over
hFig.UserData = st;

st.lblKeep.Text = sprintf('accepted: %d of %d', nnz(accepted), ...
    numel(accepted));

switch st.ddShow.Value
    case 'accepted', iRow = find(accepted);
    case 'removed',  iRow = find(~accepted);
    otherwise,       iRow = (1 : numel(accepted))';
end
tbl = viewTable(st, accepted, iRow);

ud = st.hPanel.UserData;
if isstruct(ud) && isfield(ud, 'setDataFcn')
    ud.setDataFcn(tbl);
elseif ~isempty(iRow)
    % first call. An empty table has no categories to build the filter panels
    % from, so the widget waits until there is something to draw.
    guiTbl_xy(st.tst * 1000, tbl, 'Parent', st.hPanel, 'yVar', 'lfp', ...
        'tileVar', 'state', 'grpVar', 'cluster', 'xLbl', 'time (ms)');
end

end     % refresh


function doSave(hFig)
% Persist the mask plus everything needed to resume it.
st = hFig.UserData;
[accepted, sel, states] = acceptMask(st);

saveCurated(st.files.evt, accepted, st.cid, sel, buildSpec(st), states);
buildStates(st.basepath, st.basename, st.ed, accepted);
gui_notify(hFig, sprintf('Saved: %d events from %d clusters (+ edStates)', ...
    nnz(accepted), numel(sel)), 'success');

end     % doSave


% =========================================================================
%  VIEW
% =========================================================================
function tbl = viewTable(st, accepted, iRow)
% The rows IROW asks for, with the three things worth pivoting on. Events the
% filter dropped keep the cluster label 'out' rather than being hidden, so a
% per-state view still shows what detection proposed.
idx = st.cid;
idx(isnan(idx)) = 0;
clust = categorical(idx, 0 : st.nClust, ...
    [{'out'}, cellstr(string(1 : st.nClust))]);

tbl = table(st.wv(iRow, :), clust(iRow), st.state(iRow), ...
    categorical(accepted(iRow), [false true], {'removed', 'kept'}), ...
    'VariableNames', {'lfp', 'cluster', 'state', 'status'});

end     % viewTable


function s = plotState(ed)
% ed.state as a tiling variable: <undefined> is promoted to its own 'unscored'
% level. A categorical comparison never matches <undefined>, so without this
% the unscored events would get no tile and vanish from the view without a
% word. A session that was never scored has no .state at all.
if ~isfield(ed, 'state') || isempty(ed.state)
    s = categorical(repmat({'unscored'}, numel(ed.peakTime), 1));
    return;
end
s = removecats(ed.state(:));
if any(isundefined(s))
    s = addcats(s, {'unscored'});
    s(isundefined(s)) = 'unscored';
end

end     % plotState


% =========================================================================
%  CONTROLS
% =========================================================================
function [accepted, sel, states] = acceptMask(st)
% An event is a discharge if its CLUSTER is ticked AND its STATE is ticked.
% Clusters say which shape; states are there to drop a whole stretch of the
% recording - movement artifact in WAKE, say - without touching the shape
% decision. Both start ticked, so the mask begins as the whole pool and
% curation removes from it. SEL and STATES come back with the mask because
% Save records the choice as well as its result; deriving them apart is how
% the two drift.
sel    = find(ticked(st.chk));
states = st.stateCats(ticked(st.chkState));
states = states(:)';
accepted = ismember(st.cid, sel) & ismember(st.state(:), states);

end     % acceptMask


function buildStateChecks(hFig, keepCats)
% One checkbox per vigilance state. Built ONCE - the states of a session do not
% change - so these ticks survive every re-cluster. KEEPCATS restores a saved
% selection; [] (never curated) means all ticked, an empty CELL means none.
st = hFig.UserData;
cats = categories(removecats(st.state));

g = uigridlayout(st.gState, [numel(cats), 1], 'Padding', 2, ...
    'RowHeight', repmat({'fit'}, 1, numel(cats)), 'RowSpacing', 1);
st.chkState = gobjects(numel(cats), 1);
for iCat = 1 : numel(cats)
    val = ~iscell(keepCats) || ismember(cats{iCat}, keepCats);
    st.chkState(iCat) = uicheckbox(g, 'Value', val, 'Text', ...
        sprintf('%s  (n = %d)', cats{iCat}, nnz(st.state == cats{iCat})), ...
        'ValueChangedFcn', @(~,~) refresh(hFig));
end
st.stateCats = cats;
hFig.UserData = st;

end     % buildStateChecks


function tf = ticked(h)
% Which checkboxes of an array are ticked, as a logical row.
tf = false(1, numel(h));
for iChk = 1 : numel(h)
    tf(iChk) = h(iChk).Value;
end

end     % ticked


function qa = buildSpec(st)
% The filter spec the two knobs currently describe.
qa.ranges = struct('fastZ', [st.edFast.Value, Inf], ...
    'isoZ', [st.edIso.Value, Inf]);

end     % buildSpec


function d = specDefaults(qa)
% The knob values a spec implies; an absent bound is an open one.
d = struct('fastZ', -Inf, 'isoZ', -Inf);
if ~isfield(qa, 'ranges'), return; end
if isfield(qa.ranges, 'fastZ'), d.fastZ = qa.ranges.fastZ(1); end
if isfield(qa.ranges, 'isoZ'),  d.isoZ  = qa.ranges.isoZ(1);  end

end     % specDefaults


% =========================================================================
%  DATA
% =========================================================================
function [wv, tst] = loadMaps(file, ed)
% Per-event waveforms behind the clustering and the view.
if ~isfile(file)
    error('ed_curate:noMaps', ...
        'no edMaps file; re-run detection with flgSave.');
end
S = load(file, 'edMaps');
if size(S.edMaps.lfp, 1) ~= numel(ed.peakTime)
    error('ed_curate:staleMaps', ...
        'edMaps does not match the event list; re-run detection.');
end
wv  = double(S.edMaps.lfp);
tst = S.edMaps.tstamps;

end     % loadMaps


function cats = prevStates(ed)
% The state selection a previous session saved: a cellstr when one exists, []
% when none does. The distinction matters - an empty CELL means "every state
% was rejected", which is not the same as "never curated".
cats = [];
if isfield(ed, 'info') && isfield(ed.info, 'clustStates') ...
        && iscell(ed.info.clustStates)
    cats = ed.info.clustStates;
end

end     % prevStates


% =========================================================================
%  PERSISTENCE
% =========================================================================
function saveCurated(file, accepted, clustId, sel, qa, states)
% Back up, then overwrite the mask + everything needed to resume: the labels,
% which clusters and states were accepted, and the thresholds behind them.
backup_file(file);
S = load(file);
S.ed.accepted         = logical(accepted(:));
S.ed.clustId          = clustId(:);
S.ed.info.qa          = qa;
S.ed.info.clustSel    = sel;
S.ed.info.clustStates = states;
save(file, '-struct', 'S', '-v7.3');

end     % saveCurated


function buildStates(basepath, basename, ed, accepted)
% Rebuild + save the per-bout rate table for the current mask (cheap; no
% signal). Skips silently when sleep states are unavailable.
win = ed.info.win;
w0  = win(1);
if ~isfinite(w0), w0 = 0; end

v = basepaths2vars('basepaths', {basepath}, 'vars', {'sleep_states'});
boutTimes = evt_boutTimes(v, win, win(2) - win(1));
if isempty(boutTimes)
    return;
end
evt_states(ed.times - w0, ed.peakTime - w0, boutTimes, ...
    'accepted', logical(accepted(:)), 'basepath', basepath, ...
    'basename', basename, 'flgSave', true, 'flgPlot', false, ...
    'name', 'ed', 'lbl', 'ED');

end     % buildStates
