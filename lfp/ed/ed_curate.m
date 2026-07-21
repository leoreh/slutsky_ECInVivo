function [ed, hFig] = ed_curate(basepath, varargin)
% ED_CURATE Curate discharges by waveform TYPE, over an adjustable filter.
%
%   [ed, hFig] = ED_CURATE(basepath, varargin)
%
%   SUMMARY:
%       Stage 2 of the ED pipeline. Loads <basename>.ed.mat, applies the two
%       noise thresholds as a POOL, groups the pool into waveform clusters
%       (ed_clust), and lets you accept whole clusters. Saving writes
%       .accepted, .clustId and the choice into ed.info.
%
%       Why clusters. A 24 h recording proposes thousands of candidates and
%       holds a few dozen discharges. Judging that one event at a time is a day
%       of work; judging it by a MEAN waveform is worse than useless, because
%       an average over discharges, sharp waves and step artifacts is a curve
%       that is none of them - which is exactly how the discharges got buried
%       when this pipeline was first calibrated. Split by shape first and every
%       tile shows a real waveform.
%
%       Nothing here knows what a discharge looks like. The clustering is blind
%       and you name the clusters, so a mouse whose discharges differ from the
%       raMCU3/4/5 shape still gets them in a cluster of their own - which is
%       the point, since polarity and sharpness are layer-dependent.
%
%       WHAT IS ACCEPTED AND WHAT IS SHOWN ARE SEPARATE CONTROLS, deliberately.
%       Mixing them means you cannot inspect the events you rejected without
%       rejecting or accepting something by accident.
%
%       ACCEPT (checkboxes, left). Both start TICKED: curation here is
%       REJECTION, so everything the filter passed is accepted until you rule
%       something out.
%           clusters  which waveform types are discharges. Untick one whose
%                     median waveform is a sharp wave, a step or noise.
%           states    which vigilance states count. Untick one to drop a
%                     stretch of the recording wholesale, e.g. movement
%                     artifact in WAKE, without touching the shape decision.
%       An event is accepted when its cluster AND its state are ticked. That
%       is the mask Save writes.
%
%       SHOW (dropdown, left): 'both' | 'accepted' | 'removed'. Chooses which
%       rows reach the plot and nothing else - it cannot change the mask.
%
%       THE VIEW ITSELF is guiTbl_xy over those rows, carrying four variables
%       you can pivot on:
%           lfp      the waveform (Y)
%           cluster  shape group, or 'out' for events the filter dropped
%           state    vigilance state at the peak
%           status   accepted / removed under the current choice
%       So "Plot By (Tiles)" switches between a per-CLUSTER and a per-STATE
%       view, "Group By (Colors)" overlays the other, and Dispersion + Median
%       give a robust central trace rather than a mean.
%
%       REJECT THEN RE-CLUSTER IS A REFINEMENT LOOP. Re-cluster fits only the
%       events currently ACCEPTED (intersected with the thresholds, so the
%       knobs still bite). Having thrown out WAKE, or a cluster of step
%       artifacts, you do not want twelve groups spent describing events you
%       already rejected - you want twelve groups over what is left, which
%       splits the survivors finer each round. The accepted SET does not change
%       when you press it: the input is what you had accepted, and every new
%       cluster starts ticked.
%
%       Because it only narrows, 'Reset to filter' goes back to the whole pool
%       the two thresholds imply. Without it a mis-click would be
%       unrecoverable short of reopening the session.
%
%       The two thresholds are live knobs; the label updates as you type so you
%       can see the size you are choosing before paying for a fit.
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
%       260721 waveform clustering replaces the mean: accept types, not events.
%       260721b the knobs and the per-state view came back. Hand-drawn cluster
%              tiles were dropped for guiTbl_xy, which already does tiles,
%              grouping, category selection and a median-with-spread trace -
%              and which lets one view be pivoted instead of two being built.
%              The peri-event MUA panel went with them: a dozen events per
%              cluster is too few to read, and it shared no time axis with the
%              waveforms.
%       260721c accept and view split apart. Re-cluster used to silently clear
%              the accepted clusters; state became an acceptance criterion of
%              its own rather than only a way to tile; and which rows are drawn
%              moved to its own dropdown, so looking at what you rejected can
%              no longer change what you kept.
%       260721d clusters start ACCEPTED. Curation is rejection: the eye is much
%              better at spotting the two or three tiles that are obviously not
%              discharges than at confirming the ten that are.
%       260721e Re-cluster fits the ACCEPTED events, not the whole pool, so
%              rejecting and re-clustering refines: the groups no longer get
%              spent describing events already thrown out. 'Reset to filter'
%              undoes the narrowing.

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
    saveCurated(files.evt, pool, nan(numel(pool), 1), [], met.qa);
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
st.edK = gui_labeledControl(gCtrl, 'editnum', 'clusters', ...
    'Value', met.clust.nClust);
st.lblPool = gui_labeledControl(gCtrl, 'label', '');
gui_labeledControl(gCtrl, 'button', '', 'Text', 'Re-cluster', ...
    'ButtonPushedFcn', @(~,~) onCluster(hFig, false));
gui_labeledControl(gCtrl, 'button', '', 'Text', 'Reset to filter', ...
    'ButtonPushedFcn', @(~,~) onCluster(hFig, true));

% ACCEPT: a discharge is a cluster AND a state. Both start TICKED - curation
% here is rejection, so the pool is accepted until you rule a shape or a state
% out. Untick a cluster whose median waveform is not a discharge, or a state
% carrying something like movement artifact in WAKE.
st.gClust = gui_labeledControl(gCtrl, 'panel', 'accept clusters', ...
    'RowHeight', '1x');
st.gState = gui_labeledControl(gCtrl, 'panel', 'accept states', ...
    'RowHeight', 'fit');
st.lblKeep = gui_labeledControl(gCtrl, 'label', '');

% VIEW: independent of the above. Which rows reach the plot, nothing else.
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

buildStateChecks(hFig);
onCluster(hFig, true);

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
% Cluster, and rebuild the checkbox list and the view.
%
% Re-cluster fits only what is CURRENTLY ACCEPTED (intersected with the
% thresholds, so the knobs still bite). That is the useful move: having thrown
% out WAKE, or a cluster of step artifacts, you do not want twelve groups spent
% describing events you already rejected - you want twelve groups over what is
% left. Rejecting then re-clustering is therefore a refinement loop.
%
% It only ever narrows, so 'Reset to filter' goes back to the whole pool the
% thresholds imply. Without it a mis-click would be unrecoverable short of
% reopening the session.
st = hFig.UserData;
c  = st.met.clust;

gate = evt_gate(st.ed, buildSpec(st));
if flgReset
    st.pool = gate;
else
    st.pool = gate & st.acc;
end
iPool = find(st.pool);

k = st.edK.Value;
if k < 2, k = c.nClust; end

st.cid = nan(numel(st.pool), 1);
nClust = 0;
if ~isempty(iPool)
    [cid, cInfo] = ed_clust(st.wv(iPool, :), st.tst, 'win', c.win, ...
        'nPC', c.nPC, 'nClust', k, 'scalar', scalarFeat(st.ed, iPool));
    st.cid(iPool) = cid;
    nClust = cInfo.nClust;
end
st.nClust = nClust;
if nClust > 0, st.edK.Value = nClust; end

% Every cluster starts ticked. Carrying the old ticks over would be wrong
% twice: the partition is new, so index 3 no longer means what it did, and
% the input is BY CONSTRUCTION the set that was already accepted - so the
% accepted set is unchanged by the re-cluster, which is the property that
% matters. Rejections live in the input, not in the tick pattern.
st.selRestore = 1 : nClust;
st.lblPool.Text = sprintf('gate %d | clustered %d -> %d groups', ...
    nnz(gate), numel(iPool), nClust);
hFig.UserData = st;

buildChecks(hFig);
refresh(hFig);

end     % onCluster


function buildChecks(hFig)
% One checkbox per cluster, labelled with its size, restoring what was ticked
% when the count did not change. A pool too small to hold types (ed_clust
% returns nothing) simply gets none.
st = hFig.UserData;
delete(st.gClust.Children);
st.chk = gobjects(0);

if st.nClust > 0
    g = uigridlayout(st.gClust, [st.nClust + 1, 1], 'Padding', 2, ...
        'RowHeight', [repmat({'fit'}, 1, st.nClust), {'1x'}], ...
        'RowSpacing', 1, 'Scrollable', 'on');
    st.chk = gobjects(st.nClust, 1);
    for iK = 1 : st.nClust
        st.chk(iK) = uicheckbox(g, 'Value', ismember(iK, st.selRestore), ...
            'Text', sprintf('%d   (n = %d)', iK, nnz(st.cid == iK)), ...
            'ValueChangedFcn', @(~,~) refresh(hFig));
    end
end
hFig.UserData = st;

end     % buildChecks


function buildStateChecks(hFig)
% One checkbox per vigilance state, all ticked. Built ONCE - the states of a
% session do not change - so these ticks survive every re-cluster.
st = hFig.UserData;
cats = categories(removecats(st.state));

g = uigridlayout(st.gState, [numel(cats), 1], 'Padding', 2, ...
    'RowHeight', repmat({'fit'}, 1, numel(cats)), 'RowSpacing', 1);
st.chkState = gobjects(numel(cats), 1);
for iCat = 1 : numel(cats)
    st.chkState(iCat) = uicheckbox(g, 'Value', true, 'Text', ...
        sprintf('%s  (n = %d)', cats{iCat}, nnz(st.state == cats{iCat})), ...
        'ValueChangedFcn', @(~,~) refresh(hFig));
end
st.stateCats = cats;
hFig.UserData = st;

end     % buildStateChecks


function refresh(hFig)
% Recompute the kept mask and push the table to the view.
%
% The widget is built ONCE and fed rows thereafter. guiTbl_xy freezes its Y and
% Plot By / Group By item lists at construction, and the variable NAMES never
% change here - only the cluster categories do, which its setDataFcn
% reconciles by name. So even a re-cluster keeps whatever view the user set,
% which rebuilding would throw away on every press.
st = hFig.UserData;
accepted = acceptMask(st);
st.acc = accepted;              % what a Re-cluster will be fitted over
hFig.UserData = st;

st.lblKeep.Text = sprintf('accepted: %d of %d', nnz(accepted), ...
    numel(accepted));

% the SHOW dropdown chooses rows to draw and nothing else - it cannot change
% what is accepted, which is the whole point of separating the two
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
    guiTbl_xy(st.tst * 1000, tbl, 'Parent', st.hPanel, 'yVar', 'lfp', ...
        'tileVar', 'state', 'grpVar', 'cluster', 'xLbl', 'time (ms)');
end

end     % refresh


function accepted = acceptMask(st)
% An event is a discharge if its CLUSTER is ticked AND its STATE is ticked.
% Clusters say which shape; states are there to drop a whole stretch of the
% recording - movement artifact in WAKE, say - without touching the shape
% decision. Both start ticked, so the mask begins as the whole pool and
% curation removes from it.
accepted = ismember(st.cid, selectedClusters(st));
if isfield(st, 'chkState') && ~isempty(st.chkState)
    keepCat = st.stateCats(arrayfun(@(h) h.Value, st.chkState));
    accepted = accepted & ismember(st.state(:), keepCat);
end

end     % acceptMask


function doSave(hFig)
% Persist the mask, the labels and which clusters were chosen.
st = hFig.UserData;
sel = selectedClusters(st);
accepted = acceptMask(st);

saveCurated(st.files.evt, accepted, st.cid, sel, buildSpec(st));
buildStates(st.basepath, st.basename, st.ed, accepted);
gui_notify(hFig, sprintf('Saved: %d events from %d clusters (+ edStates)', ...
    nnz(accepted), numel(sel)), 'success');

end     % doSave


% =========================================================================
%  VIEW TABLE
% =========================================================================
function tbl = viewTable(st, accepted, iRow)
% The rows IROW asks for, with the three things worth pivoting on. Events the
% filter dropped keep the cluster label 'out' rather than being hidden, so a
% per-state view still shows what detection proposed.
lbl = [{'out'}, arrayfun(@(k) sprintf('%d', k), 1 : st.nClust, ...
    'uni', false)];
idx = st.cid;
idx(isnan(idx)) = 0;
clust = categorical(idx, 0 : st.nClust, lbl);

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
%  SPEC <-> CONTROLS
% =========================================================================
function qa = buildSpec(st)
qa.ranges = struct('fastZ', [st.edFast.Value, Inf], ...
    'isoZ', [st.edIso.Value, Inf]);

end     % buildSpec


function d = specDefaults(qa)
% The control values a spec implies; an absent bound is an open one.
d = struct('fastZ', -Inf, 'isoZ', -Inf);
if ~isfield(qa, 'ranges'), return; end
fn = fieldnames(d);
for iFld = 1 : numel(fn)
    if isfield(qa.ranges, fn{iFld})
        d.(fn{iFld}) = qa.ranges.(fn{iFld})(1);
    end
end

end     % specDefaults


function sel = selectedClusters(st)
% Cluster indices currently ticked.
sel = [];
if ~isempty(st.chk) && all(isgraphics(st.chk))
    sel = find(arrayfun(@(h) h.Value, st.chk));
end
sel = sel(:)';

end     % selectedClusters


function s = scalarFeat(ed, idx)
% The per-event shape measures that join the waveform components. They are
% shape descriptors already computed, so withholding them from the clustering
% would only throw information away.
s = [ed.fastZ(idx), ed.isoZ(idx), ed.posZ(idx), ed.amp(idx), ed.dur(idx)];

end     % scalarFeat


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


% =========================================================================
%  PERSISTENCE
% =========================================================================
function saveCurated(file, accepted, clustId, sel, qa)
% Back up, then overwrite the mask + the labels in the saved struct.
backup_file(file);
S = load(file);
S.ed.accepted      = logical(accepted(:));
S.ed.clustId       = clustId(:);
S.ed.info.qa       = qa;
S.ed.info.clustSel = sel;
save(file, '-struct', 'S', '-v7.3');

end     % saveCurated


function buildStates(basepath, basename, ed, accepted)
% Rebuild + save the per-bout rate table for the current mask (cheap; no
% signal). Skips silently when sleep states are unavailable.
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
