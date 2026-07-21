function [accepted, hFig] = evt_curate(evt, wv, tstamps, cfg)
% EVT_CURATE Curate events by waveform TYPE: accept or reject whole shapes.
%
%   [accepted, hFig] = EVT_CURATE(evt, wv, tstamps, cfg)
%
%   SUMMARY:
%       The curation stage shared by both event pipelines (ed_curate,
%       ripp_curate). It applies the metric filter as a POOL, groups the pool
%       into waveform clusters (evt_clust), and lets a human accept whole
%       clusters. Saving writes .accepted, .clustId and the choice into
%       <var>.info.
%
%       Detection proposes thousands of events and a session holds one or two
%       populations worth keeping, so the unit of curation is a TYPE, not an
%       event. Judging by a MEAN waveform is worse than useless - an average
%       over ripples, sharp transients and step artifacts is a curve that is
%       none of them. Split by shape first and every tile shows a real shape.
%
%       WHAT IS ACCEPTED AND WHAT IS SHOWN ARE SEPARATE CONTROLS. Mixing them
%       means you cannot inspect the events you rejected without rejecting or
%       accepting something by accident.
%
%       ACCEPT (checkboxes). An event is accepted when its CLUSTER and its
%       STATE are both ticked - that is the mask Save writes. Both lists start
%       TICKED: curation here is REJECTION, so everything the filter passed is
%       accepted until you rule something out.
%           clusters  which waveform types are real. Untick one whose median
%                     waveform is a step, a spike-bleed transient or noise. A
%                     shape rejection STICKS across a Re-cluster - undo it with
%                     'Reset to filter'.
%           states    which vigilance states count. Untick one to drop a
%                     stretch of the recording wholesale (movement artifact in
%                     WAKE, say) without touching the shape decision. State is
%                     a SCOPE, not a verdict: re-tick it and its events are
%                     eligible again on the next Re-cluster.
%
%       SHOW (dropdown): 'both' | 'accepted' | 'removed'. Chooses which rows
%       reach the plot and nothing else - it cannot change the mask.
%
%       THE VIEW is guiTbl_xy over those rows, carrying every waveform passed
%       in as a selectable Y plus three things to pivot on: cluster, state and
%       status. So "Y Var" switches between the raw and the filtered trace,
%       "Plot By (Tiles)" between a per-CLUSTER and a per-STATE view, "Group By
%       (Colors)" overlays the other, and Dispersion + Median give a robust
%       central trace rather than a mean.
%
%       The metric knobs are live: the label updates as you type, so a pool can
%       be sized before paying for a fit. Reject then Re-cluster is a
%       refinement loop (see onCluster), and reopening resumes the saved
%       curation rather than clustering afresh (see restoreSaved).
%
%       COST. A ripple pool is 5-20x an ED pool, and both the fit and the view
%       are bounded rather than left to scale with it: met.clust.nFit caps the
%       events evt_clust estimates on (the rest are projected and assigned),
%       and met.clust.nView caps the rows handed to guiTbl_xy per cluster.
%       Neither cap touches the MASK - every event is labelled and every event
%       is saved. The counts on the checkboxes are always the true ones.
%
%   INPUTS:
%       evt     - <struct> the loaded event struct (per-event fields aligned to
%                          .peakTime; needs .state for the state scope).
%       wv      - <mat|struct> [nEv x nSamp] waveforms, or a struct of them
%                          ({name -> [nEv x nSamp]}). Every field becomes a Y
%                          option in the view; the FIRST is what is clustered.
%       tstamps - <vec>   [1 x nSamp] window time base (s), from evt_maps.
%       cfg     - <struct> what this pipeline is:
%           .met      - <struct> .qa (the pool filter, see evt_gate) and .clust
%                                (evt_clust arguments, plus .scalar - a cellstr
%                                of per-event fields joining the shape - and
%                                .nView, the per-cluster row cap for the view).
%           .file     - <char>  the .mat holding the event struct.
%           .var      - <char>  its variable name there ('ripp' | 'ed'); also
%                               the evt_states file token.
%           .basepath - <char>  session directory (for evt_states).
%           .basename - <char>  file stem (for evt_states + the window title).
%           .lbl      - <char>  human label for plots ('Ripple' | 'ED').
%           .flgGui   - <log>   open the GUI (true) or filter headless.
%           .Visible  - <char>  'on' | 'off', for headless GUI tests.
%           .onSaved  - <fh>    called as onSaved(changed) after every save;
%                               [] = nothing to do.
%
%   OUTPUTS:
%       accepted - <log> [nEv x 1] the mask as of the call (headless: the
%                        filter; GUI: the seed, before the user touches it).
%       hFig     - <handle> the GUI figure ([] when headless).
%
%   DEPENDENCIES:
%       evt_gate, evt_clust, evt_states, evt_boutTimes, basepaths2vars,
%       backup_file; GUI: gui_layout, gui_labeledControl, gui_notify,
%       guiTbl_xy.
%
%   HISTORY:
%       260722 created by lifting the ED cluster-curation GUI (lfp/ed/ed_curate,
%              260721) into lfp/events so the ripple pipeline gets the same
%              tool instead of a second copy of it. Three things became
%              generic in the move: the metric knobs are BUILT FROM
%              met.qa.ranges rather than named in code, the waveforms are an
%              INPUT rather than a file this function knows how to find, and
%              the vigilance-state selection is recorded in the saved qa spec
%              (info.clustStates is read for old files, never written). Two
%              caps were added for ripple-scale pools - see COST above.

%% ========================================================================
%  SETUP
%  ========================================================================
if ~isstruct(wv), wv = struct('lfp', wv); end
met = cfg.met;
hFig = [];

%% ========================================================================
%  HEADLESS
%  ========================================================================
% The automatic gate: the FULL spec, states included, exactly as a batch run
% would apply it with no human. That mask is not an answer - it is the pool.
% A saved PARTITION is left alone. Headless makes no shape judgement, so it has
% nothing to say about the clusters - and nulling them would let one batch run
% destroy a session's manual curation with no warning and no way back.
if ~cfg.flgGui
    accepted = evt_gate(evt, met.qa);
    changed = saveCurated(cfg, accepted, [], [], met.qa);
    buildStates(cfg, evt, accepted);
    if ~isempty(cfg.onSaved), cfg.onSaved(changed); end
    return;
end

%% ========================================================================
%  GUI
%  ========================================================================
st = struct();
st.evt   = evt;
st.met   = met;
st.cfg   = cfg;
st.wv    = wv;
st.tst   = tstamps;
st.wvFld = fieldnames(wv);
st.state = plotState(evt);

hFig = uifigure('Name', sprintf('%s curation: %s', cfg.lbl, cfg.basename), ...
    'Position', [60 60 1600 850], 'Visible', cfg.Visible);
[~, gPlot, gCtrl, gActions] = gui_layout(hFig, 'CtrlWidth', 240);

st.knob = struct('h', {}, 'fld', {}, 'iBnd', {});
hFig.UserData = st;
st = buildKnobs(hFig, gCtrl);

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
st.cid    = nan(numel(evt.peakTime), 1);
st.nClust = 0;
st.acc    = true(numel(evt.peakTime), 1);
st.gate   = false(numel(evt.peakTime), 1);  % set by the first fit / restore
st.rej    = false(numel(evt.peakTime), 1);  % sticky shape rejections
st.pool   = false(numel(evt.peakTime), 1);
hFig.UserData = st;

buildStateChecks(hFig, prevStates(evt, met.qa));
if ~restoreSaved(hFig, evt)
    onCluster(hFig, true);
end
accepted = hFig.UserData.acc;

end     % EOF


% =========================================================================
%  GUI CALLBACKS
% =========================================================================
function onKnob(hFig)
% A metric knob moved: show the pool it implies. The clustering is NOT redone -
% a fit is not something to run on every keystroke, and the point of the knob is
% to choose a pool size before paying for it.
st = hFig.UserData;
gate = evt_gate(st.evt, poolSpec(st));
st.lblPool.Text = sprintf('gate %d | %d accepted  (press Re-cluster)', ...
    nnz(gate), nnz(gate & st.acc));

end     % onKnob


function onCluster(hFig, flgReset)
% Fit waveform clusters over the pool and adopt the result.
%
% Re-cluster fits what is left after the knobs, the SHAPE rejections and the
% state scope. Having thrown out a cluster of step artifacts, you do not want
% the groups spent describing events already rejected - you want them over what
% is left, which splits the survivors finer each round. The accepted SET does
% not move: every new cluster starts ticked, so a rejection lives in the input
% rather than in the tick pattern. Carrying the old ticks over would be wrong
% anyway - the partition is new, so index 3 no longer means what it did.
%
% THE TWO KINDS OF REJECTION BEHAVE DIFFERENTLY, and they have to.
%   SHAPE (a cluster untick) is STICKY, recorded in st.rej. Refitting strips
%     the labels of the events it drops, and an unlabelled event is
%     indistinguishable from one that was never judged - so without this record
%     a rejected cluster would walk straight back in on the next round.
%   STATE (a state untick) is NOT sticky. It is a scope, flipped back and forth
%     while working, so re-ticking a state makes its events eligible again and
%     the next Re-cluster gives them labels. Feeding the fit through the
%     accepted mask instead made it a one-way door: out of the fit meant no
%     label, no label meant not accepted, and not accepted meant it could never
%     re-enter the fit.
%
% FLGRESET clears the shape rejections and goes back to the whole gate; without
% it a mis-click would be unrecoverable short of reopening the session. It does
% not touch the knobs or the state scope, which are live either way.
st = hFig.UserData;
c  = st.met.clust;
[~, sel, states] = acceptMask(st);

% The ranges the pool - and therefore the labels, and therefore the mask - was
% actually built from. Save records THESE, not the live boxes: a knob typed into
% after a fit only previews (see onKnob), so writing its value would claim a
% filter the saved mask does not obey.
qaLive = buildSpec(st);
st.gate = evt_gate(st.evt, poolSpec(st));
st.qaFit = qaLive.ranges;
if flgReset
    st.rej = false(size(st.rej));
else
    st.rej = st.rej | (~isnan(st.cid) & ~ismember(st.cid, sel));
end
st.pool = st.gate & ~st.rej & ismember(st.state(:), states);
iPool = find(st.pool);

% 0 (or anything below 2) in the box asks for the rule, which is evt_clust's
% empty default - NOT met.clust.nClust, which is the shipped count the box
% already opened on
k = st.edK.Value;
if k < 2, k = []; end

st.cid = nan(numel(st.pool), 1);
st.nClust = 0;
if ~isempty(iPool)
    [cid, cInfo] = evt_clust(st.wv.(st.wvFld{1})(iPool, :), st.tst, ...
        'win', c.win, 'nPC', c.nPC, 'nClust', k, 'nFit', c.nFit, ...
        'scalar', scalarMat(st.evt, c.scalar, iPool), ...
        'detrend', c.detrend, 'norm', c.norm, 'wSize', c.wSize);
    st.cid(iPool) = cid;
    st.nClust = cInfo.nClust;

    % evt_clust refuses a pool too small to hold TYPES and returns no labels.
    % An unlabelled event cannot be accepted, so a refinement that worked -
    % one that narrowed the pool down to the handful of events worth keeping -
    % would throw every one of them away, with Reset the only way back. Below
    % that floor the honest answer is one group, not none.
    if st.nClust == 0
        st.cid(iPool) = 1;
        st.nClust = 1;
    end

    % An event whose waveform is all NaN (a peak within the map window of the
    % recording edge) comes back unlabelled and can never be accepted. Record
    % it as rejected: left out, it would re-enter the pool on every Re-cluster
    % and keep the "unsorted - press Re-cluster" hint alive against a button
    % that provably cannot clear it.
    st.rej = st.rej | (st.pool & isnan(st.cid));
end
hFig.UserData = st;

adoptClust(hFig, 1 : st.nClust, ...
    sprintf('gate %d | clustered %d -> %d groups', nnz(st.gate), ...
    numel(iPool), st.nClust));

end     % onCluster


function ok = restoreSaved(hFig, evt)
% Put back the partition and the choices a previous session saved, instead of
% clustering afresh. The saved labels ARE the partition that was judged: a
% re-fit would cost a fit and hand back a different one, leaving the ticks
% pointing at groups nobody looked at. Anything that no longer lines up with
% the event list is refused, so a re-detection falls through to a fresh fit
% rather than drawing labels that belong to other events.
ok = isfield(evt, 'clustId') && numel(evt.clustId) == numel(evt.peakTime) ...
    && ~all(isnan(evt.clustId));
if ~ok, return; end

st = hFig.UserData;
st.cid    = evt.clustId(:);
st.pool   = ~isnan(st.cid);
st.nClust = max(st.cid);

% An EMPTY saved selection is a real answer - "every cluster was rejected" -
% and must come back as such. Only the ABSENCE of the field means no choice was
% recorded; the headless gate cannot produce one, because it leaves the
% partition untouched and a file with no partition never reaches here.
sel = 1 : st.nClust;
if isfield(evt.info, 'clustSel')
    sel = evt.info.clustSel;
end

setKnobs(st, evt.info.qa);
st.gate = evt_gate(st.evt, poolSpec(st));
qaLive = buildSpec(st);
st.qaFit = qaLive.ranges;

% Rebuild the shape rejections the saved partition implies: an event the gate
% passed and the saved state scope included, yet which carries no label, was
% dropped by a cluster untick in some earlier round. Without this a reopen
% would quietly undo the refinement and offer those events again.
scope = prevStates(evt, st.met.qa);
if isempty(scope), scope = st.stateCats; end
st.rej = st.gate & isnan(st.cid) & ismember(st.state(:), scope);
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

% one checkbox per cluster, labelled with its TRUE size - the view may show a
% capped sample of a cluster, but the count a decision is made on never is
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
[accepted, ~, states] = acceptMask(st);
st.acc = accepted;
hFig.UserData = st;

% An event needs a label to be accepted, so re-ticking a state cannot bring
% its events back on its own - they left the last fit and have none. Say so,
% rather than letting the tick look like it did nothing.
txt = sprintf('accepted: %d of %d', nnz(accepted), numel(accepted));
nPend = nnz(st.gate & ~st.rej & isnan(st.cid) & ...
    ismember(st.state(:), states));
if nPend > 0
    txt = sprintf('%s  |  %d unsorted - press Re-cluster', txt, nPend);
end

switch st.ddShow.Value
    case 'accepted', iRow = find(accepted);
    case 'removed',  iRow = find(~accepted);
    otherwise,       iRow = (1 : numel(accepted))';
end
nAll = numel(iRow);
tbl = viewTable(st, accepted, iRow);
if height(tbl) < nAll
    txt = sprintf('%s  |  view sampled %d of %d rows', txt, height(tbl), nAll);
end
st.lblKeep.Text = txt;

ud = st.hPanel.UserData;
if isstruct(ud) && isfield(ud, 'setDataFcn')
    ud.setDataFcn(tbl);
elseif ~isempty(iRow)
    % first call. An empty table has no categories to build the filter panels
    % from, so the widget waits until there is something to draw. The map is
    % cut wider than this on purpose - the flanks are where one shape separates
    % from another, and they are kept in the file - but the view opens on the
    % clustering window, which is where the decisions are made. Zoom out to see
    % the rest.
    guiTbl_xy(st.tst * 1000, tbl, 'Parent', st.hPanel, ...
        'yVar', st.wvFld{1}, 'tileVar', 'cluster', 'grpVar', 'state', ...
        'xLbl', 'time (ms)', 'xLim', st.met.clust.win * 1000);
end

end     % refresh


function doSave(hFig)
% Persist the mask plus everything needed to resume it.
st = hFig.UserData;
[accepted, sel] = acceptMask(st);

% the state scope is live (it feeds the mask directly), the metric ranges are
% the ones the last fit used (a knob typed into since then only previewed)
qa = buildSpec(st);
if isfield(st, 'qaFit'), qa.ranges = st.qaFit; end

changed = saveCurated(st.cfg, accepted, st.cid, sel, qa);
buildStates(st.cfg, st.evt, accepted);
msg = sprintf('Saved: %d events from %d clusters (+ %sStates)', ...
    nnz(accepted), numel(sel), st.cfg.var);
if ~isempty(st.cfg.onSaved)
    extra = st.cfg.onSaved(changed);
    if ~isempty(extra), msg = sprintf('%s. %s', msg, extra); end
end
gui_notify(hFig, msg, 'success');

end     % doSave


% =========================================================================
%  VIEW
% =========================================================================
function tbl = viewTable(st, accepted, iRow)
% The rows IROW asks for, with every waveform and the three things worth
% pivoting on. Events the filter dropped keep the cluster label 'out' rather
% than being hidden, so a per-state view still shows what detection proposed.
%
% NVIEW caps the rows per cluster. guiTbl_xy draws one Line object per row in
% 'Traces' mode and sets two properties on each in a loop, so a 30k-event
% ripple pool locks the figure; and no eye reads 3000 overlaid traces anyway.
% The cap is EVENLY SPACED rather than random - deterministic, and spread over
% the recording rather than clumped. It moves rows only: the mask, the counts
% on the checkboxes and what Save writes are all computed on every event.
iRow = capRows(iRow, st.cid, st.met.clust.nView);

idx = st.cid;
idx(isnan(idx)) = 0;
clust = categorical(idx, 0 : st.nClust, ...
    [{'out'}, cellstr(string(1 : st.nClust))]);

tbl = table(clust(iRow), st.state(iRow), ...
    categorical(accepted(iRow), [false true], {'removed', 'kept'}), ...
    'VariableNames', {'cluster', 'state', 'status'});
for iFld = 1 : numel(st.wvFld)
    tbl.(st.wvFld{iFld}) = st.wv.(st.wvFld{iFld})(iRow, :);
end

end     % viewTable


function iRow = capRows(iRow, cid, nView)
% At most NVIEW rows per cluster, evenly spaced through each.
if isempty(nView) || ~isfinite(nView) || numel(iRow) <= nView, return; end
grp = cid(iRow);
grp(isnan(grp)) = 0;
keep = false(numel(iRow), 1);
uGrp = unique(grp)';
for iGrp = 1 : numel(uGrp)
    iG = find(grp == uGrp(iGrp));
    if numel(iG) > nView
        iG = iG(round(linspace(1, numel(iG), nView)));
    end
    keep(iG) = true;
end
iRow = iRow(keep);

end     % capRows


function s = plotState(evt)
% evt.state as a tiling variable: <undefined> is promoted to its own 'unscored'
% level. A categorical comparison never matches <undefined>, so without this the
% unscored events would get no tile and vanish from the view without a word. A
% session that was never scored has no .state at all.
if ~isfield(evt, 'state') || isempty(evt.state)
    s = categorical(repmat({'unscored'}, numel(evt.peakTime), 1));
    return;
end
s = removecats(evt.state(:));
if any(isundefined(s))
    s = addcats(s, {'unscored'});
    s(isundefined(s)) = 'unscored';
end

end     % plotState


% =========================================================================
%  CONTROLS
% =========================================================================
function [accepted, sel, states] = acceptMask(st)
% An event is real if its CLUSTER is ticked AND its STATE is ticked. Clusters
% say which shape; states are there to drop a whole stretch of the recording -
% movement artifact in WAKE, say - without touching the shape decision. Both
% start ticked, so the mask begins as the whole pool and curation removes from
% it. SEL and STATES come back with the mask because Save records the choice as
% well as its result; deriving them apart is how the two drift.
sel    = find(ticked(st.chk));
states = st.stateCats(ticked(st.chkState));
states = states(:)';
accepted = ismember(st.cid, sel) & ismember(st.state(:), states);

end     % acceptMask


function st = buildKnobs(hFig, gCtrl)
% One edit box per FINITE bound in met.qa.ranges, so the pool filter is
% described in one place (the methods file) and this function never learns a
% metric's name. An open bound gets no box - there is nothing to tune - and a
% metric the GUI shows no box for still rides through buildSpec untouched.
st = hFig.UserData;
qa = st.met.qa;
if ~isfield(qa, 'ranges') || isempty(qa.ranges)
    hFig.UserData = st;
    return;
end
flds = fieldnames(qa.ranges);
op = {'>=', '<='};
for iFld = 1 : numel(flds)
    bnd = qa.ranges.(flds{iFld});
    for iBnd = 1 : 2
        if ~isfinite(bnd(iBnd)), continue; end
        k = numel(st.knob) + 1;
        st.knob(k).h = gui_labeledControl(gCtrl, 'editnum', ...
            sprintf('%s %s', flds{iFld}, op{iBnd}), 'Value', bnd(iBnd), ...
            'ValueChangedFcn', @(~,~) onKnob(hFig));
        st.knob(k).fld  = flds{iFld};
        st.knob(k).iBnd = iBnd;
    end
end
hFig.UserData = st;

end     % buildKnobs


function setKnobs(st, qa)
% Put a saved spec back into the boxes. A bound the spec does not carry leaves
% its box alone, so a partial spec is not an error.
if ~isfield(qa, 'ranges'), return; end
for iK = 1 : numel(st.knob)
    if ~isfield(qa.ranges, st.knob(iK).fld), continue; end
    v = qa.ranges.(st.knob(iK).fld)(st.knob(iK).iBnd);
    if isfinite(v), st.knob(iK).h.Value = v; end
end

end     % setKnobs


function buildStateChecks(hFig, keepCats)
% One checkbox per vigilance state. Built ONCE - the states of a session do not
% change - so these ticks survive every re-cluster. KEEPCATS restores a saved
% selection; empty means none was saved, and everything starts ticked.
st = hFig.UserData;
cats = categories(removecats(st.state));

g = uigridlayout(st.gState, [numel(cats), 1], 'Padding', 2, ...
    'RowHeight', repmat({'fit'}, 1, numel(cats)), 'RowSpacing', 1);
st.chkState = gobjects(numel(cats), 1);
for iCat = 1 : numel(cats)
    val = isempty(keepCats) || ismember(cats{iCat}, keepCats);
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
% The full evt_gate spec the controls currently describe: the metric knobs plus
% the state scope. Started from met.qa so a range with no box (both bounds open,
% or a metric added to the methods file later) is carried through rather than
% silently dropped.
qa = st.met.qa;
for iK = 1 : numel(st.knob)
    qa.ranges.(st.knob(iK).fld)(st.knob(iK).iBnd) = st.knob(iK).h.Value;
end

sel = st.stateCats(ticked(st.chkState));
qa.unscored = ismember('unscored', sel);
qa.states = setdiff(sel, {'unscored'}, 'stable');
qa.states = qa.states(:)';
% evt_gate reads an empty state list as "any state", so "no real state is
% ticked" has to be said with a label that matches none of them. Unconditional:
% on a session with no scoring at all, the only box is 'unscored', and
% unticking it still has to mean nothing is kept.
if isempty(qa.states)
    qa.states = {'<none>'};
end

end     % buildSpec


function qa = poolSpec(st)
% The spec that defines the POOL: the metric knobs only. State is a scope
% applied alongside the cluster choice (see onCluster), never inside the gate -
% putting it here made unticking a state a one-way door, because the events left
% the fit, lost their labels, and an unlabelled event can never be accepted.
qa = buildSpec(st);
qa.states = [];

end     % poolSpec


function cats = prevStates(evt, qaDflt)
% The state selection a previous session saved, or {} when there is none.
%
% An EMPTY saved selection counts as none. It would otherwise mean "every state
% was rejected", which restores a GUI that accepts nothing and offers no clue
% why - and that is exactly what a file written by the headless gate looks like.
% The case it gives up on (deliberately saving a curation that keeps nothing) is
% not worth the one it breaks.
%
% QADFLT is the shipped spec, used only when nothing was saved: it is how a
% methods file says which states an uncurated session counts.
cats = {};
if ~isfield(evt, 'info'), return; end
if isfield(evt.info, 'qa')
    cats = stateNames(evt.info.qa);
end
% pre-260722 ED files recorded the scope in its own field
if isempty(cats) && isfield(evt.info, 'clustStates') ...
        && iscell(evt.info.clustStates)
    cats = evt.info.clustStates(:)';
end
if isempty(cats), cats = stateNames(qaDflt); end

end     % prevStates


function names = stateNames(qa)
% The state labels a spec keeps: [] / absent is "any state", which is not a
% selection. Numeric AccuSleep indices resolve through the shared evt_gate rule.
names = {};
if ~isstruct(qa), return; end
if isfield(qa, 'states') && ~isempty(qa.states)
    if isnumeric(qa.states)
        cfg = as_loadConfig([]);
        idx = qa.states(qa.states >= 1 & qa.states <= numel(cfg.names));
        names = cfg.names(idx);
    else
        names = cellstr(qa.states);
    end
    names = names(:)';
end
% appended outside the block: a spec that keeps ONLY the unlabelled events has
% an empty real-state list, and its selection is still a selection
if isfield(qa, 'unscored') && qa.unscored, names = [names, {'unscored'}]; end

end     % stateNames


function S = scalarMat(evt, flds, iPool)
% The per-event measures that join the shape components, as one matrix. They
% are shape descriptors already computed, so withholding them from the
% clustering would only throw information away. A field the struct does not
% carry is skipped rather than fatal, so a methods file can name a measure a
% given session lacks.
S = [];
if isempty(flds), return; end
for iFld = 1 : numel(flds)
    if ~isfield(evt, flds{iFld}), continue; end
    v = evt.(flds{iFld});
    if ~isnumeric(v) || numel(v) ~= numel(evt.peakTime), continue; end
    S = [S, double(v(iPool))]; %#ok<AGROW>
end

end     % scalarMat


% =========================================================================
%  PERSISTENCE
% =========================================================================
function changed = saveCurated(cfg, accepted, clustId, sel, qa)
% Back up, then overwrite the mask + everything needed to resume: the labels,
% which clusters were accepted, and the spec (knobs + state scope) behind them.
% Loads and rewrites the WHOLE file, so anything stored alongside survives.
%
% CHANGED says whether the mask on disk actually moved - the caller's hook uses
% it to decide whether downstream products are now stale.
backup_file(cfg.file);
S = load(cfg.file);
accepted = logical(accepted(:));
old = S.(cfg.var);
changed = ~isfield(old, 'accepted') || ...
    ~isequal(logical(old.accepted(:)), accepted);

S.(cfg.var).accepted = accepted;
if ~isfield(S.(cfg.var), 'info') || ~isstruct(S.(cfg.var).info)
    S.(cfg.var).info = struct();
end
S.(cfg.var).info.qa = qa;
% an empty CLUSTID means "this caller made no shape judgement" (the headless
% gate), which is not the same as "every cluster was rejected" - so the stored
% partition and the choice over it are left exactly as they were
if ~isempty(clustId)
    S.(cfg.var).clustId       = clustId(:);
    S.(cfg.var).info.clustSel = sel;
end
% accepted-aligned, so it belongs to a mask that no longer exists; the analyze
% stage rebuilds it
if isfield(S.(cfg.var), 'spks')
    S.(cfg.var) = rmfield(S.(cfg.var), 'spks');
end
save(cfg.file, '-struct', 'S', '-v7.3');

end     % saveCurated


function buildStates(cfg, evt, accepted)
% Rebuild + save the per-bout rate/density table for the current mask (cheap;
% no signal, no spikes). Skips silently when sleep states are unavailable.
win = [0 Inf];
if isfield(evt, 'info') && isfield(evt.info, 'win'), win = evt.info.win; end
w0 = win(1);
if ~isfinite(w0), w0 = 0; end

% basepaths2vars warns and skips a var it cannot find, so a session without
% scoring simply arrives without the field
v = basepaths2vars('basepaths', {cfg.basepath}, ...
    'vars', {'session', 'sleep_states'});
if isinf(win(2)) && isfield(v, 'session') && ~isempty(v.session)
    win(2) = v.session.extracellular.nSamples / v.session.extracellular.srLfp;
end

boutTimes = evt_boutTimes(v, win, win(2) - win(1));
if isempty(boutTimes)
    return;                     % no scoring -> no rate/density table
end
evt_states(evt.times - w0, evt.peakTime - w0, boutTimes, ...
    'accepted', logical(accepted(:)), 'basepath', cfg.basepath, ...
    'basename', cfg.basename, 'flgSave', true, 'flgPlot', false, ...
    'name', cfg.var, 'lbl', cfg.lbl);

end     % buildStates
