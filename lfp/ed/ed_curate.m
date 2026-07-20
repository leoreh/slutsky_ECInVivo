function [ed, hFig] = ed_curate(basepath, varargin)
% ED_CURATE Curate discharges by waveform TYPE, not one event at a time.
%
%   [ed, hFig] = ED_CURATE(basepath, varargin)
%
%   SUMMARY:
%       Stage 2 of the ED pipeline. Loads <basename>.ed.mat, applies met.qa as
%       a noise filter, groups what survives into waveform clusters (ed_clust),
%       and lets you accept whole clusters. Saving writes .accepted, .clustId
%       and the choice into ed.info.
%
%       Why clusters. A 24 h recording proposes thousands of candidates and
%       holds a few dozen discharges. Judging that one event at a time is a day
%       of work; judging it by a MEAN waveform is worse than useless, because
%       an average over discharges, sharp waves and step artifacts is a curve
%       that is none of them - which is exactly how the discharges got buried
%       when this pipeline was first calibrated. Split by shape first and every
%       tile shows a real waveform, with a median and an IQR band rather than a
%       mean, so one outlier cannot set the picture.
%
%       Nothing here knows what a discharge looks like. The clustering is blind
%       and you name the clusters, so a mouse whose discharges differ from the
%       raMCU3/4/5 shape still gets them in a cluster of their own - which is
%       the whole point, since polarity and sharpness are layer-dependent.
%
%       The MUA row is CONTEXT, never a criterion. A discharge is followed by a
%       prolonged drop in population firing (0.18-0.62 of baseline in the
%       curated mice, versus ~1.0 for artifacts and ripples - see
%       dev/ed_pipeline_rebuild.md), so it is a useful second opinion on a
%       cluster you are unsure about. It is drawn when spikes exist and
%       silently skipped when they do not; nothing depends on it.
%
%       Headless (flgGui = false) applies only the noise filter, for a batch
%       run that has no human. That mask is NOT an answer - it is the pool the
%       GUI sorts.
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
%       evt_files, evt_gate, evt_states, evt_boutTimes, evt_spkPrep,
%       ed_methods, ed_clust, basepaths2vars, backup_file; GUI: gui_layout,
%       gui_labeledControl, gui_notify.
%
%   HISTORY:
%       260720 created as the ED twin of ripp_curate (threshold knobs over a
%              kept-vs-removed mean waveform).
%       260721 rebuilt around waveform clustering. The knobs are gone: they
%              asked the human to express "which shape is a discharge" as three
%              numbers, which is the wrong question put the wrong way round.

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

pool = evt_gate(ed, met.qa);            % the noise filter

%% ========================================================================
%  HEADLESS
%  ========================================================================
if ~flgGui
    ed.accepted = pool;
    saveCurated(files.evt, pool, nan(numel(pool), 1), [], met.qa);
    buildStates(basepath, basename, ed, pool);
    if verbose
        fprintf('[ED_CURATE] %s : %d / %d pass the noise filter\n', ...
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
st.pool     = pool;
[st.wv, st.tst] = loadMaps(files.maps, ed);
st.mua      = loadMua(basepath, ed);

hFig = uifigure('Name', ['ED curation: ' basename], ...
    'Position', [60 60 1600 800], 'Visible', vis);
[~, gPlot, gCtrl, gActions] = gui_layout(hFig, 'CtrlWidth', 230);

st.lblPool = gui_labeledControl(gCtrl, 'label', ...
    sprintf('%d of %d pass the filter', nnz(pool), numel(pool)));
st.edK = gui_labeledControl(gCtrl, 'editnum', 'clusters', ...
    'Value', met.clust.nClust);
gui_labeledControl(gCtrl, 'button', '', 'Text', 'Re-cluster', ...
    'ButtonPushedFcn', @(~,~) onCluster(hFig));
st.gClust = gui_labeledControl(gCtrl, 'panel', 'accept as discharge', ...
    'RowHeight', '1x');
st.lblKeep = gui_labeledControl(gCtrl, 'label', '');

gui_labeledControl(gActions, 'button', '', 'Text', 'Save', ...
    'ButtonPushedFcn', @(~,~) doSave(hFig));

st.gPlot = uigridlayout(gPlot, [1 1], 'Padding', 0);
st.chk = gobjects(0);
hFig.UserData = st;

onCluster(hFig);

end     % EOF


% =========================================================================
%  GUI CALLBACKS
% =========================================================================
function onCluster(hFig)
% Cluster the pool, then rebuild the checkbox list and the tiles.
st = hFig.UserData;
c = st.met.clust;

k = st.edK.Value;
if k < 2, k = c.nClust; end

iPool = find(st.pool);
[cid, cInfo] = ed_clust(st.wv(iPool, :), st.tst, 'win', c.win, ...
    'nPC', c.nPC, 'nClust', k, 'scalar', scalarFeat(st.ed, iPool));

st.iPool  = iPool;
st.cid    = cid;
st.nClust = cInfo.nClust;
hFig.UserData = st;

% a pool too small to hold types at all - a control mouse, usually. Say so
% rather than drawing an empty grid; Save still writes an all-false mask.
if st.nClust == 0
    delete(st.gClust.Children);
    delete(st.gPlot.Children);
    st.chk = gobjects(0);
    st.ax  = gobjects(0);
    hFig.UserData = st;
    st.lblKeep.Text = sprintf('%d events: too few to cluster', numel(iPool));
    return;
end

st.edK.Value = cInfo.nClust;
hFig.UserData = st;

buildChecks(hFig);
drawTiles(hFig);

end     % onCluster


function buildChecks(hFig)
% One checkbox per cluster, labelled with its size.
st = hFig.UserData;
delete(st.gClust.Children);
g = uigridlayout(st.gClust, [st.nClust + 1, 1], 'Padding', 2, ...
    'RowHeight', repmat({'fit'}, 1, st.nClust + 1), 'RowSpacing', 1);

st.chk = gobjects(st.nClust, 1);
for iK = 1 : st.nClust
    st.chk(iK) = uicheckbox(g, 'Value', false, 'Text', ...
        sprintf('%d   (n = %d)', iK, nnz(st.cid == iK)), ...
        'ValueChangedFcn', @(~,~) onPick(hFig));
end
hFig.UserData = st;

end     % buildChecks


function onPick(hFig)
% A cluster was ticked: recolour its tile and update the count.
st = hFig.UserData;
sel = arrayfun(@(h) h.Value, st.chk);
for iK = 1 : st.nClust
    if isgraphics(st.ax(iK))
        st.ax(iK).Color = tileColor(sel(iK));
    end
end
st.lblKeep.Text = sprintf('accepted: %d events', ...
    nnz(ismember(st.cid, find(sel))));

end     % onPick


function doSave(hFig)
% Persist the mask, the labels and which clusters were chosen.
st = hFig.UserData;
sel = find(arrayfun(@(h) h.Value, st.chk));
sel = sel(:)';

accepted = false(numel(st.pool), 1);
accepted(st.iPool(ismember(st.cid, sel))) = true;
clustId = nan(numel(st.pool), 1);
clustId(st.iPool) = st.cid;

saveCurated(st.files.evt, accepted, clustId, sel, st.met.qa);
buildStates(st.basepath, st.basename, st.ed, accepted);
gui_notify(hFig, sprintf('Saved: %d accepted from %d clusters (+ edStates)', ...
    nnz(accepted), numel(sel)), 'success');

end     % doSave


% =========================================================================
%  DRAWING
% =========================================================================
function drawTiles(hFig)
% One column per cluster: median waveform + IQR band on top, peri-event MUA
% below when spikes exist.
st = hFig.UserData;
delete(st.gPlot.Children);

hasMua = st.mua.ok;
nRow = 1 + hasMua;
g = uigridlayout(st.gPlot, [nRow, st.nClust], 'Padding', 4, ...
    'RowHeight', repmat({'1x'}, 1, nRow));

st.ax = gobjects(st.nClust, 1);
for iK = 1 : st.nClust
    inK = st.cid == iK;
    ax = uiaxes(g); hold(ax, 'on');
    st.ax(iK) = ax;

    W = double(st.wv(st.iPool(inK), :));
    q = prctile(W, [25 50 75], 1);
    fill(ax, [st.tst, fliplr(st.tst)] * 1000, [q(1, :), fliplr(q(3, :))], ...
        [0.3 0.3 0.3], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(ax, st.tst * 1000, q(2, :), 'k', 'LineWidth', 1.5);
    xline(ax, 0, ':');
    title(ax, sprintf('%d  (n = %d)', iK, nnz(inK)), 'FontSize', 9);
    ax.Color = tileColor(false);
    if iK == 1, ylabel(ax, 'LFP (uV)'); end
    xlim(ax, [-100 100]);

    if hasMua
        axM = uiaxes(g); hold(axM, 'on');
        r = muaRate(st.mua, st.ed.peakTime(st.iPool(inK)));
        plot(axM, st.mua.ctrs * 1000, r, 'Color', [0 0.4 0.8], ...
            'LineWidth', 1.2);
        yline(axM, 1, ':'); xline(axM, 0, ':');
        xlim(axM, [-400 400]); ylim(axM, [0 2]);
        if iK == 1, ylabel(axM, 'MUA / baseline'); end
        xlabel(axM, 'time (ms)');
    else
        xlabel(ax, 'time (ms)');
    end
end
hFig.UserData = st;
onPick(hFig);

end     % drawTiles


function c = tileColor(isSel)
if isSel, c = [0.87 0.95 0.87]; else, c = [1 1 1]; end
end


% =========================================================================
%  DATA
% =========================================================================
function [wv, tst] = loadMaps(file, ed)
% Per-event waveforms behind the clustering and the tiles.
if ~isfile(file)
    error('ed_curate:noMaps', ...
        'no edMaps file; re-run detection with flgSave.');
end
S = load(file, 'edMaps');
if size(S.edMaps.lfp, 1) ~= numel(ed.peakTime)
    error('ed_curate:staleMaps', ...
        'edMaps does not match the event list; re-run detection.');
end
wv  = S.edMaps.lfp;
tst = S.edMaps.tstamps;

end     % loadMaps


function mua = loadMua(basepath, ed)
% Pooled multi-unit times + the psth bins. Empty when there are no spikes,
% which is not an error - the MUA row is context, not a criterion.
BIN = 0.020; HALF = 0.5;
mua = struct('ok', false, 't', [], 'ctrs', [], 'edges', []);

v = basepaths2vars('basepaths', {basepath}, ...
    'vars', {'spikes', 'spktimes', 'session'}, 'flgPrnt', false);
if ~isfield(v, 'session') || ~isstruct(v.session) ...
        || ~isfield(v.session, 'extracellular')
    return;
end
fsSpk = v.session.extracellular.sr;
[~, muTimes] = evt_spkPrep(v, [0 Inf], ed.info.sigDur, fsSpk);
if isempty(muTimes) || isempty(muTimes{1}), return; end

mua.t     = muTimes{1};
mua.edges = -HALF : BIN : HALF;
mua.ctrs  = mua.edges(1 : end - 1) + BIN / 2;
mua.ok    = true;

end     % loadMua


function r = muaRate(mua, peakTime)
% Peri-event multi-unit rate, over its own baseline. Binary search per event
% keeps this instant even with a million spikes.
BASE = [-0.5 -0.2];
cnt = zeros(1, numel(mua.ctrs));
lo = discretize(peakTime + mua.edges(1), [-inf; mua.t(:); inf]);
hi = discretize(peakTime + mua.edges(end), [-inf; mua.t(:); inf]);
for iE = 1 : numel(peakTime)
    if isnan(lo(iE)) || hi(iE) <= lo(iE), continue, end
    d = mua.t(lo(iE) : min(hi(iE) - 1, numel(mua.t))) - peakTime(iE);
    cnt = cnt + histcounts(d, mua.edges);
end
r = cnt / max(1, numel(peakTime));
bl = mean(r(mua.ctrs >= BASE(1) & mua.ctrs <= BASE(2)));
if bl > 0, r = r / bl; else, r = nan(size(r)); end

end     % muaRate


function s = scalarFeat(ed, idx)
% The per-event shape measures that join the waveform components. They are
% shape descriptors already computed, so withholding them from the clustering
% would only throw information away.
s = [ed.fastZ(idx), ed.isoZ(idx), ed.posZ(idx), ed.amp(idx), ed.dur(idx)];

end     % scalarFeat


% =========================================================================
%  PERSISTENCE
% =========================================================================
function saveCurated(file, accepted, clustId, sel, qa)
% Back up, then overwrite the mask + the labels in the saved struct.
backup_file(file);
S = load(file);
S.ed.accepted    = logical(accepted(:));
S.ed.clustId     = clustId(:);
S.ed.info.qa     = qa;
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
