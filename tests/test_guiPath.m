function test_guiPath(sessionPath)
% TEST_GUIPATH Smoke test for the guiPath curation viewer.
%
%   test_guiPath() runs the data-free unit checks only.
%   test_guiPath(sessionPath) also opens the viewer on a real session and
%   checks the end-to-end path (int16 stack, tick lanes, amplitude, region).
%
%   The unit checks call guiPath_draw and guiPath_shape directly, which the
%   draw-layer extraction and the data / view split are what make possible;
%   the integration checks need a session with a <basename>.ripp.mat and a
%   <basename>.lfp, so they are skipped when no path is given or found.
%
%   HISTORY:
%       260716 - created alongside the traces / amplitude / region work.
%       260719 - presets addressed by file token (ripp, sleep_states).

if nargin < 1, sessionPath = ''; end
fprintf('== test_guiPath ==\n');

%% ------------------------------------------------------------------ unit
% 1. eventTicks draws the curated set's ACCEPTED events as one full-height line
% in the set's colour; rejected are hidden (rejecting removes the tick). Events
% live on the panel input; data.curate names the edited set (its live mask).
inp = uinput('eventTicks');
inp.name = 'ripp'; inp.clr = [0 0.2 0.8];
inp.data = struct('peakTime', [10; 20; 30; 40; 50]);
data = struct('curate', 'ripp', 'accepted', logical([1; 1; 0; 1; 0]));
f = figure('Visible', 'off'); ax = axes(f); hold(ax, 'on');
guiPath_draw(ax, inp, data, 0, 60, 1);
L = findobj(ax, 'Type', 'line');
assert(isscalar(L), 'eventTicks: expected one line (accepted only), got %d', numel(L));
y = L.YData(~isnan(L.YData));
assert(abs(min(y)) < 1e-9 && abs(max(y) - 1) < 1e-9, 'eventTicks: ticks not full height');
assert(isequal(L.Color, [0 0.2 0.8]), 'eventTicks: ticks not the set colour');
assert(numel(unique(L.XData(~isnan(L.XData)))) == 3, 'eventTicks: expected 3 accepted ticks');
close(f);
fprintf('  ok  eventTicks accepted-only, full-height, set colour\n');

% 2a. trace amplitude: yAdjust tightens the y-limits about centre.
sig = sin(2*pi*8*(0:1/1250:2)');           % 8 Hz, +/-1
tr = uinput('trace'); tr.data = sig; tr.fs = 1250; tr.ylim = [-1 1];
f = figure('Visible', 'off'); ax = axes(f); hold(ax, 'on');
tr.yAdjust = 1;    guiPath_draw(ax, tr, struct(), 0, 2, 1); base = diff(ax.YLim);
cla(ax); hold(ax, 'on');
tr.yAdjust = 2;    guiPath_draw(ax, tr, struct(), 0, 2, 1); tight = diff(ax.YLim);
assert(abs(tight - base/2) < 1e-9, 'trace amp: yAdjust=2 should halve the y-span');
close(f);
fprintf('  ok  trace amplitude (y-limits)\n');

% 2b. stack amplitude: yAdjust is a GAIN - one channel's deflection doubles at
% gain 2, while the y-limits stay put (so channels do not clip a shrinking axis).
st = uinput('traces'); st.data = sig; st.fs = 1250;   % single-channel stack
st.chInfo = struct('spacing', 4, 'base', 0, 'labels', 1);
f = figure('Visible', 'off'); ax = axes(f); hold(ax, 'on');
st.yAdjust = 1; guiPath_draw(ax, st, struct(), 0, 2, 1);
yl1 = ax.YLim; a1 = localSpan(ax);
cla(ax); hold(ax, 'on');
st.yAdjust = 2; guiPath_draw(ax, st, struct(), 0, 2, 1);
yl2 = ax.YLim; a2 = localSpan(ax);
assert(isequal(yl1, yl2), 'stack amp: y-limits must stay fixed under gain');
assert(abs(a2 - 2*a1) < 1e-6, 'stack amp: gain=2 should double the drawn deflection');
close(f);
fprintf('  ok  stack amplitude (gain, fixed limits)\n');

% 3. ylim resolution (via guiPath_shape): a scalar percentile clips tighter
% than prc.
rng(0); s = randn(2e5, 1); s(1000) = 500;   % one artifact
a = guiPath_shape(struct('type', 'trace', 'data', s, 'fs', 1250, 'ylim', 'prc', 'labels', []));
b = guiPath_shape(struct('type', 'trace', 'data', s, 'fs', 1250, 'ylim', 5,     'labels', []));
assert(diff(b.ylim) < diff(a.ylim), 'resolveYlim: scalar 5 should clip tighter than 0.1');
fprintf('  ok  scalar percentile ylim\n');

%% ----------------------------------------------------------- integration
if isempty(sessionPath) || ~isfolder(sessionPath)
    fprintf('  skip integration (no session: %s)\n', sessionPath);
    fprintf('== PASSED (unit only) ==\n');
    return;
end

hFig = guiPath(sessionPath, 'preset', 'ripp', 'Visible', 'off');
d = hFig.UserData;
assert(d.nEvents > 0, 'integration: no events loaded');

% stack is int16 in memory, labelled with real channels, drawn one line / ch
ix = find(strcmp({d.inputs.type}, 'traces'), 1);
assert(~isempty(ix), 'integration: no traces panel');
q = d.inputs(ix);
assert(isa(q.data, 'int16'), 'integration: stack must stay int16');
assert(~isempty(q.chInfo) && ~isempty(q.chInfo.labels), 'integration: stack has no channel labels');
fprintf('  ok  int16 stack, channels %s\n', mat2str(q.chInfo.labels));

% Tend is the session length, not nCh x it
assert(d.Tend_s < 48*3600, 'integration: Tend inflated (computeTend counted the matrix)');

% a row-vector signal (emg_rms) stays a full time series, not averaged to one
% scalar (the bin-averaging must fire only for a bin: matrix)
je = find(strcmp({d.inputs.name}, 'emgRms'), 1);
assert(~isempty(je) && numel(d.inputs(je).data) > 1000, ...
    'integration: emg_rms collapsed to a scalar (trace-averaging misfired)');
fprintf('  ok  emg_rms is a full trace (%d samples)\n', numel(d.inputs(je).data));

% event ticks: one full-height blue line, no y-ticks (the Ripples set is 'ripp')
jt = find(strcmp({d.wideP.source}, 'ripp'), 1);
axT = d.wideP(jt).ax;
Lt = findobj(axT, 'Type', 'line');
assert(isscalar(Lt) && (max(Lt.YData(~isnan(Lt.YData))) - min(Lt.YData(~isnan(Lt.YData)))) > 0.99, ...
    'integration: event ticks not a single full-height line');
assert(isequal(Lt.Color, [0 0.2 0.8]), 'integration: ticks not blue');
assert(isempty(axT.YTick), 'integration: event strip should have no y-ticks');
fprintf('  ok  ticks blue, full-height, no y-ticks\n');

% the event strip is the bottom Top panel, so it carries the Time axis (hours);
% panels above it show no x-tick numbers
assert(strcmp(d.wideP(end).source, 'ripp'), ...
    'integration: the event strip should be the last Top panel');
assert(~isempty(d.wideP(end).ax.XTickLabel), ...
    'integration: the bottom Top panel (event strip) should carry the Time axis');
assert(isempty(d.wideP(1).ax.XTickLabel), 'integration: a non-bottom panel should show no x numbers');
fprintf('  ok  Time axis on the bottom Top panel (event strip)\n');

% amplitude key round-trip on the active panel
d.activeSrc = q.name; hFig.UserData = d;
evt.Modifier = {'shift'}; evt.Key = 'equal';
hFig.WindowKeyPressFcn(hFig, evt);
assert(hFig.UserData.inputs(ix).yAdjust > 1, 'integration: shift+= should raise yAdjust');
evt.Key = '0'; hFig.WindowKeyPressFcn(hFig, evt);
assert(hFig.UserData.inputs(ix).yAdjust == 1, 'integration: shift+0 should reset yAdjust');
fprintf('  ok  amplitude key round-trip (bottom stack)\n');

% amplitude also reaches Top panels (a trace and the spectrogram), redrawing
% them in place without error
for topType = {'trace', 'spec'}
    tIx = find(strcmp({d.inputs.type}, topType{1}) & ...
        ismember({d.inputs.name}, {d.wideP.source}), 1);
    if isempty(tIx), continue; end
    nm = d.inputs(tIx).name;
    hFig.UserData.activeSrc = nm;
    evt.Key = 'equal'; hFig.WindowKeyPressFcn(hFig, evt);
    assert(hFig.UserData.inputs(tIx).yAdjust > 1, ...
        'integration: Top %s should respond to shift+=', topType{1});
    fprintf('  ok  amplitude on Top %s (%s)\n', topType{1}, nm);
end

% region auto-switch on a Bottom click (no Region dropdown; the Panels header
% follows the click)
d = hFig.UserData; axBot = d.narrowP(1).ax;
axBot.ButtonDownFcn(axBot, []);
assert(strcmp(hFig.UserData.cfgRegion, 'narrow'), 'integration: Bottom click should activate narrow');
assert(contains(hFig.UserData.hPanelsLbl.Text, 'Bottom'), ...
    'integration: Panels header should read Bottom');
fprintf('  ok  region auto-switch\n');

% panel-list edits: Add appends, up/down reorders, X deletes (active = narrow)
nBefore = numel(hFig.UserData.narrowP);
fire(hFig.UserData.hSrcGrid, 'Add');
assert(numel(hFig.UserData.narrowP) == nBefore + 1, 'integration: Add did not append a panel');
first0 = hFig.UserData.narrowP(1).source;
fire(hFig.UserData.hSrcGrid, char(9660));   % move row 1 down
assert(~strcmp(hFig.UserData.narrowP(1).source, first0), 'integration: down did not reorder');
fire(hFig.UserData.hSrcGrid, 'X');          % delete a row
assert(numel(hFig.UserData.narrowP) == nBefore, 'integration: X did not remove a panel');
fprintf('  ok  panel list add / reorder / delete\n');

% two event sets at once; CURATE elects one (or None) without disturbing the
% other. A second set (loaded inline) must not steal curation from 'ripp'.
d = hFig.UserData;
ev2 = struct('peakTime', d.ed.peakTime(1:min(5, d.nEvents)) + 0.01);
d.loadCoreFcn(var_recipe('value', 'data', ev2), 'eventTicks', 'top', ...
    'ripp2', 'ripp2', [], []);
di = hFig.UserData;
assert(sum(strcmp({di.inputs.type}, 'eventTicks')) == 2, ...
    'curate: expected two event sets after Load');
assert(strcmp(di.curate, 'ripp'), 'curate: a 2nd set must not steal the target');
assert(any(strcmp('ripp2', di.hCurateDD.Items)) && any(strcmp('None', di.hCurateDD.Items)), ...
    'curate: selector missing None / ripp2');
fprintf('  ok  two event sets loaded, CURATE lists both + None\n');

% curating a set auto-adds its Bottom lines panel, which is an OVERLAY: it sits
% in narrowP but takes NO tile (its .ax is empty) and draws across the signals.
jb = find(strcmp({di.narrowP.source}, 'ripp'), 1);
assert(~isempty(jb), 'overlay: curating ripp should add a Bottom ripp panel');
assert(~isgraphics(di.narrowP(jb).ax), 'overlay: a lines panel must take no tile');
assert(numel(di.axNarrow) == sum(arrayfun(@(p) isgraphics(p.ax), di.narrowP)), ...
    'overlay: tiled axes must match the non-overlay panels');
fprintf('  ok  Bottom lines panel is an overlay (no tile)\n');

% only Bottom lines-panels mark the window: ripp (curated, auto-added) marks;
% ripp2 does not until it is a Bottom panel too. Curate ripp2 - its lines panel
% is added and ripp's persists, so now BOTH mark: blue (ripp) + orange (ripp2).
di.jumpToTimeFcn(di.ed.peakTime(1));
assert(isempty(findobj(hFig.UserData.axNarrow(1), 'Type', 'line', 'Color', [0.85 0.33 0.10])), ...
    'window marks: ripp2 must not mark until it is a Bottom panel');
di.setCurateFcn('ripp2');
assert(strcmp(hFig.UserData.curate, 'ripp2'), 'curate: switch to ripp2 failed');
dm = hFig.UserData; dm.jumpToTimeFcn(dm.ed.peakTime(1));
axN = hFig.UserData.axNarrow(1);
hasBlue   = ~isempty(findobj(axN, 'Type', 'line', 'Color', [0.00 0.20 0.80]));
hasOrange = ~isempty(findobj(axN, 'Type', 'line', 'Color', [0.85 0.33 0.10]));
assert(hasBlue && hasOrange, 'window marks: both Bottom sets should mark (blue + orange)');
fprintf('  ok  window marks follow the Bottom list (ripp + ripp2)\n');

% Ops render override: switch the CURATED set's Bottom lines to a strip. It now
% takes a tile AND its spanning lines are gone - "ticks" must clear the lines,
% including the current-event emphasis (part of the lines, not an always-on
% cursor). ripp (still a lines panel) keeps marking; ripp2 (strip) does not.
dm = hFig.UserData;
jr = find(strcmp({dm.narrowP.source}, 'ripp2'), 1);
dm.narrowP(jr).render = 'strip';
hFig.UserData = dm; dm.rebuildFcn();
dm = hFig.UserData; dm.jumpToTimeFcn(dm.ed.peakTime(1));
assert(isgraphics(hFig.UserData.narrowP(jr).ax), ...
    'ops: render=strip should give the panel a tile');
assert(isempty(findobj(hFig.UserData.axNarrow(1), 'Type', 'line', 'Color', [0.85 0.33 0.10])), ...
    'ops: ticks must clear the set''s spanning lines (incl. the emphasis)');
assert(~isempty(findobj(hFig.UserData.axNarrow(1), 'Type', 'line', 'Color', [0.00 0.20 0.80])), ...
    'ops: another set still shown as lines must keep marking');
fprintf('  ok  Ops ticks clears the lines (tile + no spanning marks)\n');

% None: view only, target cleared
hFig.UserData.setCurateFcn('');
assert(isempty(hFig.UserData.curate) && strcmp(hFig.UserData.hCurateDD.Value, 'None'), ...
    'curate: None should clear the target');
fprintf('  ok  CURATE None clears the target\n');

% None is a real view mode: tick strips STAY drawn (not blanked by the old
% single-target guard), and Prev/Next step the window itself.
dv = hFig.UserData;
assert(strcmp(dv.mode, 'view'), 'view: None should enter view mode');
axRipp = dv.wideP(strcmp({dv.wideP.source}, 'ripp')).ax;
assert(~isempty(findobj(axRipp, 'Type', 'line')), ...
    'view: the ripp tick strip must stay drawn at None');
t0Before = dv.t0;
evk = struct('Modifier', {{}}, 'Key', 'rightarrow');
hFig.WindowKeyPressFcn(hFig, evk);
dv = hFig.UserData;
assert(abs((dv.t0 - t0Before) - dv.win * 0.8) < 1e-6, ...
    'view: Next should step the window by ~0.8 * win');
fprintf('  ok  view mode: ticks kept, Next steps the window\n');

% PRESET switch accumulates: everything loaded stays available (a preset changes
% the arrangement, not what is loaded). After Ripples -> States, CURATE lists
% ripp AND states AND ripp2; the panel dropdown keeps a Ripples-only signal; and
% the switch does NOT change the target (independent - still None from above).
hFig.UserData.loadPresetFcn('sleep_states');
ds = hFig.UserData;
cur = ds.hCurateDD.Items;
assert(all(cellfun(@(n) any(strcmp(n, cur)), {'ripp', 'ripp2', 'states'})), ...
    'preset switch: CURATE must keep every loaded set (ripp, ripp2, states)');
assert(any(strcmp('rippStack', {ds.inputs.name})), ...
    'preset switch: a Ripples signal must stay available (panel dropdown)');
assert(isempty(ds.curate), 'preset switch: must not change the curate target');
fprintf('  ok  preset switch accumulates sets + signals, target unchanged\n');

close(hFig);
fprintf('== PASSED ==\n');
end

% -------------------------------------------------------------------------
function inp = uinput(type)
% a minimal panel input for calling guiPath_draw in isolation
inp = struct('name', type, 'type', type, 'data', [], 'fs', NaN, 'ylim', [], ...
    'clr', 'k', 'label', '', 'height', 1, 'chInfo', [], 'yAdjust', 1);
end

function fire(grid, txt)
% click the first enabled uibutton in grid whose Text is txt
b = findobj(grid, 'Type', 'uibutton', '-and', 'Text', txt, '-and', 'Enable', 'on');
assert(~isempty(b), 'no enabled button "%s" in the panel list', txt);
b = b(1); b.ButtonPushedFcn(b, []);
end

function s = localSpan(ax)
% max vertical extent of the drawn lines in an axis (ignoring NaN separators)
L = findobj(ax, 'Type', 'line');
lo = inf; hi = -inf;
for k = 1:numel(L)
    y = L(k).YData(~isnan(L(k).YData));
    if ~isempty(y), lo = min(lo, min(y)); hi = max(hi, max(y)); end
end
s = hi - lo;
end
