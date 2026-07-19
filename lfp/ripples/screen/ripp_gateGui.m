function hFig = ripp_gateGui(basepath, varargin)
% RIPP_GATEGUI Interactively tune the ripple QA gate: kept vs removed waveforms.
%
%   hFig = RIPP_GATEGUI(basepath, varargin)
%
%   SUMMARY:
%       A bench for the false-positive gate. Takes the all-events output of
%       ripp_detect (every detected ripple, kept and rejected, with its LFP map
%       and QA metrics) and plots the mean waveform of the events a threshold
%       KEEPS against those it REMOVES - live, as the MUA-gain / EMG / prominence
%       thresholds are moved. A good gate makes the kept mean a clean oscillation
%       and the removed mean a large non-oscillatory deflection; the running count
%       shows how many events each threshold keeps. Built on guiTbl_xy (grouped
%       mean +/- spread over the shared graphics/gui layer). Re-detection is
%       optional and never re-implemented: pass the in-memory [ripp, aux] from a
%       prior ripp_detect to skip it, else it calls ripp_detect once. The saved
%       maps hold accepted events only, so a kept-vs-removed view needs the
%       all-events maps that live in aux.rippMaps - not the .rippMaps.mat on disk.
%       Restricted to NREM (the regime the gate targets). Writes nothing.
%
%   INPUTS:
%       basepath - <char> session directory (used only if ripp/aux are absent).
%       varargin - Parameter/Value:
%           'ripp'    - <struct> all-events ripp from ripp_detect. {detect here}
%           'aux'     - <struct> aux from ripp_detect (needs .rippMaps,
%                                .nremTimes). {detect here}
%           'met'     - <struct> ripp_methods config; sets the initial thresholds
%                                and the detect fallback. {ripp_methods('default')}
%           'win'     - <vec>    window for the detect fallback (s). {[0 3*3600]}
%           'Visible' - <char>   'on' | 'off' for headless use. {'on'}
%
%   OUTPUT:
%       hFig - <handle> the GUI figure.
%
%   DEPENDENCIES:
%       ripp_detect, ripp_methods, guiTbl_xy.
%
%   HISTORY:
%       260719 replaces the static ripp_gateFig with an interactive gate on
%              guiTbl_xy; drops the re-detection chain and the rho scatter.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'ripp', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'aux', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'met', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'win', [0 3 * 3600], @isnumeric);
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(char(x), {'on', 'off'})));
parse(p, basepath, varargin{:});
ripp = p.Results.ripp;
aux  = p.Results.aux;
met  = p.Results.met;
win  = p.Results.win;
vis  = char(p.Results.Visible);

if isempty(met)
    met = ripp_methods('default');
end

%% ========================================================================
%  EVENTS (from a prior detection, or detect once)
%  ========================================================================

if isempty(ripp) || isempty(aux)
    [ripp, aux] = ripp_detect(basepath, 'met', met, 'win', win, ...
        'mapDur', [-0.06 0.06], 'verbose', true);
end
maps = aux.rippMaps;

% NREM events only (the regime the gate targets)
sel = inAny(ripp.peakTime, aux.nremTimes);
if ~any(sel)
    warning('ripp_gateGui:noNrem', 'No NREM events; showing all events.');
    sel = true(numel(ripp.peakTime), 1);
end

st = struct();
st.xt      = maps.tstamps * 1000;               % ms
st.lfpMap  = maps.lfp(sel, :);
st.filtMap = [];
if isfield(maps, 'filt')
    st.filtMap = maps.filt(sel, :);
end
st.gain = ripp.spkGain(sel);
st.emg  = ripp.emg(sel);
st.prom = ripp.peakProm(sel);
st.nEvt = nnz(sel);

%% ========================================================================
%  LAYOUT: threshold bar on top, guiTbl_xy below
%  ========================================================================

hFig = uifigure('Name', 'Ripple QA gate', 'Position', [80 80 1500 850], ...
    'Visible', vis);
gMain = uigridlayout(hFig, [2 1], 'RowHeight', {'fit', '1x'}, ...
    'Padding', 6, 'RowSpacing', 6);

gBar = uigridlayout(gMain, [1 8], 'Padding', 2, 'ColumnSpacing', 6, ...
    'ColumnWidth', {90, 65, 60, 65, 60, 65, 260, '1x'});
gBar.Layout.Row = 1;

% edit fields wired to the redraw at creation (spares a post-assignment that
% would look like an unused struct write)
cbk = @(~, ~) refresh(hFig);
uilabel(gBar, 'Text', 'MUA gain >=', 'HorizontalAlignment', 'right');
st.edGain = uieditfield(gBar, 'numeric', 'Value', met.gainThr, ...
    'ValueChangedFcn', cbk);
uilabel(gBar, 'Text', 'EMG <=', 'HorizontalAlignment', 'right');
st.edEmg = uieditfield(gBar, 'numeric', 'Value', met.thrEmg, ...
    'ValueChangedFcn', cbk);
uilabel(gBar, 'Text', 'prom >=', 'HorizontalAlignment', 'right');
st.edProm = uieditfield(gBar, 'numeric', 'Value', -Inf, ...
    'ValueChangedFcn', cbk);
st.lblCount = uilabel(gBar, 'Text', '', 'FontWeight', 'bold');

st.hPanel = uipanel(gMain, 'BorderType', 'none');
st.hPanel.Layout.Row = 2;

hFig.UserData = st;
refresh(hFig);

end     % EOF


% =========================================================================
%  LOCALS
% =========================================================================
function refresh(hFig)
% Recompute the kept/removed split from the current thresholds and redraw. The
% widget is rebuilt (not mutated) because guiTbl_xy owns its own data copy; the
% rebuild is a sub-second tiledlayout pass and only fires on threshold commit.
st = hFig.UserData;

% a NaN metric means the criterion is unavailable, so it does not reject
% (mirrors evt_qa): the event is kept on that axis
keep = (isnan(st.gain) | st.gain >= st.edGain.Value) & ...
       (isnan(st.emg)  | st.emg  <= st.edEmg.Value)  & ...
       (isnan(st.prom) | st.prom >= st.edProm.Value);

status = repmat("removed", numel(keep), 1);
status(keep) = "kept";
status = categorical(status, {'removed', 'kept'});

tbl = table(st.lfpMap, 'VariableNames', {'lfp'});
if ~isempty(st.filtMap)
    tbl.filt = st.filtMap;
end
tbl.status = status;

% preserve the Y-variable choice across the rebuild
curY = 'lfp';
ud = st.hPanel.UserData;
if isstruct(ud) && isfield(ud, 'yVar')
    curY = ud.yVar;
end

delete(allchild(st.hPanel));
st.hPanel.UserData = [];
guiTbl_xy(st.xt, tbl, 'Parent', st.hPanel, 'yVar', curY, ...
    'grpVar', 'status', 'xLbl', 'time (ms)');

st.lblCount.Text = sprintf('kept %d  |  removed %d  (of %d NREM)', ...
    nnz(keep), nnz(~keep), st.nEvt);
hFig.UserData = st;

end     % refresh


function tf = inAny(t, wins)
% true for each t inside any [start end] row of wins
tf = false(numel(t), 1);
if isempty(wins)
    return;
end
for iWin = 1:size(wins, 1)
    tf = tf | (t >= wins(iWin, 1) & t <= wins(iWin, 2));
end
end     % inAny
