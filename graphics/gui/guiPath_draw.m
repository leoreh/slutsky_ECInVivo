function guiPath_draw(ax, inp, data, a, b, xf)

% Draws one panel into an axis, dispatched on the panel input's type.
%
% Every draw is pure: it renders into the axis it is handed and touches no
% shared state. The host (guiPath) owns the axis, clears it, and sets the
% x-limits and label around this call; adding a panel type means adding a case
% here and a typeDefaults entry in guiPath_panel, nothing else.
%
% Time is handed over in two units at once, which is the one subtlety here.
% The range [a, b] and any signal's own time base are in SECONDS, while the
% axis is drawn in DISPLAY units (hours for the Top, seconds for the Bottom).
% xf converts: display = seconds / xf. So a draw slices its data with a and b,
% then divides the x it plots by xf.
%
% INPUTS
% - ax              <handle> target axis. Already cleared and held by the host.
% - inp             <struct> the panel's input record: .type .data .fs .clr
%                   .ylim. See guiPath_load.
% - data            <struct> the host's live state. Only the curation types
%                   read it: eventTicks / stateStrip compare inp.name to
%                   data.curate to tell the edited set (the live accept mask /
%                   labels) from a read-only one (drawn from its own loaded
%                   data). The edited event set draws in its own colour too.
%                   traces instead reads .chInfo (stack spacing / baseline /
%                   labels) off inp; trace / traces / spec read inp.yAdjust, the
%                   live amplitude factor (shift+scroll, shift+/-/0), applied as
%                   y-limits / per-channel gain / brightness respectively.
% - a, b            <num> the window to draw, in seconds.
% - xf              <num> seconds per display unit (3600 for hours, 1 for s).
%
% SEE ALSO
% - guiPath
% - guiPath_panel
% - guiPath_load
% - guiPath_doc
%
% HISTORY
% - 260716          extracted from guiPath's PANEL DRAW block, unchanged, so
%                   the entry point stops growing and a new panel type has a
%                   file of its own to land in.

switch inp.type
    case 'trace',      drawTrace(ax, inp, a, b, xf);
    case 'traces',     drawTraces(ax, inp, a, b, xf);
    case 'spec',       drawSpec(ax, inp, xf);
    case 'hypnogram',  drawHypno(ax, inp, xf);
    case 'eventTicks', drawTicks(ax, inp, data, xf);
    case 'stateStrip', drawStateStrip(ax, inp, data, xf);
    case 'raster',     drawRaster(ax, inp, a, b, xf);
end

end     % MAIN

%% ========================================================================
%  SIGNAL PANELS
%  ========================================================================

function drawTrace(ax, inp, a, b, xf)
% windowed raw trace, decimated for display, x in display units
sig = inp.data; fs = inp.fs;
s1 = max(1, floor(a * fs) + 1);
s2 = min(numel(sig), ceil(b * fs) + 1);
if s2 < s1, return; end          % window outside this signal's extent -> nothing to draw
rng = s1:s2;
t = ((rng - 1) / fs) / xf;
np = numel(rng); maxPts = 20000;
if np > maxPts
    st = ceil(np / maxPts);
    plot(ax, t(1:st:end), sig(rng(1:st:end)), 'Color', inp.clr);
else
    plot(ax, t, sig(rng), 'Color', inp.clr);
end
% only a valid, increasing, finite range (a flat / NaN signal gives lo==hi)
if numel(inp.ylim) == 2 && all(isfinite(inp.ylim)) && inp.ylim(2) > inp.ylim(1)
    ax.YLim = ampScale(inp.ylim, inp);
end
end

function drawTraces(ax, inp, a, b, xf)
% a vertical stack of channels: channel 1 on top, each centred on its own
% baseline and offset down by chInfo.spacing. Native (int16) whole-session
% data; only the window slice is cut and cast, so navigation is a pure index.
%
% The amplitude factor is a GAIN on each channel's deflection - spacing and
% y-limits stay fixed, so a bigger factor makes the wiggles bigger IN PLACE
% (they may cross into a neighbour; they never clip a shrinking axis). This is
% deliberately unlike a single trace, where amplitude tightens the y-limits:
% for a stack that would just push the outer channels off-panel.
sig = inp.data; fs = inp.fs;
nCh = size(sig, 2);
s1 = max(1, floor(a * fs) + 1);
s2 = min(size(sig, 1), ceil(b * fs) + 1);
if s2 < s1 || nCh == 0, return; end
rng = s1:s2;
np = numel(rng); maxPts = 20000;
if np > maxPts, rng = rng(1:ceil(np / maxPts):end); end
t = ((rng - 1) / fs) / xf;

ci = inp.chInfo;
sp = ci.spacing; g = ampFactor(inp);
for iCh = 1:nCh
    y = (double(sig(rng, iCh)) - ci.base(iCh)) * g - (iCh - 1) * sp;
    plot(ax, t, y, 'Color', inp.clr);
end
ax.YLim = [-(nCh - 1) * sp - sp, sp];              % fixed; gain grows the wiggles
ax.YTick = fliplr(-((1:nCh) - 1) * sp);            % ascending for the axis
ax.YTickLabel = string(fliplr(ci.labels));
end

function ya = ampFactor(inp)
% the panel's live amplitude factor (1 = as loaded), set by shift+scroll and
% shift+/-/0. Each draw applies it in the axis natural to its type: a trace
% tightens its y-limits (ampScale), a stack scales its per-channel gain
% (drawTraces), a spectrogram tightens its colour range (drawSpec).
ya = 1;
if isfield(inp, 'yAdjust') && ~isempty(inp.yAdjust) && inp.yAdjust > 0
    ya = inp.yAdjust;
end
end

function yl = ampScale(yl, inp)
% tighten / loosen y-limits about their centre by the amplitude factor: > 1
% makes the signal fill more of the panel. Multiplies the loaded limits rather
% than replacing them, so reset (factor 1) restores exactly what the preset
% resolved, with no stored original. (For a single trace only; a stack scales
% gain instead - see drawTraces.)
ya = ampFactor(inp);
if ya == 1, return; end
c = mean(yl);
yl = c + (yl - c) / ya;
end

function drawSpec(ax, inp, xf)
plot_spec(inp.data, 'axh', ax, 'saveFig', false, 'xtime', xf);
% same amplitude gesture, applied to brightness: tighten the colour range from
% the top (bigger factor = brighter). plot_spec leaves CLim on auto, so CLim
% here is the data range.
ya = ampFactor(inp);
if ya ~= 1 && numel(ax.CLim) == 2 && diff(ax.CLim) > 0
    cl = ax.CLim;
    ax.CLim = [cl(1), cl(1) + diff(cl) / ya];
end
end

function drawHypno(ax, inp, xf)
% bout times arrive in hours; convert to display units. Pin sstates to the
% number of bout-cells provided so the strip is independent of cfg.nstates.
bt = cellfun(@(x) x * 3600 / xf, inp.data, 'uni', false);
plot_hypnogram('boutTimes', bt, 'sstates', 1:numel(bt), 'style', 'strip', 'hAx', ax);
end

function drawRaster(ax, inp, a, b, xf)
% inward ticks so the unit numbers stay but no tick marks protrude left
spk = cellfun(@(s) s(s >= a & s <= b) / xf, inp.data, 'uni', false);
if all(cellfun(@isempty, spk)), return; end   % skip (avoids plot_raster's empty warning)
plot_raster(spk, 'hAx', ax, 'xLim', [a, b] / xf, 'flgLbls', false, 'tickDir', 'in');
end

%% ========================================================================
%  CURATION PANELS (these read the host's live state)
%  ========================================================================

function drawTicks(ax, inp, data, xf)
% full-height ticks, one per event, in the set's own colour so several event
% sets read apart (the overview mirror of the Bottom-window marks). The curated
% set uses the live accept mask, so rejecting removes its tick; any other set
% uses its own stored mask if it has one, else shows every event. No y-ticks (a
% horizontal y-label is the only label); the x-axis stays, so the strip shows
% the Time axis when it is the region's bottom panel.
peaks = [];
if isstruct(inp.data) && isfield(inp.data, 'peakTime'), peaks = inp.data.peakTime(:); end
if isfield(data, 'curate') && strcmp(inp.name, data.curate)
    if numel(data.accepted) == numel(peaks), peaks = peaks(data.accepted); end
elseif isstruct(inp.data) && isfield(inp.data, 'accepted') ...
        && numel(inp.data.accepted) == numel(peaks)
    peaks = peaks(logical(inp.data.accepted));
end
drawTickLine(ax, peaks / xf, inp.clr);
ylim(ax, [0, 1]); ax.YTick = [];   % no y-ticks; the x-axis stays (Time, when last)
end

function drawTickLine(ax, x, clr)
% one full-height vertical line per x
if isempty(x), return; end
x = x(:)';
X = [x; x; nan(1, numel(x))];
Y = repmat([0; 1; NaN], 1, numel(x));
line(ax, X(:), Y(:), 'Color', clr);
end

function drawStateStrip(ax, inp, data, xf)
% per-epoch coloured label strip drawn as one truecolor image (fast to redraw).
% the curated strip reads the live labels + epoch centres; any other state set
% is read-only, drawn from its own loaded data. undefined (> nstates) render
% gray.
if isfield(data, 'curate') && strcmp(inp.name, data.curate)
    T = data.ed.peakTime(:)'; L = data.labels(:)';
    ns = data.nstates; colors = data.stateColors;
else
    D = inp.data;
    if ~isstruct(D) || ~isfield(D, 'labels') || isempty(D.labels), return; end
    T = D.epochT(:)'; L = D.labels(:)';
    if isfield(D, 'nstates') && ~isempty(D.nstates)
        ns = D.nstates;
    elseif isfield(D, 'names') && ~isempty(D.names)
        ns = numel(D.names);
    else
        ns = max(1, max(L));
    end
    if isfield(D, 'colors'), colors = D.colors; else, colors = {}; end
end
n = min(numel(T), numel(L));
if n == 0, return; end
T = T(1:n); L = L(1:n);
cmap = stateCmap(ns, colors);
Lc = min(max(round(L), 1), size(cmap, 1));
cdata = reshape(cmap(Lc, :), [1, n, 3]);
if n == 1, xl = [T(1) - 0.5, T(1) + 0.5]; else, xl = [T(1), T(end)]; end
image(ax, 'XData', xl / xf, 'YData', [0, 1], 'CData', cdata);
ax.YLim = [0, 1]; ax.YTick = [];   % no y-ticks; the x-axis stays (Time, when last)
end

function cmap = stateCmap(ns, colors)
% (nstates+1) x 3 state colormap; the extra row (undefined) is gray
ns = max(1, ns);
cmap = repmat([0.6 0.6 0.6], ns + 1, 1);
for iState = 1:min(ns, numel(colors))
    c = colors{iState};
    if numel(c) >= 3, cmap(iState, :) = c(1:3); end
end
end

% EOF
