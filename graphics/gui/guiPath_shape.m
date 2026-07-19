function inp = guiPath_shape(inp, nstates)

% Shape a raw varMap value into the drawn input, dispatched on the panel type.
%
% The loader (var_load) is view-blind: it returns raw .data / .fs. This is the
% view half - it turns that raw value into exactly what guiPath_draw expects:
% event / state structs, an hours-cell hypnogram, a channel stack's display
% stats, a resolved y-limit. Called once per input when a preset lands, not on
% every redraw. A stateStrip arrives already composed (by the states preset),
% so it passes through.
%
% INPUTS
% - inp             <struct> a partial input record with .type, raw .data, .fs,
%                   .labels (bin channel numbers) and .ylim (a spec, for a trace).
% - nstates         <num>(opt) state count, to pad a hypnogram's cell.
%
% OUTPUTS
% - inp             <struct> .data shaped, .fs adjusted, .chInfo set (traces),
%                   .ylim resolved (trace).
%
% DEPENDENCIES
% - none. See guiPath_draw for what consumes the shaped input.
%
% HISTORY
% - 260719          created (per-type shaping, moved verbatim from guiPath_load
%                   when var_load became the single, view-blind loader).

if nargin < 2, nstates = []; end

switch inp.type
    case 'trace'
        inp.data = double(inp.data(:));
        inp.ylim = resolveYlim(inp.ylim, inp.data);

    case 'traces'
        inp.chInfo = traceStack(inp.data, inp.labels);   % native class kept

    case 'spec'
        inp.data = specAdapter(inp.data);
        inp.fs   = NaN;

    case 'hypnogram'
        inp.data = toHoursCell(inp.data, nstates);
        inp.fs   = NaN;

    case 'raster'
        spk = inp.data;
        if ~iscell(spk), spk = {spk(:)}; end
        inp.data = spk(:);
        inp.fs   = NaN;

    case 'eventTicks'
        inp.data = eventsFrom(inp.data);
        inp.fs   = NaN;

    case 'stateStrip'
        % already composed by the states preset (an inline strip struct)

    otherwise
        error('guiPath_shape:type', 'unknown panel type "%s"', inp.type);
end

end


% =========================================================================
%  SHAPERS (pure; moved from guiPath_load)
% =========================================================================

function yl = resolveYlim(spc, data)
% Resolve a panel's ylim spec against its data:
%       [lo hi]     absolute limits, taken as given
%       <scalar> p  percentile clip to [p, 100-p]; 0 <= p < 50, 0 = full range
%       'prc'       percentile clip at prcDflt
%       'full', []  autoscale
%
% The percentile is estimated on a subsample (~1e5 points) rather than the whole
% signal: on a full-session trace the clip is visually identical but the sort is
% orders of magnitude cheaper. A degenerate range (flat / NaN, or p >= 50) is
% dropped so the axis autoscales.
prcDflt = 0.1;

if isnumeric(spc) && numel(spc) == 2, yl = spc; return; end
yl = [];

prc = [];
if ischar(spc) && strcmp(spc, 'prc'),  prc = prcDflt; end
if isnumeric(spc) && isscalar(spc),    prc = spc;     end
if isempty(prc) || isempty(data), return; end

n = numel(data);
if n > 2e5, s = double(data(1:ceil(n / 1e5):end)); else, s = double(data(:)); end
yl = prctile(s, [prc, 100 - prc]);
if numel(yl) ~= 2 || ~all(isfinite(yl)) || yl(2) <= yl(1), yl = []; end
end


function info = traceStack(val, labels)
% fixed display stats for a channel stack, from a subsample of the whole signal:
%   .spacing  vertical gap between channels (data units)
%   .base     per-channel baseline (median), subtracted so each channel is
%             centred on its own row regardless of DC offset
%   .labels   channel numbers for the y-axis
% Computed once and held for the panel's life so the stack neither drifts nor
% rescales as the window moves. Spacing is robust (MAD-based) so one large
% channel does not blow the stack apart.
nCh = size(val, 2);
nr  = size(val, 1);
if nr > 2e5, s = double(val(1:ceil(nr / 1e5):end, :)); else, s = double(val); end
base = median(s, 1);
sd = median(abs(s - base), 1) / 0.6745;            % per-channel robust SD
sp = 6 * median(sd(isfinite(sd)));
if isempty(sp) || ~isfinite(sp) || sp <= 0, sp = 1; end
if ~isempty(labels) && numel(labels) == nCh, labs = labels(:)'; else, labs = 1:nCh; end
info = struct('spacing', sp, 'base', base(:)', 'labels', labs);
end


function adapter = specAdapter(raw)
% the {s, freq, tstamps} plot_spec adapter from the raw sleep_sig spec fields
adapter = struct('s', raw.spec, 'freq', raw.spec_freq(:), ...
    'tstamps', raw.spec_tstamps(:));
end


function boutHr = toHoursCell(bt, nstates)
% sleep-state bouts (seconds cell) -> hours cell of length nstates
if ~iscell(bt), bt = {bt}; end
if isempty(nstates), nstates = numel(bt); end
boutHr = cell(1, nstates);
for iState = 1 : nstates
    if iState <= numel(bt) && ~isempty(bt{iState})
        boutHr{iState} = bt{iState} / 3600;
    else
        boutHr{iState} = zeros(0, 2);
    end
end
end


function ev = eventsFrom(val)
% compact events struct from a detection struct (ed / ripp) or an Nx1/2/3 matrix
if isstruct(val)
    if ~isfield(val, 'peakTime') || isempty(val.peakTime)
        error('guiPath_shape:events', 'an events struct needs a non-empty peakTime');
    end
    ev = struct('peakTime', val.peakTime(:));
    if isfield(val, 'times') && ~isempty(val.times),  ev.times = val.times; end
    if isfield(val, 'accepted') && ~isempty(val.accepted)
        ev.accepted = logical(val.accepted(:));
    end
    if isfield(val, 'state') && ~isempty(val.state),  ev.state = val.state(:); end
    return
end
M = double(val);
if ~isnumeric(M) || isempty(M)
    error('guiPath_shape:events', 'events must be a struct or numeric matrix');
end
if isvector(M)
    ev = struct('peakTime', M(:));
elseif size(M, 2) == 3
    ev = struct('peakTime', M(:, 2), 'times', [M(:, 1), M(:, 3)]);
elseif size(M, 2) == 2
    ev = struct('peakTime', mean(M, 2), 'times', M);
else
    error('guiPath_shape:events', 'events matrix must be Nx1, Nx2, or Nx3');
end
end

% EOF
