function h = gui_statePanel(parent, api, names, colors)
% GUI_STATEPANEL  Epoch-state assignment stepper (sibling of gui_eventPanel).
%
%   h = gui_statePanel(parent, api, names, colors) builds a control block for
%   scoring vigilance states epoch by epoch: an epoch-index field (+ "/ N"),
%   Prev/Next, one colored button per state (1..N) plus Undefined (N+1), Save,
%   and a current-state readout. It drops into the pinned action area of
%   gui_layout, swapping in for gui_eventPanel when a viewer curates
%   states rather than events.
%
%   The host owns all state. API is a struct of callbacks invoked on action:
%       api.prev()      - step to the previous epoch
%       api.next()      - step to the next epoch
%       api.assign(k)   - set the current epoch's state to k (N+1 = undefined)
%       api.save()      - persist the labels
%       api.setIdx(v)   - jump to 1-based epoch index v
%   Any missing field is ignored (the control becomes a no-op).
%
%   names/colors are 1xN cell arrays of state names and RGB triplets. The
%   returned struct H carries the control handles and:
%       H.refresh(idx, total, curLabel)
%   which the host calls after every change to sync the index/total and show the
%   current epoch's state (curLabel is the integer label; > N reads "undefined").
%
%   HISTORY:
%       05 Jul 2026 - created for guiPath_curate state-label curation.

narginchk(2, 4)
if ~isstruct(api), error('gui_statePanel:api', 'api must be a struct of callbacks'); end
if nargin < 3 || isempty(names),  names  = {}; end
if nargin < 4 || isempty(colors), colors = {}; end
ns = numel(names);

nBtnRows = max(1, ceil(ns / 2));
nRows    = 2 + nBtnRows + 3;              % idx | prev/next | states | undef+save | status
g = uigridlayout(parent, [nRows, 2], 'RowHeight', repmat({'fit'}, 1, nRows), ...
    'ColumnWidth', {'1x', '1x'}, 'Padding', 4, 'RowSpacing', 4, 'ColumnSpacing', 4);

% row 1: epoch index + total
hIdx = uieditfield(g, 'numeric', 'Limits', [1, Inf], 'RoundFractionalValues', 'on', ...
    'ValueDisplayFormat', '%.0f', ...              % plain integer (not 1.244e+04)
    'Value', 1, 'ValueChangedFcn', @(s, ~) safecall(api, 'setIdx', round(s.Value)));
hIdx.Layout.Row = 1; hIdx.Layout.Column = 1;
hTot = uilabel(g, 'Text', '/ 0', 'VerticalAlignment', 'center');
hTot.Layout.Row = 1; hTot.Layout.Column = 2;

% row 2: prev / next  (arrows match the host keyboard shortcuts)
bPrev = uibutton(g, 'Text', 'Prev (<-)', 'ButtonPushedFcn', @(~, ~) safecall(api, 'prev'));
bPrev.Layout.Row = 2; bPrev.Layout.Column = 1;
bNext = uibutton(g, 'Text', 'Next (->)', 'ButtonPushedFcn', @(~, ~) safecall(api, 'next'));
bNext.Layout.Row = 2; bNext.Layout.Column = 2;

% state buttons 1..ns (number key = state), two per row, colored
hStates = gobjects(1, ns);
for k = 1:ns
    r   = 2 + ceil(k / 2);
    c   = 2 - mod(k, 2);                  % odd k -> col 1, even k -> col 2
    clr = [0.6 0.6 0.6];
    if k <= numel(colors) && numel(colors{k}) >= 3, clr = colors{k}(1:3); end
    b = uibutton(g, 'Text', sprintf('%d %s', k, names{k}), 'BackgroundColor', clr, ...
        'FontColor', txtColor(clr), 'ButtonPushedFcn', @(~, ~) safecall(api, 'assign', k));
    b.Layout.Row = r; b.Layout.Column = c;
    hStates(k) = b;
end

% undefined + save
rU = 2 + nBtnRows + 1;
bUndef = uibutton(g, 'Text', 'Undef (X)', 'ButtonPushedFcn', @(~, ~) safecall(api, 'assign', ns + 1));
bUndef.Layout.Row = rU; bUndef.Layout.Column = 1;
bSave = uibutton(g, 'Text', 'Save (Ctrl+S)', 'ButtonPushedFcn', @(~, ~) safecall(api, 'save'));
bSave.Layout.Row = rU; bSave.Layout.Column = 2;

% status: current epoch's state
hStat = uilabel(g, 'Text', '', 'WordWrap', 'on', 'FontColor', [0.3 0.3 0.3]);
hStat.Layout.Row = nRows; hStat.Layout.Column = [1, 2];

h = struct('grid', g, 'idx', hIdx, 'total', hTot, 'states', hStates, ...
    'undef', bUndef, 'save', bSave, 'status', hStat);
h.names   = {names};                      % wrap so the scalar struct keeps the cell
h.refresh = @(idx, total, curLabel) refreshPanel(h, idx, total, curLabel);

end     % MAIN

% ------------------------------------------------------------------------
function refreshPanel(h, idx, total, curLabel)
% sync the index field + total, and show the current epoch's state name
if ~isvalid(h.idx), return; end
total = max(0, round(total));
h.idx.Limits = [1, max(1, total)];
if total > 0, h.idx.Value = min(max(1, round(idx)), total); end
h.total.Text = sprintf('/ %d', total);
names = h.names{1};
if nargin < 4 || isempty(curLabel)
    h.status.Text = ''; return;
end
if curLabel >= 1 && curLabel <= numel(names)
    h.status.Text = sprintf('state: %s', names{curLabel});
else
    h.status.Text = 'state: undefined';
end
end

% ------------------------------------------------------------------------
function c = txtColor(rgb)
% white text on dark buttons, black on light (luminance split)
if numel(rgb) < 3, c = 'k'; return; end
if (0.299 * rgb(1) + 0.587 * rgb(2) + 0.114 * rgb(3)) < 0.55, c = 'w'; else, c = 'k'; end
end

% ------------------------------------------------------------------------
function safecall(api, name, varargin)
% invoke api.(name) if present and non-empty; otherwise no-op
if isfield(api, name) && ~isempty(api.(name))
    api.(name)(varargin{:});
end
end

% EOF
